// Copyright (C) 2013 Lile Wang
// For full notice, see "main.cpp" and "COPYING".

#include "driver.h"

////////////////////////////////////////////////////////////
// Constructor and destructor

driver::driver( const std::string & par_file )
{
    start_t = omp_get_wtime(  );
    this->par_file_name = par_file;
}

driver::~driver(  ){  }

////////////////////////////////////////////////////////////
// Read from parameter file

void driver::read_from_par(  )
{
    input read_par( par_file_name );
    read_par.read(  );
    read_par.find_key( "corr_stat", corr_stat, 1 );
    read_par.find_key( "num_threads", num_threads, 1 );
    read_par.find_key( "bin_count_type", bin_count_type, 1 );
    read_par.find_key( "s_max", s_max, 30. );
    read_par.find_key( "s_min", s_min, 1. );
    read_par.find_key( "s_bin_num", s_num, 20 );
    read_par.find_key( "phi_bin_num", phi_num, 40 );
    read_par.find_key( "log_bin", log_bin, 0 );
    read_par.find_key( "regular_phi_bin", regular_phi_bin, 0 );
    read_par.find_key( "weighted_bin", weighted_bin, 0 );    
    read_par.find_key( "file_data", data_file_name, "catalog" );
    read_par.find_key( "file_rand", rand_file_name, "random" );
    read_par.find_key( "out_name_base", out_name_base, "" );
    read_par.find_key( "lambda", lambda, 0.7 );
    read_par.find_key( "z_max", z_max, 5. );    
    read_par.find_key( "jackknife_depth", jk_depth, 4 );
    jk_depth = jk_depth <= 0 ? 4 : jk_depth;
    jk_num = ( 1 << jk_depth );
    
    return;
}

////////////////////////////////////////////////////////////
// Read from the data file and build trees

void driver::build_trees(  )
{
    read_data read;
    read.set_cosmology( lambda, z_max );
    read.set_par( corr_stat == 0, weighted_bin == 1 );
    
    std::cout << "Reading data from files... " << std::flush;
    read.read_from_file( data_file_name, data );
    read.read_from_file( rand_file_name, rand );
    
    std::cout << "Done.\n" << "Building trees..." << std::flush;
    kdtree::set_jackknife_depth( jk_depth );
    data.build_tree(  );
    rand.build_tree(  );
    std::cout << " Done.\n";
    show_wall_t( "Preprocessing", start_t, precomp_t );

    return;
}

////////////////////////////////////////////////////////////
// Conduct calculations

void driver::cal(  )
{
    const bool get_ls = ( bin_count_type == 1 );
    const bool get_all_bins = ( bin_count_type < 2 );
    const bool get_dd = ( bin_count_type == 2 )
        || ( bin_count_type == 5 ) || ( bin_count_type == 6 );
    const bool get_rr = ( bin_count_type == 3 )
        || ( bin_count_type == 6 ) || ( bin_count_type == 7 );
    const bool get_dr = ( bin_count_type == 4 )
        || ( bin_count_type == 5 ) || ( bin_count_type == 7 );

    std::string dd_name( out_name_base );
    std::string rr_name( out_name_base );
    std::string dr_name( out_name_base );
    if( out_name_base.size(  ) < 1 )
    {
        dd_name = data_file_name;
        rr_name = rand_file_name;
        dr_name = data_file_name + "_" + rand_file_name;
    }

    parallel para_corr;
    para_corr.set_num_threads( num_threads );
    correlate::set_par( s_max, s_min, s_num, log_bin > 0,
			phi_num, regular_phi_bin > 0,
			corr_stat, jk_num );

    if( get_all_bins || get_dd )
    {
        para_corr.cal_corr( data, data );
        dd    = correlate::bin_count_ref(  );
        dd_jk = correlate::bin_count_jk_ref(  );
        correlate::output( dd_name + "_ddbins" );
    }
    if( get_all_bins || get_rr )
    {
        para_corr.cal_corr( rand, rand );
        rr    = correlate::bin_count_ref(  );
        rr_jk = correlate::bin_count_jk_ref(  );
        correlate::output( rr_name + "_rrbins" );
    }
    show_wall_t( "Auto-correlation", precomp_t, auto_t );
    if( get_all_bins || get_dr )
    {
        para_corr.cal_corr( data, rand );
        dr    = correlate::bin_count_ref(  );
        dr_jk = correlate::bin_count_jk_ref(  );
        correlate::output( dr_name + "_drbins" );
    }
    show_wall_t( "Cross-correlation", auto_t, cross_t );

    if( get_ls )
        cal_ls(  );
    return;
}

void driver::cal_ls(  )
{
    std::vector<double> cor, ecor;
    cor.resize( dd.size(  ), 0. );
    ecor.resize( dd.size(  ), 0. );
    for( unsigned i = 0; i < dd.size(  ); ++ i )
    {
        rr[ i ] == 0. ? rr[ i ] = 1 : rr[ i ] = rr[ i ];
        // Note that dd and rr are not multiplied by 2.
        cor[ i ] = 1. + ( dd[ i ] - dr[ i ] ) / rr[ i ];
        double ecor_local_sum( 0. );
        for( int j = 0; j < jk_num; ++ j )
        {
            const int k = j + i * jk_num;
            double ecor_local = dd[ i ] - dd_jk[ k ];
            ecor_local -= dr[ i ] - dr_jk[ k ];
            ecor_local /= rr[ i ] - rr_jk[ k ];
            ecor_local += 1.;
	    if( std::isinf( ecor_local ) ||
		std::isnan( ecor_local ) )
		continue;
            ecor_local_sum += pow( ecor_local - cor[ i ], 2 )
                * jk_num / ( jk_num - 1. );
        }
        ecor[ i ] = sqrt( ecor_local_sum );
    }

    std::string out_file_name = out_name_base + "_corr";
    if( out_name_base.size(  ) < 1 )
        out_file_name = data_file_name + "_"
            + rand_file_name + "_corr";
    
    std::ofstream fout( out_file_name.c_str(  ) );
    for( int i = 0; i < s_num; ++ i )
    {
        typedef correlate corr;
        const double s_bin_min = corr::s_val( i, 0. );
        const double s_bin_cen = corr::s_val( i );
        const double s_bin_max = corr::s_val( i, 1. );
        if( corr_stat == 2 )
            for( int j = 0; j < phi_num; ++ j )
            {
                const double phi_bin_min = corr::phi_val( j, 0. );
                const double phi_bin_cen = corr::phi_val( j );
                const double phi_bin_max = corr::phi_val( j, 1. );
                const int k = phi_num * i + j;
                fout << s_bin_min << '\t' << s_bin_cen
                     << '\t' << s_bin_max << '\t'
                     << phi_bin_min << '\t' << phi_bin_cen
                     << '\t' << phi_bin_max << '\t'
                     << cor[ k ] << '\t' << ecor[ k ] << '\n';
            }
        else
            fout << s_bin_min << '\t' << s_bin_cen << '\t'
                 << s_bin_max << '\t' << cor[ i ] << '\t'
                 << ecor[ i ] << '\n';
    }
    return;
}

////////////////////////////////////////////////////////////
// Timer; showing the wall time

void driver::show_wall_t( std::string title,
                          double & start, double & end )
{
    end = omp_get_wtime(  );
    std::cout << std::setw( 25 ) << std::left << title
              << ": " << std::setw( 8 ) << std::left 
              << end - start << " sec." << std::endl;
    return;
}

////////////////////////////////////////////////////////////
// Get everything!

void driver::go(  )
{
    read_from_par(  );
    build_trees(  );
    cal(  );
    show_wall_t( "Total time", start_t, cross_t );
    return;
}

