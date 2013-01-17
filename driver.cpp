// Copyright (C) 2013 Lile Wang
// For full notice, see "main.cpp" and "COPYING".

#include "driver.h"
#include "kdtree.h"
#include "correlate.h"
#include "read_data.h"
#include "input.h"
#include "parallel.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <ctime>
#include <omp.h>
#include <cmath>

void show_wall_t( std::string title, double start,
                  double finish )
{
    std::cout << std::setw( 25 ) << std::left << title+": "
              << std::setw( 8 ) << std::left 
              << finish - start << " sec."<< std::endl;
    return;
}

void driver( const std::string & par_file_name )
{
    const double start_t = omp_get_wtime(  );
    
    input read_par( par_file_name );
    kdtree data, rand;
    parallel para_corr;
    read_data read;

    int corr_stat( 0 ), num_threads( 0 );
    double s_max( 0. ), s_min( 0. );
    int s_bin( 0 ), phi_bin( 0 ), log_bin( 0 );
    std::string data_file_name, rand_file_name;
    double lambda( 0. ), z_max( 0. );
    read_par.read(  );
    read_par.find_key( "corr_stat", corr_stat );
    read_par.find_key( "num_threads", num_threads );
    read_par.find_key( "s_max", s_max );
    read_par.find_key( "s_min", s_min );
    read_par.find_key( "s_bin", s_bin );
    read_par.find_key( "phi_bin", phi_bin );
    read_par.find_key( "log_bin", log_bin );    
    read_par.find_key( "file_data", data_file_name );
    read_par.find_key( "file_rand", rand_file_name );
    read_par.find_key( "lambda", lambda );
    read_par.find_key( "z_max", z_max );
    
    read.set_cosmology( lambda, z_max );
    read.set_ang_cor( corr_stat == 0 );
    std::cout << "Reading data from files..." << std::flush;
    read.read_from_file( data_file_name, data );
    read.read_from_file( rand_file_name, rand );
    std::cout << "Done." << std::endl;
    
    std::cout << "Building trees..." << std::flush;
    omp_set_num_threads( 2 );
    #pragma omp parallel sections
    {
        #pragma omp section
        data.build_tree(  );
        #pragma omp section
        rand.build_tree(  );
    }
    const double precomp_t = omp_get_wtime(  );
    std::cout << "Done." << std::endl;
    
    para_corr.set_num_threads( num_threads );
    correlate::set_par( s_max, s_min, s_bin, phi_bin,
                        log_bin > 0, corr_stat );
    
    para_corr.cal_corr( data, data );
    correlate::output( data_file_name + "_ddbins" );
    para_corr.cal_corr( rand, rand );
    correlate::output( rand_file_name + "_rrbins" );
    const double autocal_t = omp_get_wtime(  );
    para_corr.cal_corr( data, rand );
    correlate::output( data_file_name + "_" + 
                       rand_file_name + "_drbins" );
    const double crosscal_t = omp_get_wtime(  );

    show_wall_t( "Total time", start_t, crosscal_t );
    show_wall_t( "Preprocessing", start_t, precomp_t );
    show_wall_t( "Auto-correlation", precomp_t, autocal_t );
    show_wall_t( "Cross-correlation", autocal_t, crosscal_t );
    return;
}

