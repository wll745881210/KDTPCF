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
	double start( 0 ), precomp( 0 );
	double autocal( 0 ), crosscal( 0 );
	start = omp_get_wtime(  );
	
	input read_par( par_file_name );
	read_par.read(  );
	
	kdtree data;
	kdtree rand;
	parallel para_corr;
	read_data read;

	double lambda( 0. ), z_max( 0. );
	read_par.find_key( "lambda", lambda );
	read_par.find_key( "z_max", z_max );
	read.set_cosmology( lambda, z_max );
	
	std::string data_file_name, rand_file_name;
	read_par.find_key( "file_data", data_file_name );
	read_par.find_key( "file_rand", rand_file_name );	
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
	precomp = omp_get_wtime(  );
	std::cout << "Done." << std::endl;
	
	int num_threads( 0 );
	read_par.find_key( "num_threads", num_threads );
	para_corr.set_num_threads( num_threads );
	double s_max( 0. ), s_min( 0. );
	int s_bin( 0 );
	read_par.find_key( "s_max", s_max );
	read_par.find_key( "s_min", s_min );
	read_par.find_key( "s_bin", s_bin );
	para_corr.set_dist_bin( s_max, s_min, s_bin );
	
	para_corr.cal_corr( data, data );
	para_corr.output( data_file_name + "_ddbins" );
	para_corr.cal_corr( rand, rand );
	para_corr.output( rand_file_name + "_rrbins" );
	autocal = omp_get_wtime(  );
	para_corr.cal_corr( data, rand );
	para_corr.output( data_file_name + "_" + 
				 rand_file_name + "_drbins" );
	crosscal = omp_get_wtime(  );

	show_wall_t( "Total time", start, crosscal );
	show_wall_t( "Preprocessing", start, precomp );
	show_wall_t( "Auto-correlation", precomp, autocal );
	show_wall_t( "Cross-correlation", autocal, crosscal );
	return;
}

