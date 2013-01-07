// Copyright (C) 2013 Lile Wang
// For full notice, see "main.cpp" and "COPYING".

#include "driver.h"
#include "kdtree.h"
#include "correlate.h"
#include "read_data.h"
#include "input.h"
#include <iostream>
#include <string>
#include <ctime>

void driver( const std::string & par_file_name )
{
	clock_t start( 0 ), precomp( 0 );
	clock_t autocal( 0 ), crosscal( 0 );
	start = clock(  );
	
	input read_par( par_file_name );
	read_par.read(  );
	read_par.get_init(  );
	
	kdtree data;
	kdtree rand;
	correlate corr;
	read_data read;

	read.set_cosmology( read_par.lambda, read_par.z_max );
	corr.set_dist_bin( read_par.s_max, read_par.s_min,
					   read_par.s_bin );
	
	read.read_from_file( read_par.data_file_name, data );
	data.build_tree(  );
	read.read_from_file( read_par.rand_file_name, rand );
	rand.build_tree(  );
	precomp = clock(  );
	
	corr.gen_bin_counts_auto( data );
	corr.output( read_par.data_file_name + "_ddbins" );
	corr.gen_bin_counts_auto( rand );
	corr.output( read_par.rand_file_name + "_rrbins" );
	autocal = clock(  );
	
	corr.gen_bin_counts_cross( data, rand );
	corr.output( read_par.data_file_name + "_" + 
				 read_par.rand_file_name + "_drbins" );
	crosscal = clock(  );

	std::cout << "Time consumption: "
			  << float( crosscal - start ) / CLOCKS_PER_SEC
			  << " sec."<< std::endl;
	std::cout << "Preprocessing: "
			  << float( precomp - start ) / CLOCKS_PER_SEC
			  << " sec."<< std::endl;
	std::cout << "Auto-correlation: "
			  << float( autocal - precomp ) / CLOCKS_PER_SEC
			  << " sec."<< std::endl;
	std::cout << "Cross-correlation: "
			  << float( crosscal - autocal ) / CLOCKS_PER_SEC
			  << " sec."<< std::endl;
	return;
}

