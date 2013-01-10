// Copyright (C) 2013 Lile Wang
// For full notice, see "main.cpp" and "COPYING".

#include "correlate.h"
#include <omp.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <ctime>

////////////////////////////////////////////////////////////
// Constructor, destructor and initializer

correlate::correlate(  )
{
	this->num_threads = 1;
}

correlate::~correlate(  )
{
	
}

void correlate::set_dist_bin( double s_max, double s_min,
							  int num_bins )
{
	this->s_max    = s_max;
	this->s_min    = s_min;
	this->num_bins = num_bins;
	this->ds = ( s_max - s_min ) / double ( num_bins - 1 );
	return;
}

void correlate::set_num_threads( int num_threads )
{
	this->num_threads = num_threads;
	return;
}

void correlate::clear(  )
{
	auto_cor = false;
	bin_counts.clear(  );
	bin_counts.resize( num_bins, 0. );
	return;
}

////////////////////////////////////////////////////////////
// Compare trees ( nodes )

void correlate::compare_node
( const kdtree_node * node0, const kdtree_node * node1 )
{
	if( node0 == NULL || node1 == NULL )
		return;
	if( auto_cor )		// If auto-correlation
	{
		if( node0->idx_start >= node1->idx_end )
			return;
		if( node0 == node1 )
		{
			compare_node( node0->left, node0->left );
			compare_node( node0->right, node0->right );
			compare_node( node0->right, node0->left );
			compare_node( node0->left, node0->right );
			return;
		}
	}
		
	const int idx_bin = dist_bin( node0, node1 );
	switch( idx_bin )
	{
	case -2:
		return;
	case -1:			
		if( node0->idx_end - node0->idx_start >
			node1->idx_end - node1->idx_start )
		{
			compare_node( node0->left, node1 );
			compare_node( node0->right, node1 );		
		}
		else
		{
			compare_node( node0, node1->left );
			compare_node( node0, node1->right );
		}
		break;
	default:
		bin_counts[ idx_bin ]
			+= ( node0->idx_end - node0->idx_start + 1 )
			* ( node1->idx_end - node1->idx_start + 1 );
	}
	return;
}

void correlate::gen_bin_counts_auto( const kdtree & tree0 )
{
	clear(  );
	auto_cor = true;
	std::cout << "Auto-corr pair counting... " << std::flush;

	const kdtree_node * root = tree0.get_root_node(  );
	compare_node( root, root );
	// get_node_vec( root );
	// omp_set_num_threads( this->num_threads );
	// #pragma omp parallel for
	// for( int i = 0; i < int( work_node_vec.size(  ) ); ++ i )
	// 	compare_node( work_node_vec[ i ], root );
	
	std::cout << "Done." << std::endl;
	return;
}

void correlate::gen_bin_counts_cross
( const kdtree & tree0, const kdtree & tree1 )
{
	clear(  );
	auto_cor = false;
	std::cout << "Conducting cross pair counting... ";
	std::cout.flush(  );
	const kdtree_node * root0 = tree0.get_root_node(  );
	const kdtree_node * root1 = tree1.get_root_node(  );
	compare_node( root0, root1 );

	// std::vector<kdtree_node *> parallel_nodes;
	// parallel_nodes.push_back( root0->left );
	// parallel_nodes.push_back( root0->right );

	// omp_set_num_threads( this->num_threads );
	// #pragma omp parallel for
	// for( int i = 0; i < int( work_node_vec.size(  ) ); ++ i )
	// 	compare_node( work_node_vec[ i ], root1 );
	
	std::cout << "Done." << std::endl;
	return;
}

////////////////////////////////////////////////////////////
// Distance bin index calculation

int correlate::dist_bin_val( double d[  ] )
{
	const double s
		= sqrt( pow( d[ 0 ], 2 ) + pow( d[ 1 ], 2 )
				+ pow( d[ 2 ], 2 ) );
	return int( ( s - s_min ) / ds );
}

int correlate::dist_bin( const kdtree_node * node0,
						 const kdtree_node * node1 )
{
	double d0( 0. ), d1( 0. );
	double d_max[ 3 ], d_min[ 3 ];

	for( int i = 0; i < 3; ++ i )
	{
		d0 = node0->max[ i ] - node1->min[ i ];
		d1 = node0->min[ i ] - node1->max[ i ];

		d_max[ i ] = fabs( d0 ) > fabs( d1 ) ? d0 : d1;
		if( d0 * d1 > 0 )
			d_min[ i ] = fabs( d0 ) < fabs( d1 ) ? d0 : d1;
		else
			d_min[ i ] = 0.;
	}

	const int bin_max = dist_bin_val( d_max );
	const int bin_min = dist_bin_val( d_min );

	if( bin_min > num_bins - 1 || bin_max < 0 )
		return -2;
	else if( bin_min != bin_max )
		return -1;
	else
		return bin_min;
}

////////////////////////////////////////////////////////////
// Brute-force method

void correlate::brute_force_ac( const kdtree & tree0 )
{
	clear(  );
	auto_cor = true;
	std::cout << "Brute-force AC..." << std::flush;
	const std::vector<galaxy_point> & src_buf
		= tree0.source_list;

	double d[ 3 ];
	for( unsigned i = 0; i < src_buf.size(  ); ++ i )
		for( unsigned j = i + 1; j < src_buf.size(  ); ++ j )
		{
			for( unsigned k = 0; k < 3; ++ k )
				d[ k ] = src_buf[ i ][ k ] - src_buf[ j ][ k ];
			const int dist_bin = dist_bin_val( d );
			if( dist_bin < 0 || dist_bin > num_bins - 1 )
				continue;
			
			bin_counts[ dist_bin ] += 1;
		}
	std::cout << "Done." << std::endl;
	return;
}

void correlate::brute_force_cc( const kdtree & tree0,
								const kdtree & tree1 )
{
	clear(  );
	auto_cor = false;
	std::cout << "Brute-force CC..." << std::flush;
	const std::vector<galaxy_point> & src_buf0
		= tree0.source_list;
	const std::vector<galaxy_point> & src_buf1
		= tree1.source_list;

	double d[ 3 ];
	for( unsigned i = 0; i < src_buf0.size(  ); ++ i )
		for( unsigned j = 0; j < src_buf1.size(  ); ++ j )
		{
			for( unsigned k = 0; k < 3; ++ k )
				d[ k ] = src_buf0[ i ][k] - src_buf1[ j ][k];
			const int dist_bin = dist_bin_val( d );
			if( dist_bin < 0 || dist_bin > num_bins - 1 )
				continue;
			
			bin_counts[ dist_bin ] += 1;
		}
	std::cout << "Done." << std::endl;
	return;
}

////////////////////////////////////////////////////////////
// Output bin counting results

void correlate::output( std::string file_name )
{
	std::ofstream fout( file_name.c_str(  ) );
	int mult_factor( 0 );
	if( auto_cor )
		mult_factor = 2;
	else
		mult_factor = 1;

	for( int i = 0; i < num_bins; ++ i )
		fout << s_min + ds * ( i + 0.5 ) << '\t'
			 << bin_counts[ i ] * mult_factor << '\n';
	fout.flush(  );
	return;
}

////////////////////////////////////////////////////////////
// Load balancing

void correlate::get_node_vec( const kdtree_node * root )
{
	int max_depth( 0 );
	if( ( ( num_threads ) & ( num_threads - 1 ) ) == 0 )
		max_depth = int( log( num_threads ) / log( 2. ) + 0.5 );
	else
		max_depth = ceil( log( num_threads ) / log( 2. ) ) + 2.;
	
	work_node_vec.clear(  );
	add_work_node( root, max_depth );
	return;
}

void correlate::add_work_node( const kdtree_node * node,
							   int depth_remain )
{
	if( node == NULL )
		return;
	
	if( depth_remain > 0 )
	{
		add_work_node( node->left, depth_remain - 1 );
		add_work_node( node->right, depth_remain - 1 );		
	}
	else
		work_node_vec.push_back( node );
	return;
}

