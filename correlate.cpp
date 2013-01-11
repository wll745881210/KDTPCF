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
	bin_counts.resize( num_bins + 10, 0. );
	num_total_cal = 0;
	return;
}

////////////////////////////////////////////////////////////
// Compare trees ( nodes )

void correlate::brute_force_sec
( const kdtree_node * node0, const kdtree_node * node1 )
{
	const std::vector<galaxy_point> & vec0 = *(node0->p_vec);
	const std::vector<galaxy_point> & vec1 = *(node1->p_vec);
	double d[ 3 ], e[ 3 ], d0( 0. ), d1( 0. );
	const int inner_start = node1->idx_start;
	const int inner_end = node1->idx_end;

	if( node0 == node1 )
	{
		for( int i = node0->idx_start;
			 i <= node0->idx_end; ++ i )
			for( int j = i + 1;
				 j <= inner_end; ++ j)
			{
				for( int k = 0; k < 3; ++ k )
					d[ k ] = vec0[ i ].x[ k ] - vec0[ j ].x[ k ];
				const int bin_idx = dist_bin_val( d );
				if( bin_idx < 0 || bin_idx > num_bins - 1 )
					continue;
				++ bin_counts[ bin_idx ];
			}
	}
	else
	{
		for( int i = node0->idx_start;
			 i <= node0->idx_end; ++ i )
		{
			for( int k = 0; k < 3; ++ k )
			{
				d0 = vec0[ i ].x[ k ] - node1->max[ k ];
				d1 = vec0[ i ].x[ k ] - node1->min[ k ];
				const bool min_is_d0 = fabs( d0 ) < fabs( d1 );
				e[ k ] = min_is_d0 ? d1 : d0;
				if( d0 * d1 > 0 )
					d[ k ] = min_is_d0 ? d0 : d1;
				else
					d[ k ] = 0.;
			}
			const int min_box_idx = dist_bin_val( d );
			const int max_box_idx = dist_bin_val( e );
			if( min_box_idx > num_bins - 1 || max_box_idx < 0 )
				continue;
			
			if( min_box_idx == max_box_idx )
			{
				bin_counts[ min_box_idx ] += node1->idx_end
					- node1->idx_start + 1;
				continue;
			}
			for( int j = inner_start; j <= inner_end; ++ j)
			{
				for( int k = 0; k < 3; ++ k )
					d[ k ] = vec0[ i ].x[ k ] - vec1[ j ].x[ k ];

				const int bin_idx = dist_bin_val( d );				
				if( bin_idx < 0 || bin_idx > num_bins - 1 )
					continue;
				++ bin_counts[ bin_idx ];
			}
		}
	}
	return;
}

void correlate::compare_node
( const kdtree_node * node0, const kdtree_node * node1 )
{
	if( node0 == NULL || node1 == NULL )
		return;
	if( auto_cor )		// If auto-correlation
	{
		if( node0->idx_start >= node1->idx_end )
			return;
		else if( node0 == node1 )
		{
			if( node0->left == NULL )
			{
				brute_force_sec( node0, node0 );
				return;
			}
			compare_node( node0->left, node0->left );
			compare_node( node0->right, node0->right );
			compare_node( node0->right, node0->left );
			compare_node( node0->left, node0->right );
			return;
		}
	}
		
	const int idx_bin = dist_bin( node0, node1 );
	if( idx_bin == -1 )
	{
		if( node0->left == NULL && node1->left == NULL )
		{
			brute_force_sec( node0, node1 );
			return;
		}
			
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
	}
	else if( idx_bin > -1 )
		bin_counts[ idx_bin ]
			+= ( node0->idx_end - node0->idx_start + 1 )
			* ( node1->idx_end - node1->idx_start + 1 );
	return;
}

void correlate::gen_bin_counts_auto( const kdtree & tree0 )
{
	clear(  );
	auto_cor = true;
	std::cout << "Auto-corr... " << std::flush;

	const kdtree_node * root = tree0.get_root_node(  );
	get_node_vec( root );
	omp_set_num_threads( this->num_threads );
	#pragma omp parallel for
	for( int i = 0; i < int( work_node_vec.size(  ) ); ++ i )
		compare_node( work_node_vec[ i ], root );
	std::cout << "Done." << std::endl;
	return;
}

void correlate::gen_bin_counts_cross
( const kdtree & tree0, const kdtree & tree1 )
{
	clear(  );
	auto_cor = false;
	std::cout << "Cross-corr... " << std::flush;

	const kdtree_node * root0 = tree0.get_root_node(  );
	const kdtree_node * root1 = tree1.get_root_node(  );
	get_node_vec( root0 );
	omp_set_num_threads( this->num_threads );
	#pragma omp parallel for
	for( int i = 0; i < int( work_node_vec.size(  ) ); ++ i )
		compare_node( work_node_vec[ i ], root1 );
	
	std::cout << "Done. " << std::endl;
	return;
}

////////////////////////////////////////////////////////////
// Distance bin index calculation

inline int correlate::dist_bin_val( double d[  ] )
{
	const double s
		= sqrt( d[ 0 ]*d[ 0 ] + d[ 1 ]*d[ 1 ]
				+ d[ 2 ]*d[ 2 ] ) - s_min;
	return int( s / ds );
}

int correlate::dist_bin( const kdtree_node * node0,
						 const kdtree_node * node1 )
{
	double d0, d1;
	double d_max[ 3 ], d_min[ 3 ];
	for( int i = 0; i < 3; ++ i )
	{		
		d0 = node0->max[ i ] - node1->min[ i ];
		d1 = node1->max[ i ] - node0->min[ i ];
		const bool min_is_d0 = fabs( d0 ) < fabs( d1 );	    	
		d_max[ i ] = min_is_d0 ? d1 : d0;
		if( d0 * d1 < 0 )
			d_min[ i ] = min_is_d0 ? d0 : d1;
		else
			d_min[ i ] = 0.;
	}

	const int bin_min = dist_bin_val( d_min );
	if( bin_min > num_bins - 1 )
		return -2;
	const int bin_max = dist_bin_val( d_max );
	return bin_min == bin_max ? bin_min : -1;
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

