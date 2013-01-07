// Copyright (C) 2013 Lile Wang
// For full notice, see "main.cpp" and "COPYING".

#include "correlate.h"
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

void correlate::clear(  )
{
	cor_type = '\0';
	bin_counts.clear(  );
	bin_counts.resize( num_bins, 0. );
	return;
}

////////////////////////////////////////////////////////////
// Compare trees ( nodes )

void correlate::compare_node_auto
( const kdtree_node * node0, const kdtree_node * node1 )
{
	if( node0 == NULL || node1 == NULL || node0 < node1 )
		return;
	if( node0 == node1 )
	{
		compare_node_auto( node0->left, node0->right );
		compare_node_auto( node0->right, node0->left );
		compare_node_auto( node0->left, node0->left );
		compare_node_auto( node0->right, node0->right );
		return;
	}
		
	const int idx_bin = dist_bin( node0, node1 );
	switch( idx_bin )
	{
	case -2:
		return;
	case -1:
		if( node0->num_nodes > node1->num_nodes )
		{
			compare_node_auto( node0->left, node1 );
			compare_node_auto( node1, node0->left );
			compare_node_auto( node0->right, node1 );		
			compare_node_auto( node1, node0->right );
		}
		else
		{
			compare_node_auto( node1->left, node0 );
			compare_node_auto( node0, node1->left );
			compare_node_auto( node1->right, node0 );
			compare_node_auto( node0, node1->right );
		}
		break;
	default:
		bin_counts[ idx_bin ] += node0->num_nodes
			* node1->num_nodes;
	}
	return;
}

void correlate::compare_node_cross
( const kdtree_node * node0, const kdtree_node * node1 )
{
	if( node0 == NULL || node1 == NULL )
		return;
		
	const int idx_bin = dist_bin( node0, node1 );
	switch( idx_bin )
	{
	case -2:
		return;
	case -1:
		if( node0->num_nodes > node1->num_nodes )
		{
			compare_node_cross( node0->left, node1 );
			compare_node_cross( node0->right, node1 );		
		}
		else
		{
			compare_node_cross( node1->left, node0 );
			compare_node_cross( node1->right, node0 );
		}
		break;
	default:
		bin_counts[ idx_bin ] += node0->num_nodes
			* node1->num_nodes;
	}
	return;
}

void correlate::gen_bin_counts_auto( const kdtree & tree0 )
{
	clear(  );
	cor_type = 'a';
	std::cout << "Conducting auto pair counting... ";
	std::cout.flush(  );
	const kdtree_node * root = tree0.get_root_node(  );
	compare_node_auto( root, root );
	std::cout << "Done." << std::endl;
	return;
}

void correlate::gen_bin_counts_cross
( const kdtree & tree0, const kdtree & tree1 )
{
	clear(  );
	cor_type = 'c';
	std::cout << "Conducting cross pair counting... ";
	std::cout.flush(  );
	const kdtree_node * root0 = tree0.get_root_node(  );
	const kdtree_node * root1 = tree1.get_root_node(  );
	compare_node_cross( root0, root1 );
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
// Output bin counting results

void correlate::output( std::string file_name )
{
	std::ofstream fout( file_name.c_str(  ) );
	int mult_factor( 0 );
	if( cor_type == 'a' )
		mult_factor = 2;
	else
		mult_factor = 1;

	for( int i = 0; i < num_bins; ++ i )
		fout << s_min + ds * ( i + 0.5 ) << '\t'
			 << bin_counts[ i ] * mult_factor << '\n';
	fout.flush(  );
	return;
}

