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
// Static variables

std::vector<unsigned> correlate::bin_counts_total;
double correlate::s_max;
double correlate::s_min;
double correlate::ds;
double correlate::s_min_ds;
int  correlate::num_bins;
bool correlate::is_auto_cor;
bool correlate::is_2d_cor;

////////////////////////////////////////////////////////////
// Constructor, destructor and initializer

correlate::correlate(  )
{

}

correlate::~correlate(  )
{
    
}

void correlate::clear(  )
{
    bin_counts.clear(  );
    if( is_2d_cor )
		bin_counts.resize( num_bins*num_bins, 0 );
	else
		bin_counts.resize( num_bins, 0 );
    return;
}

void correlate::static_clear(  )
{
    bin_counts_total.clear(  );
    if( is_2d_cor )
		bin_counts_total.resize( num_bins*num_bins, 0 );
	else
		bin_counts_total.resize( num_bins, 0 );
	return;
}

////////////////////////////////////////////////////////////
// Compare trees ( nodes )

void correlate::brute_force_1d
( const kdtree_node * node0, const kdtree_node * node1 )
{
    const std::vector<galaxy_point> & vec0 = *(node0->p_vec);
    const std::vector<galaxy_point> & vec1 = *(node1->p_vec);
    double d_min[ dim ], d_max[ dim ], d0( 0. ), d1( 0. );
    int inner_start = node1->idx_start;
    int inner_end   = node1->idx_end;
    const bool not_same_node( node0 != node1 );

    for( int i = node0->idx_start; i <= node0->idx_end; ++ i )
    {
        if( not_same_node )
        {
            for( int k = 0; k < dim; ++ k )
            {
                d0 = vec0[ i ].x[ k ] - node1->max[ k ];
                d1 = vec0[ i ].x[ k ] - node1->min[ k ];
                const bool min_is_d0 = fabs( d0 ) < fabs( d1 );
                d_max[ k ] = min_is_d0 ? d1 : d0;
                if( d0 * d1 > 0 )
                    d_min[ k ] = min_is_d0 ? d0 : d1;
                else
                    d_min[ k ] = 0.;
            }
            const int min_box_idx = dist_bin_val( d_min );
            if( min_box_idx > num_bins - 1 )
                continue;
            
            if( min_box_idx == dist_bin_val( d_max ) )
            {
                bin_counts[ min_box_idx ] += node1->idx_end
                    - node1->idx_start + 1;
                continue;
            }
        }
        else 
            inner_start = i + 1;
        
        for( int j = inner_start; j <= inner_end; ++ j )
        {
            for( int k = 0; k < dim; ++ k )
                d_min[ k ]
                    = vec0[ i ].x[ k ] - vec1[ j ].x[ k ];
            const int bin_idx = dist_bin_val( d_min );
            if( bin_idx < 0 || bin_idx > num_bins - 1 )
                continue;
            ++ bin_counts[ bin_idx ];
        }
    }
    return;
}

void correlate::compare_node
( const kdtree_node * node0, const kdtree_node * node1 )
{
    if( node0 == NULL || node1 == NULL )
        return;
    if( is_auto_cor )        // If auto-correlation
    {
        if( node0->idx_start >= node1->idx_end )
            return;
        else if( node0 == node1 )
        {
            if( node0->left == NULL )
			{
				if( is_2d_cor )
					brute_force_2d( node0, node1 );
				else
					brute_force_1d( node0, node1 );
			}
			else
            {
				compare_node( node0->left, node0->left );
				compare_node( node0->right, node0->right );
				compare_node( node0->right, node0->left );
				compare_node( node0->left, node0->right );
			}
            return;
        }
    }
        
    s_idx = dist_bin( node0, node1 );
    if( s_idx == -1 || ( is_2d_cor && s_idx > -2 ) )
    {
        if( node0->left == NULL && node1->left == NULL )
        {
			if( is_2d_cor )
				brute_force_2d( node0, node1 );
			else
				brute_force_1d( node0, node1 );
		}
        else if( node0->idx_end - node0->idx_start >
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
    else if( s_idx > -1 && ( !is_2d_cor ) )
        bin_counts[ s_idx ]
            += ( node0->idx_end - node0->idx_start + 1 )
            * ( node1->idx_end - node1->idx_start + 1 );
    return;
}

void correlate::cal_corr( const kdtree & tree0,
						  const kdtree & tree1 )
{
    static_clear(  );
    const kdtree_node * root0 = tree0.get_root_node(  );
    const kdtree_node * root1 = tree1.get_root_node(  );
	is_auto_cor = ( root0 == root1 );
    std::cout << ( is_auto_cor ? "Auto":"Cross" )
              << "-corr... " << std::flush;
    compare_node( root0, root1 );
    add_to_total(  );
    std::cout << "Done." << std::endl;
    return;
}

////////////////////////////////////////////////////////////
// Distance bin index calculation

int correlate::dist_bin_val( double d[  ] )
{
	return sqrt( pow( d[ 0 ], 2 ) + pow( d[ 1 ], 2 )
				 + pow( d[ 2 ], 2 ) ) / ds - s_min_ds;
}

int correlate::dist_bin( const kdtree_node * node0,
                         const kdtree_node * node1 )
{
    double d0, d1;
    double d_max[ dim ], d_min[ dim ];
    for( int i = 0; i < dim; ++ i )
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

void correlate::set_dist_bin
( double s_max_src, double s_min_src, int num_bins_src )
{
    s_max    = s_max_src;
    s_min    = s_min_src;
    num_bins = num_bins_src;
    ds = ( s_max - s_min ) / double ( num_bins - 1 );
    s_min_ds = s_min / ds;
    bin_counts_total.clear(  );
    bin_counts_total.resize( num_bins, 0 );
    return;
}

////////////////////////////////////////////////////////////
// Output bin counting results

void correlate::add_to_total(  )
{
    for( unsigned i = 0; i < bin_counts.size(  ); ++ i )
        bin_counts_total[ i ] += bin_counts[ i ];
    return;
}

void correlate::output( std::string file_name )
{
	double s( 0. ), rp( 0. ), pi( 0. );
    std::ofstream fout( file_name.c_str(  ) );
    const int mult_f = ( is_auto_cor ? 2 : 1 );

	if( is_2d_cor )
	{
		for( int i = 0; i < num_bins; ++ i )
		{
			rp = s_min + ds * ( i + 0.5 );
			for( int j = 0; j < num_bins; ++ j )
			{
				pi = s_min + ds * ( j + 0.5 );
				s = sqrt( rp*rp + pi*pi );
				if( s > s_max )
					continue;
				fout << rp << '\t' << pi << '\t'
					 << mult_f * bin_counts_total[ idx_2d( i, j ) ]
					 << '\n';
			}
		}
	}
	else
		for( int i = 0; i < num_bins; ++ i )
			fout << s_min + ds * ( i + 0.5 ) << '\t'
				 << bin_counts_total[ i ] * mult_f << '\n';
	
    return;
}

////////////////////////////////////////////////////////////
// Functions specified for 2D calculations

int correlate::idx_2d( int i, int j )
{
	return i * num_bins + j;
}

void correlate::extreme_bin_2d
( const double x0[  ], const kdtree_node * node )
{
	double s_2_min( 0. ), s_2_max( 0. );
	double los_2_min( 0. ), los_2_max( 0. );
	double pi_numer_min( 0. ), pi_numer_max( 0. );
	double d0( 0. ), d1( 0. ), e0( 0. ), e1( 0. );
	bool min_is_d0( false ), min_is_e0( false );

	for( int i = 0; i < dim; ++ i )
	{
		d0 = x0[ i ] - node->max[ i ];
		d1 = x0[ i ] - node->min[ i ];
		min_is_d0 = fabs( d0 ) < fabs( d1 );
		s_2_max += min_is_d0 ? d1*d1 : d0*d0;
		if( d0 * d1 > 0 )
			s_2_min += min_is_d0 ? d0*d0 : d1*d1;

		e0 = x0[ i ] + node->max[ i ];
		e1 = x0[ i ] + node->min[ i ];
		min_is_e0 = fabs( e0 ) < fabs( e1 );
		los_2_max += min_is_e0 ? e1*e1 : e0*e0;
		if( e0 * e1 > 0 )
			los_2_min += min_is_e0 ? e0*e0 : e1*e1;		

		e0 = e0 * d0;
		e1 = e1 * d1;
		min_is_e0 = e0 < e1;	// Caution: No fabs!
		pi_numer_max += min_is_e0 ? e1 : e0;
		pi_numer_min += min_is_e0 ? e0 : e1;
	}
	
	const int s_idx_min = sqrt( s_2_min ) / ds - s_min_ds;
	if( s_idx_min > num_bins - 1 )
	{
		s_idx = -2;
		return;
	}
	const int s_idx_max = sqrt( s_2_max ) / ds - s_min_ds;
	if( s_idx_max != s_idx_min )
	{
		s_idx = -1;
		return;
	}

	if( pi_numer_max * pi_numer_min < 0 )
		pi_numer_min = 0.;
	const double pi_max
		= fabs( pi_numer_max ) / sqrt( los_2_min );
	const double pi_min
		= fabs( pi_numer_min ) / sqrt( los_2_max );
	const int pi_idx_max = pi_max / ds - s_min_ds;
	const int pi_idx_min = pi_min / ds - s_min_ds;
	if( pi_idx_min != pi_idx_max )
	{
		s_idx = -1;
		return;
	}

	const int rp_idx_max = sqrt( s_2_max - pi_min*pi_min )
		/ ds - s_min_ds;
	const int rp_idx_min = sqrt( s_2_min - pi_max*pi_max )
		/ ds - s_min_ds;
    if( rp_idx_max != rp_idx_min )
	{
		s_idx = -1;
		return;
	}
    pi_idx = pi_idx_min;
    rp_idx = rp_idx_min;
    s_idx  = s_idx_min;

    return;
}

void correlate::brute_force_2d
( const kdtree_node * node0, const kdtree_node * node1 )
{
    const std::vector<galaxy_point> & vec0 = *(node0->p_vec);
    const std::vector<galaxy_point> & vec1 = *(node1->p_vec);
    int inner_start = node1->idx_start;
    int inner_end   = node1->idx_end;
    const bool not_same_node( node0 != node1 );
	double los_2( 0. ), s_2( 0. ), pi( 0. ), rp( 0. );
	double temp_s( 0. ), temp_los( 0. ), temp_pi( 0. );

    for( int i = node0->idx_start; i <= node0->idx_end; ++ i )
    {
        if( not_same_node )
        {
			extreme_bin_2d( vec0[ i ].x, node1 );
			if( s_idx < -1 )
				continue;
            if( s_idx > -1 )
            {
                bin_counts[ idx_2d( pi_idx, rp_idx ) ]
					+= node1->idx_end - node1->idx_start + 1;
                continue;
            }
        }
        else 
            inner_start = i + 1;
        
        for( int j = inner_start; j <= inner_end; ++ j )
        {
            for( int k = 0; k < dim; ++ k )
            {
				temp_s = vec0[ i ].x[ k ] - vec1[ j ].x[ k ];
				s_2	+= pow( temp_s, 2 );
				temp_los = vec0[ i ].x[ k ] + vec1[ j ].x[ k ];
				los_2 += pow( temp_los, 2 );
				temp_pi += temp_s * temp_los;
			}
			if( s_2 > s_max * s_max )
				continue;
			pi = fabs( temp_pi / sqrt( los_2 ) );
			rp = sqrt( s_2 - pi*pi );
			pi_idx = pi / ds - s_min_ds;
			rp_idx = rp / ds - s_min_ds;
            ++ bin_counts[ idx_2d( pi_idx, rp_idx ) ];
        }
    }
    return;
}

