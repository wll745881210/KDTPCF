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
bool correlate::auto_cor;

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
    bin_counts.resize( num_bins, 0 );
    return;
}

void correlate::static_clear(  )
{
    bin_counts_total.clear(  );
    bin_counts_total.resize( num_bins, 0 );
}

////////////////////////////////////////////////////////////
// Compare trees ( nodes )

void correlate::brute_force_sec
( const kdtree_node * node0, const kdtree_node * node1 )
{
    const std::vector<galaxy_point> & vec0 = *(node0->p_vec);
    const std::vector<galaxy_point> & vec1 = *(node1->p_vec);
    double d_min[ 3 ], d_max[ 3 ], d0( 0. ), d1( 0. );
    int inner_start = node1->idx_start;
    int inner_end   = node1->idx_end;
    const bool not_same_node( node0 != node1 );

    for( int i = node0->idx_start; i <= node0->idx_end; ++ i )
    {
        if( not_same_node )
        {
            for( int k = 0; k < 3; ++ k )
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
            for( int k = 0; k < 3; ++ k )
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
    if( auto_cor )        // If auto-correlation
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

void correlate::cal_corr( const kdtree & tree0,
                           const kdtree & tree1 )
{
    auto_cor = false;
    static_clear(  );
    const kdtree_node * root0 = tree0.get_root_node(  );
    const kdtree_node * root1 = tree1.get_root_node(  );
    if( root0 == root1 )
        auto_cor = true;    
    std::cout << ( auto_cor ? "Auto":"Cross" )
              << "-corr... " << std::flush;
    compare_node( root0, root1 );
    
    std::cout << "Done." << std::endl;
    return;
}

////////////////////////////////////////////////////////////
// Distance bin index calculation

int correlate::dist_bin_val( double d[  ] )
{
    const double s
        = sqrt( d[ 0 ]*d[ 0 ] + d[ 1 ]*d[ 1 ]
                + d[ 2 ]*d[ 2 ] ) / ds;
    return s - s_min_ds;
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
    for( int i = 0; i < num_bins; ++ i )
        bin_counts_total[ i ] += bin_counts[ i ];
    return;
}

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
             << bin_counts_total[ i ] * mult_factor << '\n';
    fout.flush(  );
    return;
}

