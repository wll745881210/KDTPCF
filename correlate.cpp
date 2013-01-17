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

std::vector<unsigned> correlate::bin_counts_tot;
double correlate::s_max;
double correlate::s_min;
double correlate::ds;
double * correlate::bin_lim;
int  correlate::s_bin;
int  correlate::phi_bin;
bool correlate::is_auto_cor;
bool correlate::is_2d_cor;
bool correlate::is_ang_cor;

// double correlate::time_waster_d;


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
        bin_counts.resize( s_bin*phi_bin, 0 );
    else
        bin_counts.resize( s_bin, 0 );
    return;
}

void correlate::static_clear(  )
{
    bin_counts_tot.clear(  );
    if( is_2d_cor )
        bin_counts_tot.resize( s_bin*phi_bin, 0 );
    else
        bin_counts_tot.resize( s_bin, 0 );        
    return;    
}

////////////////////////////////////////////////////////////
// Set correlation status

void correlate::set_cor_status( int stat )
{
    if( stat < 0 || stat > 2 )
        throw "Incorrect status";
    is_2d_cor  = ( stat == 2 );
    is_ang_cor = ( stat == 0 );
    return;
}

void correlate::set_auto_cor( bool auto_cor )
{
    is_auto_cor = auto_cor;
    return;
}

////////////////////////////////////////////////////////////
// Compare trees ( nodes )

void correlate::brute_force
( const kdtree_node * node0, const kdtree_node * node1 )
{
    const std::vector<galaxy_point> & vec0 = *(node0->p_vec);
    const std::vector<galaxy_point> & vec1 = *(node1->p_vec);
    double d0, d1, d_max[ dim ], d_min[ dim ];
    int inner_start = node1->idx_start;
    int inner_end   = node1->idx_end;
    const bool is_same_node( node0 == node1 );

    for( int i = node0->idx_start; i <= node0->idx_end; ++ i )
    {
        if( is_same_node )
            inner_start = i + 1;
        else
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
            if( min_box_idx > s_bin - 1 )
                continue;
            const int max_box_idx = dist_bin_val( d_max );
            if( min_box_idx == max_box_idx )
            {
                bin_counts[ min_box_idx ] += node1->idx_end
                    - node1->idx_start + 1;
                continue;
            }
        }
        
        for( int j = inner_start; j <= inner_end; ++ j )
        {
            if( is_2d_cor )
            {
                for( int k = 0; k < dim; ++ k )
                {
                    const double a = vec0[ i ].x[ k ];
                    const double b = vec1[ j ].x[ k ];
                    d_min[ k ] = a - b;
                    d_max[ k ] = a + b;
                }
                const double s_2 = dot( d_min, d_min );
                const int s_idx  = ( sqrt( s_2 ) - s_min ) / ds;
                if( s_idx > s_bin - 1 )
                    continue;
                const double los_2 = dot( d_max, d_max );
                const double phi = fabs( dot( d_max, d_min ) )
                    / sqrt( los_2 * s_2 );
                if( phi > 1. )
                    continue;
                
                const int phi_idx = acos( phi ) * phi_bin / pi_2;
                ++ bin_counts[ idx_2d( s_idx, phi_idx ) ];
            }
            else
            {
                for( int k = 0; k < dim; ++ k )
                    d_min[ k ] = vec0[ i ].x[ k ]
                        - vec1[ j ].x[ k ];
                const int bin_idx = dist_bin_val( d_min );
                if( bin_idx < 0 || bin_idx > s_bin - 1 )
                    continue;
                ++ bin_counts[ bin_idx ];
            }
        }
    }
    return;
}

void correlate::compare_node     //  Recursion is NOT slow!
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
                brute_force( node0, node1 );
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
        
    const int s_idx = dist_bin( node0, node1 );
    if( s_idx == -1 || ( is_2d_cor && ( s_idx > -2 ) ) )
    {
        if( node0->left == NULL && node1->left == NULL )
            brute_force( node0, node1 );
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
    else if( s_idx > -1 )
        bin_counts[ s_idx ]
            += ( node0->idx_end - node0->idx_start + 1 )
            * ( node1->idx_end - node1->idx_start + 1 );
    return;
}

////////////////////////////////////////////////////////////
// Distance bin index calculation

int correlate::dist_bin_val( double d[  ] )
{
    return find_idx( pow( d[ 0 ], 2 ) + pow( d[ 1 ], 2 )
                    + pow( d[ 2 ], 2 ) );
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
    if( bin_min > s_bin - 1 )
        return -2;
    const int bin_max = dist_bin_val( d_max );
    return bin_min == bin_max ? bin_min : -1;
}

void correlate::set_dist_bin
( double s_max_src, double s_min_src,
  int s_bin_src, int phi_bin_src )
{
    s_max   = s_max_src;
    s_min   = s_min_src;
    s_bin   = s_bin_src;
    phi_bin = phi_bin_src;
    ds = ( s_max - s_min ) / double ( s_bin - 1 );
    gen_bin_lim(  );
    return;
}

////////////////////////////////////////////////////////////
// Output bin counting results

void correlate::add_to_tot(  )
{
    for( unsigned i = 0; i < bin_counts.size(  ); ++ i )
        bin_counts_tot[ i ] += bin_counts[ i ];
    return;
}

void correlate::output( std::string file_name )
{
    double s( 0. ), phi( 0. );
    std::ofstream fout( file_name.c_str(  ) );
    const int f = ( is_auto_cor ? 2 : 1 );

    if( is_2d_cor )
        for( int i = 1 - s_bin; i < s_bin; ++ i )
        {
            s = s_min + ds * ( i + 0.5 );
            for( int j = 1 - phi_bin; j < phi_bin; ++ j )
            {
                phi = pi_2 * ( j + 0.5 ) / double( phi_bin );
                const int i_abs = i > 0 ? i : -i;
                const int j_abs = j > 0 ? j : -j;
                int k = idx_2d( i_abs, j_abs );
                fout << s << '\t' << phi * rad_to_deg << '\t'
                     << f * bin_counts_tot[ k ] << '\n';
            }
        }
    else
        for( int i = 0; i < s_bin; ++ i )
            fout << s_min + ds * ( i + 0.5 ) << '\t'
                 << bin_counts_tot[ i ] * f << '\n';
    return;
}

////////////////////////////////////////////////////////////
// Functions specified for 2D calculations

int correlate::idx_2d( int i, int j )
{
    return i * phi_bin + j;
}

double correlate::dot( const double a[  ],
                       const double b[  ] )
{
    return a[ 0 ]*b[ 0 ] + a[ 1 ]*b[ 1 ] + a[ 2 ]*b[ 2 ];
}

////////////////////////////////////////////////////////////
// Finding binning index more quickly

void correlate::gen_bin_lim(  )
{
    bin_lim = new double [ s_bin ];
    if( is_ang_cor )
    {
        for( int i = 0; i < s_bin; ++ i )
        {
            const double theta
                = ( s_min + i * ds ) / rad_to_deg;
            bin_lim[ i ] = pow( 2. * sin( theta / 2. ), 2 );
        }
        s_max = pow( 2. * sin( ( s_max + ds )
                               / rad_to_deg / 2. ), 2 );
    }
    else
    {
        for( int i = 0; i < s_bin; ++ i )
        {
            const double s = s_min + i * ds;
            bin_lim[ i ] = pow( s, 2 );
        }
        s_max = pow( s_max + ds, 2 );
    }
    return;
}

int correlate::find_idx( const double & a_2 )
{
    if( a_2 > s_max )
        return s_bin;
    for( int i = s_bin - 1; i > -1; -- i )
        if( bin_lim[ i ] <= a_2 )
            return i;
    return 0;
}

