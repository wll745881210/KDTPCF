// Copyright (C) 2013 Lile Wang
// For full notice, see "main.cpp" and "COPYING".

#include "correlate.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <ctime>
#include <omp.h>

////////////////////////////////////////////////////////////
// Static variables

correlate::dvec correlate::bin_counts_tot;
correlate::dvec correlate::bin_counts_tot_jk;
double correlate::s_max, correlate::s_min, correlate::ds;
double * correlate::s_bin_lim( NULL );
double * correlate::phi_bin_lim( NULL );
int correlate::s_num, correlate::phi_num, correlate::jk_num;
bool correlate::is_log_bin( false );
bool correlate::is_regular_phi_bin( false );
bool correlate::is_auto_cor( false );
bool correlate::is_2d_cor( false );
bool correlate::is_ang_cor( false );

const double correlate::pi_2      ( 1.5707963267948966);
const double correlate::rad_to_deg( 57.29577951308232 );
const int    correlate::dim       ( 3                 ); 
const double correlate::tiny      ( 1.e-2             );

////////////////////////////////////////////////////////////
// Constructor, destructor and initializer

correlate::correlate(  ){  }

correlate::~correlate(  ){  }

void correlate::clear(  )
{
    bin_counts.clear(  );
    bin_counts_jk.clear(  );
    const int size = is_2d_cor ? s_num*phi_num : s_num;
    bin_counts.resize( size, 0. );
    bin_counts_jk.resize( size * jk_num, 0. );
    return;
}

void correlate::static_clear(  )
{
    bin_counts_tot.clear(  );
    bin_counts_tot_jk.clear(  );
    const int size = is_2d_cor ? s_num*phi_num : s_num;
    bin_counts_tot.resize( size, 0. );
    bin_counts_tot_jk.resize( size * jk_num, 0. );
    return;    
}

void correlate::set_auto_cor( bool auto_cor )
{
    is_auto_cor = auto_cor;
    return;
}

void correlate::set_par
( double s_max_s, double s_min_s, int s_num_s,
  bool log_bin, int phi_num_s, bool regular_phi_bin,
  int corr_stat, int jk_n )
{
    if( corr_stat < 0 || corr_stat > 2 )
        throw "Incorrect status";
    is_2d_cor	       = ( corr_stat == 2 );
    is_ang_cor	       = ( corr_stat == 0 );
    is_log_bin	       = log_bin;
    is_regular_phi_bin = regular_phi_bin;
    
    s_max = is_log_bin ? log(s_max_s) : s_max_s;
    if( is_log_bin && s_min_s < tiny )
        s_min_s = tiny;
    s_min = is_log_bin ? log(s_min_s) : s_min_s;
    s_num = s_num_s;
    ds = ( s_max - s_min ) / double( s_num );
    s_bin_lim = new double [ s_num + 1 ];
    for( int i = 0; i < s_num + 1; ++ i )
        s_bin_lim[ i ] = s_lim_val( i );

    phi_num = phi_num_s;
    phi_bin_lim = new double [ phi_num + 1 ];
    for( int i = 0; i < phi_num + 1; ++ i )
        phi_bin_lim[ i ] = phi_lim_val( i );

    jk_num = jk_n;
    return;
}

////////////////////////////////////////////////////////////
// Compare trees ( nodes )


void correlate::brute_force        // Only when necessary
( const kdtree_node * node0, const kdtree_node * node1 )
{
    const std::vector<galaxy_point> & vec0 = *(node0->p_vec);
    const std::vector<galaxy_point> & vec1 = *(node1->p_vec);
    double d0, d1, d_max[ dim ], d_min[ dim ];
    int inner_start = node1->idx_start;
    int inner_end   = node1->idx_end;
    int idx_tot( 0 );
    const bool is_same_node( node0 == node1 );
    const int sample0( node0->jk_sample );
    const int sample1( node1->jk_sample );
    const bool is_valid_jk_node0( sample0 >= 0 );
    const bool is_valid_diff_jk
        = ( sample1 >= 0 ) && ( sample0 != sample1 );

    for( int i = node0->idx_start;
	 i <= node0->idx_end; ++ i )
    {
        if( is_same_node )
            inner_start = i + 1;
        else
        {
            for( int k = 0; k < dim; ++ k )
            {
                d0 = vec0[ i ].x[ k ] - node1->max[ k ];
                d1 = vec0[ i ].x[ k ] - node1->min[ k ];
                const bool min_is_d0
		    = fabs( d0 ) < fabs( d1 );
                d_max[ k ] = min_is_d0 ? d1 : d0;
                if( d0 * d1 > 0 )
                    d_min[ k ] = min_is_d0 ? d0 : d1;
                else	      // Overlap in this dimension
                    d_min[ k ] = 0.;
            }
            const int min_box_idx = s_idx_arr( d_min );
            if( min_box_idx > s_num - 1 )
                continue;
            
            if( min_box_idx == s_idx_arr( d_max )
                && min_box_idx > 0 && ( ! is_2d_cor ) )
            {
		const double add = node1->weight
		    * vec0[ i ].weight;
                bin_counts[ min_box_idx ] += add;
                if( is_valid_jk_node0 )
                    jk_add( min_box_idx, sample0, add );
                if( is_valid_diff_jk )
                    jk_add( min_box_idx, sample1, add );
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
                const int s_idx = s_idx_val( s_2 );
                if( s_idx > s_num - 1 || s_idx < 0 )
                    continue;
                const double los_2 = dot( d_max, d_max );
                const double mu = fabs( dot( d_max, d_min ) )
                    / sqrt( los_2 * s_2 );
                const int phi_idx = phi_idx_val( mu );
                idx_tot = idx_2d( s_idx, phi_idx );
            }
            else
            {
                for( int k = 0; k < dim; ++ k )
                    d_min[ k ] = vec0[ i ].x[ k ]
                        - vec1[ j ].x[ k ];
                idx_tot = s_idx_arr( d_min );
                if( idx_tot > s_num - 1 || idx_tot < 0 )
                    continue;
            }
	    const double add = vec0[ i ].weight
		* vec1[ j ].weight;
	    bin_counts[ idx_tot ] += add;
	    if( is_valid_jk_node0 )
		jk_add( idx_tot, sample0, add );
	    if( is_valid_diff_jk )
		jk_add( idx_tot, sample1, add );
        }
    }
    return;
}

int correlate::node_bin( const kdtree_node * node0,
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
        else		         // Overlap in this dimension
            d_min[ i ] = 0.;
    }
    
    const int bin_min = s_idx_arr( d_min );
    if( bin_min > s_num - 1 )
        return -2;
    return bin_min == s_idx_arr( d_max ) ? bin_min : -1;
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

    const int s_idx = node_bin( node0, node1 );    
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
    {
        const double add = node0->weight * node1->weight;
        bin_counts[ s_idx ] += add;
        const int sample0( node0->jk_sample );
        const int sample1( node1->jk_sample );        
        if( sample0 >= 0 )
            jk_add( s_idx, sample0, add );
        if( ( sample0 != sample1 ) && ( sample1 >= 0 ) )
            jk_add( s_idx, sample1, add );
    }
    return;
}

////////////////////////////////////////////////////////////
// Jackknife

void correlate::jk_add( int idx, int sample, double add )
{
    bin_counts_jk[ idx * jk_num + sample ] += add;
    return;
}

////////////////////////////////////////////////////////////
// Output bin counting results

void correlate::add_to_tot(  )
{
    for( unsigned i = 0; i < bin_counts.size(  ); ++ i )
        bin_counts_tot[ i ] += bin_counts[ i ];
    for( unsigned i = 0; i < bin_counts_jk.size(  ); ++ i )
        bin_counts_tot_jk[ i ] += bin_counts_jk[ i ];
    return;
}

void correlate::normalize( double factor )
{
    for( unsigned i = 0; i < bin_counts_tot.size(  ); ++ i )
	bin_counts_tot[ i ] /= factor;
    for( unsigned i = 0; i < bin_counts_tot_jk.size(  ); ++ i )
	bin_counts_tot_jk[ i ] /= factor;
    return;
}

void correlate::out_one_line( std::ofstream & fout, int idx )
{
    const int f = ( is_auto_cor ? 2 : 1 );
    fout << '\t' << f * bin_counts_tot[ idx ];
    for( int j = 0; j < jk_num; ++ j )
    {
        const double jk_res = bin_counts_tot[ idx ]
            - bin_counts_tot_jk[ j + jk_num * idx ];
        fout << '\t' << f * jk_res;
    }
    fout << '\n';
    return;
}

void correlate::output( std::string file_name )
{
    std::ofstream fout( file_name.c_str(  ) );
    for( int i = 0; i < s_num; ++ i )
    {
        const double s_bin_min = s_val( i, 0. );
        const double s_bin_cen = s_val( i );
        const double s_bin_max = s_val( i, 1. );
        if( is_2d_cor )
            for( int j = 0; j < phi_num; ++ j )
            {
                const double phi_bin_min = phi_val( j, 0. );
                const double phi_bin_cen = phi_val( j );
                const double phi_bin_max = phi_val( j, 1. );
                const int k = idx_2d( i, j );
                fout << s_bin_min << '\t' << s_bin_cen << '\t'
                     << s_bin_max << '\t' << phi_bin_min << '\t'
                     << phi_bin_cen << '\t' << phi_bin_max;
                out_one_line( fout, k );
            }
        else
        {
            fout << s_bin_min << '\t' << s_bin_cen << '\t'
                 << s_bin_max;
            out_one_line( fout, i );
        }            
    }
    return;
}

const correlate::dvec & correlate::bin_count_ref(  )
{
    return bin_counts_tot;
}

const correlate::dvec & correlate::bin_count_jk_ref(  )
{
    return bin_counts_tot_jk;
}

////////////////////////////////////////////////////////////
// Functions specified for 2D calculations

int correlate::idx_2d( int i, int j )
{
    return i * phi_num + j;
}

double correlate::dot( const double a[  ],
                       const double b[  ] )
{
    return a[ 0 ]*b[ 0 ] + a[ 1 ]*b[ 1 ] + a[ 2 ]*b[ 2 ];
}

////////////////////////////////////////////////////////////
// Find binning index more quickly: Look up the table!

double correlate::s_lim_val( int i )
{
    double s = s_min + i * ds;
    s = is_log_bin ? exp( s ) : s;
    if( is_ang_cor )
        return pow( 2. * sin( s / rad_to_deg / 2. ), 2 );
    return pow( s, 2 );
}

double correlate::phi_lim_val( int i )
{
    if( is_regular_phi_bin )
        return cos( pi_2 * i / double( phi_num ) );
    else
        return 1. - i / double( phi_num );
}

int correlate::s_idx_val( const double & a_2 )
{
    if( s_bin_lim[ s_num ] <= a_2 )
        return s_num;
    else if( s_bin_lim[ 0 ] > a_2 )
        return -1;    
    int i( 0 ), j( s_num );
    while( j - i > 1 )
    {
        const int k = ( i + j ) >> 1;
        s_bin_lim[ k ] < a_2 ? ( i = k ) : ( j = k );
    }
    return i;
}

int correlate::s_idx_arr( const double a[  ] )
{
    return s_idx_val( a[ 0 ]*a[ 0 ] + a[ 1 ]*a[ 1 ]
                      + a[ 2 ]*a[ 2 ] );
}

int correlate::phi_idx_val( const double & mu )
{    
    if( mu >= 1. )
        return phi_num - 1;
    int i( 0 ), j( phi_num );
    while( j - i > 1 )
    {
        const int k = ( i + j ) >> 1;
        phi_bin_lim[ k ] > mu ? ( i = k ) : ( j = k );
    }
    return i;
}

double correlate::s_val( int i, double offset )
{
    const int i_sgn = i >= 0 ? 1 : -1;
    const double s = s_min + ( i * i_sgn + offset ) * ds;
    return ( is_log_bin ? exp( s ) : s ) * i_sgn;
}

double correlate::phi_val( int i, double offset )
{
    const double idx_ratio
	= ( i + offset ) / double(phi_num);
    if( is_regular_phi_bin )
        return pi_2 * idx_ratio * rad_to_deg;
    else
        return rad_to_deg * acos( 1. - idx_ratio );
}
