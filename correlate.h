// Copyright (C) 2013 Lile Wang
// For full notice, see "main.cpp" and "COPYING".

#ifndef CORRELATE_H_
#define CORRELATE_H_

#include "kdtree.h"
#include <vector>
#include <string>

class correlate
{
    ////////// Con/destructor and initializer //////////
public:
    correlate(  );
    ~correlate(  );
    void clear(  );
    static void static_clear(  );

    ////////// Correlation status //////////
private:
    static bool is_auto_cor, is_2d_cor, is_ang_cor;
public:
    static void set_cor_status( int stat );
    static void set_auto_cor( bool auto_cor );

    ////////// Compare trees //////////
private:                        // Data
    std::vector<unsigned> bin_counts;
private:                        // Functions
    void brute_force( const kdtree_node * node0,
                      const kdtree_node * node1 );
public:
    void compare_node( const kdtree_node * node0,
                       const kdtree_node * node1 );

    ////////// Distance bin index //////////
private:                        // Data
    static double s_max, s_min, ds;
    static int s_bin, phi_bin;
    static const int dim = 3;
private:                        // Functions
    inline static int dist_bin_val( double d[  ] );
    int dist_bin( const kdtree_node * node0,
                  const kdtree_node * node1 );
public:
    static void set_dist_bin
    ( double s_max_src, double s_min_src,
      int s_bin_src, int phi_bin_src );

    ////////// Output //////////
private:
    static std::vector<unsigned> bin_counts_tot;
public:
    void add_to_tot(  );
    static void output( std::string file_name );    
    
    ////////// 2D distance //////////
private:                        // Function
    inline static int idx_2d( int i, int j );
    inline static double dot( const double a[  ],
                              const double b[  ] );

    ////////// Binning index //////////
private:                        // Data
    static double * bin_lim;
private:                        // Function
    static void gen_bin_lim(  );
    static int find_idx( const double & a_2 );

    ////////// Constants //////////
private:
    static const double pi_2 = 1.5707963267948966;
    static const double rad_to_deg = 57.29577951308232;
};

#endif

