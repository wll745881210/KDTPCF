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

    ////////// Compare trees //////////
private:                        // Data
    std::vector<unsigned> bin_counts;
private:                        // Functions
    void brute_force( const kdtree_node * node0,
					  const kdtree_node * node1 );
public:
    void compare_node( const kdtree_node * node0,
                       const kdtree_node * node1 );
    void cal_corr( const kdtree & tree0,
                   const kdtree & tree1 );

    ////////// Distance bin index //////////
private:                        // Data
    static double s_max, s_min;
    static double ds, s_min_ds;
    static int s_bin;
	static int theta_bin;
	static const int dim = 3;
private:                        // Functions
	inline static int dist_bin_val( double d[  ] );
    int dist_bin( const kdtree_node * node0,
                  const kdtree_node * node1 );
public:
    static void set_dist_bin
    ( double s_max_src, double s_min_src,
	  int s_bin_src, int theta_bin_src );

    ////////// Output //////////
private:
    static std::vector<unsigned> bin_counts_tot;
public:
    static bool is_auto_cor;
    void add_to_tot(  );
    static void output( std::string file_name );    
	
	////////// 2D distance bin //////////
private:						// Data
	static const double pi_2 = 1.5707963267948966;
private:						// Function
	inline static int idx_2d( int i, int j );
	inline static double dot( double a[  ], double b[  ] );
public:
	static bool is_2d_cor;
};

#endif

