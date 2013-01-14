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
    void brute_force_1d( const kdtree_node * node0,
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
    static int num_bins;
	static const int dim = 3;
private:                        // Functions
    static int dist_bin_val( double d[  ] );
    int dist_bin( const kdtree_node * node0,
                  const kdtree_node * node1 );
public:
    static void set_dist_bin
    ( double s_max_src, double s_min_src, int num_bins_src );

    ////////// Output //////////
private:
    static std::vector<unsigned> bin_counts_total;
public:
    static bool is_auto_cor;
    void add_to_total(  );
    static void output( std::string file_name );    
	
	////////// 2D distance bin //////////
private:						// Data
	int pi_idx, rp_idx, s_idx;
private:						// Functions
	inline static int idx_2d( int i, int j );
    void brute_force_2d( const kdtree_node * node0,
						 const kdtree_node * node1 );
	void extreme_bin_2d( const double x0[  ],
						 const kdtree_node * node );
public:
	static bool is_2d_cor;
};

#endif

