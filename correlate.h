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
	void set_dist_bin( double s_max, double s_min,
					   int num_bins );
	void set_num_threads( int num_threads );
	void clear(  );

	////////// Compare trees //////////
private:						// Data
	int num_threads;
	std::vector<unsigned> bin_counts;
private:						// Functions
	void compare_node( const kdtree_node * node0,
					   const kdtree_node * node1 );
public:
	void gen_bin_counts_auto( const kdtree & tree0 );
	void gen_bin_counts_cross( const kdtree & tree0,
							   const kdtree & tree1 );

	////////// Distance bin index //////////
private:						// Data
	double s_max, s_min;
	double ds;
	int num_bins;
private:						// Functions
	int dist_bin_val( double d[  ] );
	int dist_bin( const kdtree_node * node0,
				  const kdtree_node * node1 );

	////////// Brute-force //////////
public:						// Functions
	void brute_force_ac( const kdtree & tree0 );
	void brute_force_cc( const kdtree & tree0,
						 const kdtree & tree1 );

	////////// Output //////////
private:						// Data
	bool auto_cor;
public:
	void output( std::string file_name );

	
	////////// Load balancing //////////
private:						// Data
	std::vector<const kdtree_node *> work_node_vec;
private:						// Function
	void get_node_vec( const kdtree_node * root );
	void add_work_node( const kdtree_node * node,
						int depth_remain );

	////////// Algo, math & constants //////////
private:						// Data
	static const double nearly0 = 1e-6;
	
};

#endif

