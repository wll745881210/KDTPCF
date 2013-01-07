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
	void clear(  );

	////////// Compare trees //////////
private:						// Data
	std::vector<double> bin_counts;
private:						// Functions
	void compare_node_auto( const kdtree_node * node0,
							const kdtree_node * node1 );
	void compare_node_cross( const kdtree_node * node0,
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

	////////// Output //////////
private:						// Data
	char cor_type;
public:
	void output( std::string file_name );

	////////// Algo, math & constants //////////
private:						// Data
	static const double nearly0 = 1e-6;
	
};

#endif

