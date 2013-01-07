// Copyright (C) 2013 Lile Wang
// For full notice, see "main.cpp" and "COPYING".

#ifndef KDTREE_H_
#define KDTREE_H_

#include <vector>
#include <string>

////////////////////////////////////////////////////////////
// Auxillary data types

struct galaxy_point
{
	double x, y, z;
	galaxy_point & operator = ( const galaxy_point & rhs );
	double & operator[] ( const unsigned & idx );
	const double & operator[] ( const unsigned & idx ) const;
	void swap( galaxy_point & rhs );
};

struct kdtree_node
{	
	kdtree_node * left;
	kdtree_node * right;

	int num_nodes;
	galaxy_point max, min;
};

class read_data;


////////////////////////////////////////////////////////////
// K-d Tree building

class kdtree
{
	////////// Con/Destructor & Initializer //////////
public:
	kdtree(  );
	~kdtree(  );
private:
	void clear( kdtree_node * & node );

	////////// Point data buffer //////////
private:						// Data
	std::vector<galaxy_point> source_list;
	friend class read_data;

	////////// Tree structure //////////
private:						// Data	
	kdtree_node * root_node;
	int max_depth;
private:						// Function
	kdtree_node * create_node( kdtree_node * parent_node,
							   int idx_start,
							   int idx_end, int depth );
	void display_node( kdtree_node * node, int depth );
public:
	void build_tree(  );
	const kdtree_node * get_root_node(  ) const;
	void display(  );
	
	////////// Median search //////////
private:						// Data
	bool max_min_lock;
	galaxy_point coord_max, coord_min;
private:						// Function
	void max_min_compare( int idx );
	void max_min_init( int idx_start );
	int locate_pivot( int idx_start, int idx_end,
					  int idx_pivot, int axis );
	void select_kth( int k, int idx_start,
					 int idx_end, int axis );
	int select_median( int idx_start,
					   int idx_end, int axis );
};

#endif

