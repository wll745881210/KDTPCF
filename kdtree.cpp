// Copyright (C) 2013 Lile Wang
// For full notice, see "main.cpp" and "COPYING".

#include "kdtree.h"

#include <vector>
#include <iostream>
#include <cmath>

////////////////////////////////////////////////////////////
// Static variables

int kdtree::jk_depth;

////////////////////////////////////////////////////////////
// Struct: galaxy point

galaxy_point & galaxy_point::operator =
( const galaxy_point & rhs )
{
    for( int i = 0; i < 3; ++ i )
        this->x[ i ] = rhs.x[ i ];
    this->weight = rhs.weight;
    return * this;
}

void galaxy_point::swap( galaxy_point & rhs )
{
    const galaxy_point temp = rhs;
    rhs    = * this;
    * this = temp;
    return;
}

////////////////////////////////////////////////////////////
// Constructor, destructor and initializer

kdtree::kdtree(  )
{
    root_node = NULL;
    jk_sample_count = 0;
    jk_depth = 4;
}

kdtree::~kdtree(  )
{
    clear( root_node );
}

void kdtree::set_jackknife_depth( int jk_d )
{
    jk_depth = jk_d;
    return;
}

void kdtree::clear( kdtree_node * node )
{
    if( node == NULL )
        return;
    clear( node->left );
    clear( node->right );
    delete node;
    return;
}

////////////////////////////////////////////////////////////
// Build tree structure

void kdtree::set_node_weight( kdtree_node * node )
{
    node->weight = 0.;
    if( node->left != NULL && node->right != NULL )
    {
	node->weight += node->left->weight;
	node->weight += node->right->weight;
    }
    else
	for( int i = node->idx_start;
	     i <= node->idx_end; ++ i )
	    node->weight += ( * node->p_vec )[ i ].weight;
    return;
}

kdtree_node * kdtree::create_node
( kdtree_node * parent_node, int idx_start, int idx_end,
  int depth ) 
{    
    kdtree_node * current_node = new kdtree_node;
    const int axis = depth % 3;
    const int idx_median
        = select_median( idx_start, idx_end, axis );
    current_node->idx_start = idx_start;
    current_node->idx_end   = idx_end;
    if( depth == jk_depth )
        current_node->jk_sample = jk_sample_count ++;
    else if( depth > jk_depth )
        current_node->jk_sample = parent_node->jk_sample;
    else
        current_node->jk_sample = -1;
            
    for( int i = 0; i < 3; ++ i )
    {        
        current_node->max[ i ] = coord_max.x[ i ];
        current_node->min[ i ] = coord_min.x[ i ];
    }
    current_node->p_vec = & source_list;
    
    if( current_node->idx_end - current_node->idx_start
        <= leaf_node_num )
    {
        current_node->left  = NULL;
        current_node->right = NULL;
    }
    else
    {        
        current_node->left
	    = create_node( current_node, idx_start,
                           idx_median, depth + 1  );
        current_node->right
            = create_node( current_node, idx_median
                           + 1, idx_end, depth + 1  );
    }
    set_node_weight( current_node );
    return current_node;
}

void kdtree::build_tree(  )
{
    if( source_list.size(  ) < 1 )
        return;
    root_node = create_node
        ( NULL, 0, source_list.size(  ) - 1, 0 );
    return;
}

const kdtree_node * kdtree::get_root_node(  ) const
{
    return ( const kdtree_node * ) root_node;
}

////////////////////////////////////////////////////////////
// Search the median and separate sources

void kdtree::max_min_compare( int idx )
{
    galaxy_point & current_point = source_list[ idx ];
    for( unsigned i = 0; i < 3; ++ i )
    {
        if( current_point.x[ i ] > coord_max.x[ i ] )
            coord_max.x[ i ] = current_point.x[ i ];
        if( current_point.x[ i ] < coord_min.x[ i ] )
            coord_min.x[ i ] = current_point.x[ i ];
    }

    return;
}

void kdtree::max_min_init( int idx_start )
{
    max_min_lock = false;
    coord_max = source_list[ idx_start ];
    coord_min = source_list[ idx_start ];
    return;
}

int kdtree::locate_pivot( int idx_start, int idx_end,
                          int idx_pivot, int axis )
{
    source_list[ idx_pivot ].swap( source_list[ idx_end ] );
    int idx_store = idx_start;

    for( int i = idx_start; i < idx_end; ++ i )
    {
        if( ! max_min_lock )
            max_min_compare( i );
        
        if( source_list[ i ].x[ axis ]
            <= source_list[ idx_end ].x[ axis ] )
        {
            source_list[ i ].swap(source_list[ idx_store ]);
            ++ idx_store;            
        }
    }
    max_min_lock = true;
    source_list[ idx_store ].swap( source_list[ idx_end ] );
    return idx_store;
}

void kdtree::select_kth( int k, int idx_start,
                         int idx_end, int axis )
{
    int idx_pivot
        = locate_pivot( idx_start, idx_end,
                        idx_start, axis );

    if( idx_pivot < k + idx_start )
        select_kth( k + idx_start - idx_pivot - 1,
                    idx_pivot + 1, idx_end, axis );
    else if ( idx_pivot > k + idx_start )
        select_kth( k, idx_start, idx_pivot - 1, axis );

    return;
}

int kdtree::select_median( int idx_start, int idx_end,
                           int axis )
{
    if( idx_end <= idx_start )
    {
        coord_max = source_list[ idx_start ];
        coord_min = source_list[ idx_start ];
        return idx_start;
    }

    max_min_init( idx_start );
    int loc_median = ( idx_end - idx_start ) / 2;
    select_kth( loc_median, idx_start, idx_end, axis );
    return loc_median + idx_start;
}

