// Copyright (C) 2013 Lile Wang
// For full notice, see "main.cpp" and "COPYING".

#include "kdtree.h"

#include <vector>
#include <iostream>
#include <cmath>

////////////////////////////////////////////////////////////
// Class: galaxy point

galaxy_point & galaxy_point::operator =
( const galaxy_point & rhs )
{
    for( int i = 0; i < 3; ++ i )
        this->x[ i ] = rhs.x[ i ];
    return * this;
}

void galaxy_point::swap( galaxy_point & rhs )
{
    for( unsigned i = 0; i < 3; ++ i )
    {
        const double temp = this->x[ i ];
        this->x[ i ] = rhs.x[ i ];
        rhs.x[ i ] = temp;
    }

    return;
}

////////////////////////////////////////////////////////////
// Constructor, destructor and initializer

kdtree::kdtree(  )
{
    root_node = NULL;
}

kdtree::~kdtree(  )
{
    clear( root_node );
}

void kdtree::clear( kdtree_node * & node )
{
    if( node == NULL )
        return;

    clear( node->left );
    clear( node->right );
    delete node;
    node = NULL;
    return;
}

////////////////////////////////////////////////////////////
// Build tree structure

kdtree_node * kdtree::create_node
( kdtree_node * parent_node, int idx_start, int idx_end,
  int depth ) 
{
    if( idx_end < idx_start )
        return NULL;
    
    kdtree_node * current_node = new kdtree_node;
    const int axis = depth % 3;
    const int idx_median
        = select_median( idx_start, idx_end, axis );
    current_node->idx_start = idx_start;
    current_node->idx_end   = idx_end;
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
    return current_node;
}

void kdtree::display_node( kdtree_node * node, int depth )
{
    if( node == NULL )
        return;
    if( depth > 4 )
        return;
    
    for( int i = 0; i < depth; ++ i )
        std::cout << "  ";
    std::cout << "|__"
              << node->idx_end- node->idx_start + 1 << ' '
              << node->max[ 0 ] << ' ' << node->max[ 1 ] << ' '
              << node->max[ 2 ] << '\n';
    
    display_node( node->left, depth + 1 );
    display_node( node->right, depth + 1 );
    return;
}

void kdtree::build_tree(  )
{
    if( source_list.size(  ) < 1 )
        return;
    root_node
        = create_node( NULL, 0,
                       source_list.size(  ) - 1, 0 );
    return;
}

const kdtree_node * kdtree::get_root_node(  ) const
{
    return ( const kdtree_node * ) root_node;
}

void kdtree::display(  )
{
    display_node( root_node, 0 );
    return;
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

