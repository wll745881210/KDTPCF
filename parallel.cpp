#include "parallel.h"

#include <omp.h>
#include <cmath>
#include <iostream>

////////////////////////////////////////////////////////////
// Constructor and destructor

parallel::parallel(  ){  }

parallel::~parallel(  ){  }

////////////////////////////////////////////////////////////
// Generate object pool

void parallel::get_node_vec( const kdtree_node * root )
{
    const int max_depth = int( ceil( log( num_threads )
				     / log( 2. ) ) ) + 3;
    work_node_vec.clear(  );
    add_work_node( root, max_depth );

    num_objs = work_node_vec.size(  );
    corr_obj_vec.clear(  );
    corr_obj_vec.resize( num_objs );
    return;
}

void parallel::add_work_node( const kdtree_node * node,
                              int depth_remain )
{
    if( node == NULL )
        return;
    
    if( depth_remain > 0 && node->left != NULL )
    {
        add_work_node( node->left, depth_remain - 1 );
        add_work_node( node->right, depth_remain - 1 );        
    }
    else
        work_node_vec.push_back( node );
    return;
}

void parallel::set_num_threads( int num_threads_src )
{
    num_threads = num_threads_src;
    return;
}

////////////////////////////////////////////////////////////
// Conduct calculation

void parallel::cal_corr( const kdtree & tree0,
                         const kdtree & tree1 )
{    
    correlate::static_clear(  );
    const kdtree_node * root0 = tree0.get_root_node(  );
    const kdtree_node * root1 = tree1.get_root_node(  );

    const bool is_auto_cor( root0 == root1 )
    correlate::set_auto_cor( is_auto_cor );
    std::cout << ( is_auto_cor ? "Auto":"Cross" )
              << "-corr... " << std::flush;
    
    get_node_vec( root0 );
    omp_set_num_threads( num_threads );
    #pragma omp parallel for schedule(dynamic)
    for( int i = 0; i < num_objs; ++ i )
    {
        corr_obj_vec[ i ].clear(  );
        corr_obj_vec[ i ].compare_node
            ( work_node_vec[ i ], root1 );
    }
    #pragma omp parallel for
    for( int i = 0; i < num_objs; ++ i )
        corr_obj_vec[ i ].add_to_tot(  );
	
    const double norm = root0->weight * root1->weight;
    correlate::normalize( norm );

    std::cout << "Done.\n";
    return;
}

