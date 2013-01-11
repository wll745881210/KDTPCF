#ifndef PARALLEL_H_
#define PARALLEL_H_

#include "correlate.h"
#include "kdtree.h"
#include <string>

class parallel
{
    ////////// Con/destructor //////////
public:
    parallel(  );
    ~parallel(  );


    ////////// Object pool //////////
private:                        // Data
    int num_threads;
    int num_objs;
    std::vector<const kdtree_node *> work_node_vec;
    std::vector<correlate> corr_obj_vec;
private:                        // Function
    void get_node_vec( const kdtree_node * root );
    void add_work_node( const kdtree_node * node,
                        int depth_remain );
public:    
    void set_num_threads( int num_threads_src );
    void set_dist_bin( double s_max, double s_min,
                       int num_bins );

    ////////// Conduct calculation //////////
public:
    void cal_corr( const kdtree & tree0,
                   const kdtree & tree1 );
    void output( std::string file_name );
};

#endif
