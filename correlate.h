// Copyright (C) 2013 Lile Wang
// For full notice, see "main.cpp" and "COPYING".

#ifndef CORRELATE_H_
#define CORRELATE_H_

#include "kdtree.h"
#include <fstream>
#include <vector>
#include <string>

class correlate
{
    ////////// Con/destructor and initializer //////////
private:
    static bool is_auto_cor, is_2d_cor, is_ang_cor;
    typedef 
public:
    correlate(  );
    ~correlate(  );
    void clear(  );
    static void static_clear(  );
    static void set_auto_cor( bool auto_cor );
    static void set_par
    ( double s_max_s, double s_min_s, int s_num_s,
      int phi_num_s, bool log_bin, int corr_stat, int jk_n );
    
    ////////// Compare trees //////////
private:                        // Data
    std::vector<unsigned> bin_counts;
private:                        // Functions
    void brute_force( const kdtree_node * node0,
                      const kdtree_node * node1 );
    int node_bin( const kdtree_node * node0,
                  const kdtree_node * node1 );
public:
    void compare_node( const kdtree_node * node0,
                       const kdtree_node * node1 );

    ////////// Jackknife //////////
private:                        // Data
    static int jk_num;
    std::vector<unsigned> bin_counts_jk;
private:                        // Function
    inline void jk_add( int idx, int sample, int add );
    inline void jk_add( int idx, int sample );

    ////////// Output //////////
private:                        // Data
    static std::vector<unsigned> bin_counts_tot;
    static std::vector<unsigned> bin_counts_tot_jk;
private:                        // Function
    static void out_one_line( std::ofstream & fout, int idx );
public:
    void add_to_tot(  );
    static void output( std::string file_name );
    static const std::vector<unsigned> & bin_count_ref(  );
    static const std::vector<unsigned> & bin_count_jk_ref(  );
    static void get_bin_center( std::vector<double> & s_v,
                                std::vector<double> & phi_v );
    
    ////////// 2D distance //////////
private:                        // Function
    inline static int idx_2d( int i, int j );
    inline static double dot( const double a[  ],
                              const double b[  ] );

    ////////// Binning values //////////
private:                        // Data
    static double s_max, s_min, ds;
    static double * s_bin_lim, * phi_bin_lim;
    static int s_num, phi_num;
    static bool is_log_bin;
private:                        // Function
    static double s_lim_val( int i );
    static double phi_lim_val( int i );
    inline static int s_idx_val( const double & a_2 );
    inline static int s_idx_arr( const double a[  ] );
    inline static int phi_idx_val( const double & mu );
public:
    static double s_center( int i );
    static double phi_center( int i );

    ////////// Constants //////////
private:
    static const double pi_2 = 1.5707963267948966;
    static const double rad_to_deg = 57.29577951308232;
    static const int dim = 3;
    static const double tiny = 1.e-2;
};

#endif

