// Copyright (C) 2013 Lile Wang
// For full notice, see "main.cpp" and "COPYING".

#ifndef DRIVER_H_
#define DRIVER_H_

#include "kdtree.h"
#include "correlate.h"
#include "read_data.h"
#include "input.h"
#include "parallel.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <ctime>
#include <omp.h>
#include <cmath>

class driver
{
    ////////// Con/destructor & initilaizer //////////
public:
    driver( const std::string & par_file );
    ~driver(  );

    ////////// Read from par file //////////
private:                        // Data
    std::string par_file_name;
    int         corr_stat;
    int         num_threads;
    int         bin_count_type;
    double      s_max, s_min;
    int         s_num, phi_num;
    int         log_bin;
    int         jk_depth, jk_num;
    std::string data_file_name;
    std::string rand_file_name;
    std::string out_name_base;
    double      lambda, z_max;
private:                        // Function
    void read_from_par(  );

    ////////// Build the trees //////////
private:
    void build_trees(  );

    ////////// Calculate //////////
private:                        // Data
    kdtree data,      rand;
    int    data_size, rand_size;
    std::vector<long long unsigned> dd, dd_jk;
    std::vector<long long unsigned> rr, rr_jk;
    std::vector<long long unsigned> dr, dr_jk;
private:                        // Function
    void cal(  );
    void cal_ls(  );			// Landy-Sazlay

    ////////// Timer //////////
private:                        // Data
    double start_t, precomp_t, auto_t, cross_t;
private:                        // Function
    void show_wall_t( std::string title,
                      double & start, double & end );

    ////////// Get everything! //////////
public:
    void go(  );
};

#endif

