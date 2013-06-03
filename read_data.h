// Copyright (C) 2013 Lile Wang
// For full notice, see "main.cpp" and "COPYING".

#ifndef READ_DATA_H_
#define READ_DATA_H_

#include "kdtree.h"
#include <string>
#include <vector>


class read_data
{
    ////////// Con/destructor and initializer //////////
public:
    read_data(  );
    ~read_data(  );
    void set_cosmology( const double lambda,
                        const double z_max );
    
    ////////// Read from files //////////
public:
    void read_from_file( std::string file_name,
                         kdtree & tree );

    ////////// Cosmology //////////
private:                        // Data
    double lambda;
    double z_max;
    double dz;
    std::vector<double> chi_of_z_buf;
private:                        // Function
    double dchi_dz( double z );
    void rk4( double & chi, double & z );
    void get_chi_of_z(  );
    double chi_of_z( const double & z );
    void convert( galaxy_point & src );

    ////////// Specification //////////
private:                        // Data
    bool is_ang_cor, is_weighted;
public:
    void set_par( bool ang_cor, bool weighted );

    ////////// Constants //////////
private:
    static const double c  = 299792.458; // c in km/s
    static const double H1 = 100;         // in km/s/Mpc
    static const double pi = 3.14159265358979;
    static const double nearly_zero = 1e-15;
    static const double rad_to_deg = 57.29577951308232;
};

#endif

