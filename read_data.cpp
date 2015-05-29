// Copyright (C) 2013 Lile Wang
// For full notice, see "main.cpp" and "COPYING".

#include "read_data.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

////////////////////////////////////////////////////////////
// Static variables

const double read_data::c	    = 299792.458;
const double read_data::H1	    = 100;  // in km/s/Mpc
const double read_data::pi	    = 3.14159265358979;
const double read_data::nearly_zero = 1e-15;
const double read_data::rad_to_deg  = 57.29577951308232;


////////////////////////////////////////////////////////////
// Constructor, destructor and initializer

read_data::read_data(  )
{
    
}

read_data::~read_data(  )
{
    
}

void read_data::set_cosmology( const double lambda,
                               const double z_max )
{
    this->lambda = lambda;
    this->z_max  = z_max;
    dz = 0.05;
    get_chi_of_z(  );
    return;
}

////////////////////////////////////////////////////////////
// Read from files

void read_data::read_from_file( std::string file_name,
                                kdtree & tree )
{
    std::ifstream fin( file_name.c_str(  ) );
    if( !fin )
        throw "Unable to open source list.";
    
    galaxy_point temp;
    temp.weight = 1.;
    std::vector<galaxy_point> & buf = tree.source_list;
    buf.clear(  );
    while( true )
    {
        fin >> temp.x[ 0 ] >> temp.x[ 1 ];
        if( !is_ang_cor )
	    fin >> temp.x[ 2 ];
	if( is_weighted )
	    fin >> temp.weight;
        fin.ignore( 64, '\n' );
	if( fin.eof(  ) )
	    break;
        convert( temp );
        buf.push_back( temp );
    }
    return;
}

////////////////////////////////////////////////////////////
// Set cosmology

double read_data::dchi_dz( double z )
{
    return c / H1 / sqrt( ( 1 - lambda ) * pow( 1 + z, 3 )
                          + lambda );
}

void read_data::rk4( double & chi, double & z )
{
    double chi0[ 4 ];
    chi0[ 0 ] = dchi_dz( z );
    chi0[ 1 ] = dchi_dz( z + dz / 2 );
    chi0[ 2 ] = chi0[ 1 ];
    chi0[ 3 ] = dchi_dz( z + dz );

    z   += dz;
    chi += dz / 6.
        * ( chi0[ 0 ] + 4. * chi0[ 1 ] + chi0[ 2 ] );
    return;
}

void read_data::get_chi_of_z(  )
{
    double z( 0. ), chi( 0. );

    while( z < z_max )
    {
        chi_of_z_buf.push_back( chi );
        rk4( chi, z );
    }
    return;
}

double read_data::chi_of_z( const double & z )
{
    unsigned idx = unsigned( z / dz );
    if( idx + 1 > chi_of_z_buf.size(  ) )
    {
        std::cerr << z << ' ' << idx << std::endl;
        throw "Redshift z value exceeds boundary.";
    }
    
    const double chi0 = chi_of_z_buf[ idx ];
    const double chi1 = chi_of_z_buf[ idx + 1 ];
    return chi0 + ( chi1 - chi0 ) * ( z - idx * dz ) / dz;
}

////////////////////////////////////////////////////////////
// Coordinate conversion

void read_data::convert( galaxy_point & src )
{
    const double ra  = src.x[ 0 ] / rad_to_deg;
    const double dec = src.x[ 1 ] / rad_to_deg;
    const double rsh = src.x[ 2 ];
    const double chi = ( is_ang_cor ? 1. : chi_of_z( rsh ) );
    src.x[ 0 ] = cos( dec ) * cos( ra ) * chi;
    src.x[ 1 ] = cos( dec ) * sin( ra ) * chi;
    src.x[ 2 ] = sin( dec ) * chi;
    return;
}

void read_data::set_par( bool ang_cor, bool weighted )
{
    this->is_ang_cor  = ang_cor;
    this->is_weighted = weighted;
    return;
}
