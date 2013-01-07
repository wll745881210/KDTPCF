// Copyright (C) 2013 Lile Wang
// For full notice, see "main.cpp" and "COPYING".

#include "read_data.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

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
	std::cout << "Reading from file \"" << file_name
			  << "\"...";
	std::cout.flush(  );
	
	galaxy_point temp;
	std::vector<galaxy_point> & buf = tree.source_list;
	buf.clear(  );
	fin >> temp.x >> temp.y >> temp.z;
	fin.ignore( 64, '\n' );
	while( ! fin.eof(  ) )
	{
		convert( temp );
		buf.push_back( temp );
		fin >> temp.x >> temp.y >> temp.z;
		fin.ignore( 64, '\n' );
	}
	std::cout << " Done." << std::endl;
	fin.close(  );
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

void read_data::convert( galaxy_point & src )
{
	const double ra  = src.x / rad_to_deg;
	const double dec = src.y / rad_to_deg;
	const double rsh = src.z;
	const double chi = chi_of_z( rsh );
	src.x = cos( dec ) * cos( ra ) * chi;
	src.y = cos( dec ) * sin( ra ) * chi;
	src.z = sin( dec ) * chi;
	return;
}

