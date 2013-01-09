// KDTPCF: Two-point statistics using k-d tree.
// Copyright (C) 2013 Lile Wang; lilewang@lbl.gov

// This program is free software: you can redistribute it
// and/or modify it under the terms of the GNU General
// Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your
// option) any later version.

// This program is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.

#include "driver.h"
#include <iostream>

int main( int argn, char * argv[  ] )
{
	try
	{
		std::string par_file_name;
		if( argn == 1 )
			par_file_name = "par.txt";
		else if( argn == 2 )
			par_file_name = argv[ 1 ];
		else
			throw "Incorrect parameter file.";

		driver( par_file_name );
	}
	catch( const char * err )
	{
		std::cerr << "\nError: " << err << std::endl;
	}
	
	return 0;
}


