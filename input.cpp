// Copyright (C) 2013 Lile Wang
// For full notice, see "main.cpp" and "COPYING".

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include "input.h"

input::input(  ) : fin( "par.txt" ), length( 0 )
{

}

input::input( const std::string & file_name )
	: fin( file_name.c_str(  ) ), length( 0 )
{

}

input::~input(  )
{
	
}

void input::read(  )
{
	std::string item_temp;
	if( !fin )
		throw( "Unable to open input." );
	getline( fin, item_temp );
	if( item_temp != "CORR_FUNC" )
		throw( "Incorrect format of input file." );
	std::cout << item_temp << std::endl;
	
	do{
		getline( fin, item_temp );
		std::cout << item_temp;
	}while( item_temp[ 0 ] == '/' );
	std::cout << '\n';
	
	while( !fin.eof(  ) )
	{
		getline( fin, item_temp );
		std::cout << item_temp << std::endl;
		if( item_temp == "Init" )
			break;
	}
	
	return ;
}

void input::sort_item(  )
{
	std::string item_temp;
	std::string	value_temp;

	int i( 0 ), j( 0 );

	for( i = 1; i < length; ++ i )
	{
		item_temp  = item_name[ i ];
		value_temp = value[ i ];
		for( j = i; j > 0 && item_temp <
				 item_name[ j-1 ]; -- j )
		{
			item_name[ j ] = item_name[ j-1 ];
			value[ j ]	   = value[ j-1 ];
		}
		item_name[ j ] = item_temp;
		value[ j ]	   = value_temp;
	}
	
	return;
}

void input::get_general(  )
{
	std::string item_temp, line_temp;
	std::string value_temp;
	std::stringstream ss;
	
	while( !fin.eof(  ) )
	{
		getline( fin, line_temp );
		
		if( line_temp == "" || line_temp[ 0 ] == ' '
			|| line_temp[ 0 ] == '\t' )
			break;
		
		ss.clear(  );
		ss.str( line_temp );		
		ss >> item_temp;
		ss >> value_temp;
		item_name.push_back( item_temp );
		value.push_back( value_temp );
		std::cout << std::setw( 15 ) << std::left
				  << item_temp << value_temp << std::endl;
	}
	length = item_name.size(  );
	sort_item(  );	
	return;
}

void input::get_init(  )
{
	get_general(  );
	int i( 0 );
	std::stringstream ss;

	data_file_name = value[ i ];
	++ i;
	rand_file_name = value[ i ];
	++ i;

	ss.clear(  );
	ss.str( value[ i ] );
	ss >> lambda;
	++ i;

	ss.clear(  );
	ss.str( value[ i ] );
	ss >> s_bin;
	++ i;

	ss.clear(  );
	ss.str( value[ i ] );
	ss >> s_max;
	++ i;

	ss.clear(  );
	ss.str( value[ i ] );
	ss >> s_min;
	++ i;

	ss.clear(  );
	ss.str( value[ i ] );
	ss >> z_max;
	++ i;
	
	return;
}

