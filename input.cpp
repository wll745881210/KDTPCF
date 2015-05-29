#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include "input.h"

input::input(  ) : fin( "par.txt" )
{

}

input::input( const std::string & file_name )
    : fin( file_name.c_str(  ) )
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
    }while( item_temp[ 0 ] == '#' );
    std::cout << '\n';
    
    while( !fin.eof(  ) )
    {
        getline( fin, item_temp );
        std::cout << item_temp << std::endl;
        if( item_temp == "Init" )
            break;
    }
    get_items(  );
    return ;
}

void input::get_items(  )
{
    std::string item_temp, line_temp, value_temp;
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
	std::pair<std::string, std::string>
	    item ( item_temp, value_temp );
	item_map.insert( item );
	
        std::cout << std::setw( 20 ) << std::left
                  << item_temp << value_temp << std::endl;
    }
    return;
}
