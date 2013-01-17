// Copyright (C) 2013 Lile Wang
// For full notice, see "main.cpp" and "COPYING".

#ifndef INPUT_H_
#define INPUT_H_

#include <fstream>
#include <sstream>
#include <string>    
#include <vector>

class input
{
public:
    input(  );
    input( const std::string & file_name );
    ~input(  );

    void read(  );
    template <typename T>
    void find_key( std::string key_name, T & val );

private:
    std::ifstream fin;
    std::vector<std::string> item_name;
    std::vector<std::string> value;
    int    length;
    
    void get_items(  );
};

// I have to implement this function here since it
// involves template...

template <typename T>
void input::find_key( std::string key_name, T & val )
{
    std::stringstream ss;
    for( unsigned i = 0; i < item_name.size(  ); ++ i )
        if( item_name[ i ].compare( key_name ) == 0 )
        {
            ss.str( value[ i ] );
            ss >> val;
            return;
        }
    ss.str( "0" );
    ss >> val;
    return;
}

#endif

