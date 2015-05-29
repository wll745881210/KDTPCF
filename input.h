#ifndef INPUT_H_
#define INPUT_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>    
#include <vector>
#include <unordered_map>

class input
{
public:
    input(  );
    input( const std::string & file_name );
    ~input(  );

    void read(  );
    template <typename T, typename t>
    void find_key( const std::string key_name,
	           T & val, t def_val );

private:
    std::ifstream fin;
    std::unordered_map<std::string, std::string> item_map;
    
    void get_items(  );
};

template <typename T, typename t>
void input::find_key( const std::string key_name,
                      T & val, t def_val )
{
    auto p = item_map.find( key_name );
    if( p != item_map.end(  ) )
    {
	std::stringstream ss;
	ss.str( p->second );
	ss >> val;
    }
    else
    {
	std::cout << "Entry \"" << key_name
		  << "\" not found; Using default value: "
		  << def_val << std::endl;
	val = def_val;
    }
    return;
}

#endif
