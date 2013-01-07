// Copyright (C) 2013 Lile Wang
// For full notice, see "main.cpp" and "COPYING".

#ifndef INPUT_H_
#define INPUT_H_

#include <fstream>
#include <string>
#include <vector>

class input
{
public:
	input(  );
	input( const std::string & file_name );
	~input(  );

	void read(  );
	void get_init(  );

	double lambda;
	double z_max;
	double s_max;
	double s_min;
	double s_bin;
	std::string data_file_name;
	std::string rand_file_name;
	
private:
	std::ifstream fin;
	std::vector<std::string> item_name;
	std::vector<std::string> value;
	int	length;
	
	void sort_item(  );
	void get_general(  );
};

#endif

