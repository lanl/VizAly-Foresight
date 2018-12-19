/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
================================================================================*/
#pragma once

#include<stdio.h> 
#include<stdbool.h> 
#include <string>
#include <fstream>

inline bool fileExisits(char *filename) 
{
    std::ifstream ifs(filename);
    return ifs.good();
}



bool isPowerOfTwo(int n) 
{ 
  	if (n == 0) 
    	return 0; 

  	while (n != 1) 
  	{ 
      	if (n%2 != 0) 
         	return 0; 
      	n = n/2; 
  	} 
  	return 1; 
} 


std::string extractFileName(std::string inputString)
{
	std::size_t pos = inputString.find_last_of("/\\");
	return inputString.substr(pos+1); 
}