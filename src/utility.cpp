/*
 * This file is part of Susa.

 * Susa is free software: you can redistribute it and/or modify
 * it under the terms of the Lesser GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * at your option) any later version.

 * Susa is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * Lesser GNU General Public License for more details.

 * You should have received a copy of the Lesser GNU General Public License
 * along with Susa.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file utility.cpp
 * @brief The utility functions.
 * @author Behrooz, Kamary Aliabadi
 * @version 1.0.0
 */

#include "../inc/susa.h"

namespace susa {

// Elapsed time utility functions
unsigned int ___uint_stamp;
void tic() { ___uint_stamp = clock(); }
void toc_print() {std::cout << std::endl << "Elapsed time is " << ((double)(clock() - ___uint_stamp)/CLOCKS_PER_SEC) << " seconds.";}
double toc() {return((double)(clock() - ___uint_stamp)/CLOCKS_PER_SEC);}

std::string timestamp() {

  std::string str_ret("20");
  time_t rawtime;
  char char_buffer[5];
  tm *ptm;


  time(&rawtime);
  ptm = localtime(&rawtime);

  itoa(ptm->tm_year-100,char_buffer,10);
  str_ret += std::strlen(char_buffer) > 1 ? std::string(char_buffer) : std::string("0") + std::string(char_buffer);

  itoa(ptm->tm_mon+1,char_buffer,10);
  str_ret += std::strlen(char_buffer) > 1 ? std::string(char_buffer) : std::string("0") + std::string(char_buffer);

  itoa(ptm->tm_mday,char_buffer,10);
  str_ret += std::strlen(char_buffer) > 1 ? std::string(char_buffer) : std::string("0") + std::string(char_buffer);

  itoa(ptm->tm_hour,char_buffer,10);
  str_ret += std::strlen(char_buffer) > 1 ? std::string(char_buffer) : std::string("0") + std::string(char_buffer);

  itoa(ptm->tm_min,char_buffer,10);
  str_ret += std::strlen(char_buffer) > 1 ? std::string(char_buffer) : std::string("0") + std::string(char_buffer);

  itoa(ptm->tm_sec,char_buffer,10);
  str_ret += std::strlen(char_buffer) > 1 ? std::string(char_buffer) : std::string("0") + std::string(char_buffer);
  
  return str_ret;
}

// Temporary Utility for Utility
// Taken from a web resource.
	
void strreverse(char* begin, char* end) {
	char aux;
	while(end>begin) aux=*end, *end--=*begin, *begin++=aux;
}
	
void itoa(int value, char* str, int base) {
  static char num[] = "0123456789abcdefghijklmnopqrstuvwxyz";
  char* wstr=str;
  int sign;
  div_t res;

  // Validate base

  if (base<2 || base>35) { 
    *wstr='\0';
	return;
  }

  // Take care of sign
  sign = value;
  if (value < 0) value = -value;
	
  // Conversion. Number is reversed.
  do {
    res = div(value,base);
    *wstr++ = num[res.rem];
    value=res.quot;
	} while(value);
	if(sign<0) *wstr++='-';
	*wstr='\0';	
	
	// Reverse string
	strreverse(str,wstr-1);
}


} // NAMESPACE SUSA
