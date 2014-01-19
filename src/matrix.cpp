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
 * @file matrix.cpp
 * @brief matrix class.
 * @author Behrooz, Kamary Aliabadi
 * @version 1.0.0
 */


#include "../inc/susa.h"

namespace susa {

void pre_parser(std::string &str_string) {
  // Coversion from 'string' to numerical data types are done
  // by using the capabilities of STL stream objects.
  // They faciliate parsing the string into a Susa matrix.
  // This pre-parsing method is to prepare the input string for
  // STL stream objects in parser() methode.
  // At the end of this routine the string should only contain
  // neccessary spaces and semicolons.
  int int_length = str_string.length();
  
  // it is needed to append some extra spaces since in the loops below the
  // code will access elements outside the original string size them.
  str_string = str_string + std::string("  ");


  // (0) Replace tabs and "," with spaces
  for (int int_i = 0; int_i < int_length; int_i++) if (str_string[int_i] == 0x09 || str_string[int_i] == ',') str_string[int_i] = 0x20;

  // (1) Remove spaces (0x20) after ';'
  int int_counter = 0;
  for (int int_i = 0; int_i < int_length; int_i++) {
    str_string[int_counter] = str_string[int_i];
    if (str_string[int_i] == ';' && str_string[int_i+1] == ' ') {
      while(str_string[int_i+1] == ' ') int_i++;
    }
    int_counter++;
  }
  // Update the string according to its new shape/size
  str_string.erase(int_counter,int_length - 1);
  int_length = str_string.length();
  str_string = str_string + std::string("  ");
  

  // (2) This is to remove extra spaces.
  int_counter = 0;
  for (int int_i = 0; int_i < int_length; int_i++) {
    if (str_string[int_i] != 0x20 && str_string[int_i] != ']' && str_string[int_i] != '[') {
      str_string[int_counter] = str_string[int_i];
      int_counter++;
    } else if (str_string[int_i] == 0x20 && str_string[int_i + 1] != ';' && ((str_string[int_i + 1] > 0x2F || str_string[int_i + 1] == '.'))) {
		if (int_counter != 0) { // To avoid the 0x20 in the begining of the string
		  str_string[int_counter] = 0x20;
          int_counter++;
        }		  
      }
  }
    
    // (3) shrink the string to its real size (delete extra 0x20 that was created in the loop above)
    if (str_string[int_counter-1] == 0x20) str_string.erase(int_counter-1,int_length - 1);
    else str_string.erase(int_counter,int_length - 1);

    //std::cout << std::endl << str_string;

}
}      // NAMESPACE SUSA
