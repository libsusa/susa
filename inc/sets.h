/*
 * This file is part of Susa.
 *
 * Susa is free software: you can redistribute it and/or modify
 * it under the terms of the Lesser GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * at your option) any later version.
 *
 * Susa is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * Lesser GNU General Public License for more details.
 *
 * You should have received a copy of the Lesser GNU General Public License
 * along with Susa.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file sets.h
 * @brief Set types and operations.
 * @author Behrooz Kamary Aliabadi
 * @version 1.0.0
 *
 * @defgroup Set Set
 */

#ifndef SUSA_SETS_H
#define SUSA_SETS_H

#include <debug.h>

namespace susa
{

class index_set
{
  public :
    index_set(unsigned int uint_num_nodes);

    ~index_set();

    void add(unsigned int index);

    void remove(unsigned int index);
    
    unsigned int pop();

    void push(unsigned int index);

    bool exists(unsigned int index);
    
    void add_all();
    
    void remove_all();
    
    bool is_not_empty();

  private:
    unsigned int nbytes;
    unsigned int uint_set_size;
    char* entity;
};

}       // NAMESPACE SUSA
#endif  // SUSA_SETS_H
