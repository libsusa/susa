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
 * @brief The <i>set</i> type (declaration and definition).
 * @author Behrooz Kamary
 * @version 1.0.0
 *
 */

#ifndef SUSA_SETS_H
#define SUSA_SETS_H

#include <susa/debug.h>
#include <susa/memory.h>

namespace susa
{

/**
* @brief The <i>bitset</i> class.
* <i>bitset</i> is a set type aimed to store
* a set of bits in a quick (using single memory allocation)
* and memory efficient manner. You may alternatively use <i>std::bitset</i>
* where the set size is in hand at compile time (no runtime initialization).
* @ingroup TYPES
*
*/
class bitset :
public memory <unsigned char>
{

  public :

    bitset(size_t maxsize);

    void set(size_t index);

    void reset(size_t index);

    size_t pop();

    void push(size_t index);

    bool exists(size_t index);

    void set();

    void reset();

    bool any();

  private:

    size_t nbytes;
    size_t uint_set_size;

};

}       // NAMESPACE SUSA
#endif  // SUSA_SETS_H
