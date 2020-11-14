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
 * along with Susa. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file vector.h
 * @brief stream overload for STL containers (declaration and definition).
 *
 * @author Behrooz Kamary
 * @version 1.0.0
 *
 * @ingroup Math
 */

#include <iterator>
#include <vector>

template <typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &vec_arg)
{
    if (!vec_arg.empty())
    {
        out << '[';
        std::copy(vec_arg.begin(), vec_arg.end(), std::ostream_iterator<T>(out, ", "));
        out << "\b\b]";
    }
    return out;
}