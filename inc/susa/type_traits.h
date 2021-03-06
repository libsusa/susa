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
 * @file type_traits.h
 * @brief Susa type traits (declaration).
 * @author Behrooz Kamary
 *
 */

#ifndef SUSA_TYPE_TRAITS_H
#define SUSA_TYPE_TRAITS_H

#include <type_traits>

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace susa
{

template<class T> struct is_complex : public std::false_type {};

template<class T> struct is_complex<std::complex<T>> : public std::true_type {};

template<typename, typename>
struct get_allocator {};

template<typename A, template<typename> typename C, typename B>
struct get_allocator<A, C<B>> {
    using type = C<A>;
};

}
#endif // DOXYGEN_SHOULD_SKIP_THIS
#endif // SUSA_TYPE_TRAITS_H
