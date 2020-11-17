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
 * @file auxiliary.h
 * @brief auxiliary (declaration and definition).
 *
 * @author Behrooz Kamary
 *
 */

#ifndef SUSA_AUXILIARY_H
#define SUSA_AUXILIARY_H

#include <iterator>
#include <vector>
#include <tuple>
#include <complex>

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace aux
{
    template <std::size_t...>
    struct seq
    {
    };

    template <std::size_t N, std::size_t... Is>
    struct gen_seq : gen_seq<N - 1, N - 1, Is...>
    {
    };

    template <std::size_t... Is>
    struct gen_seq<0, Is...> : seq<Is...>
    {
    };

    template <class Ch, class Tr, class Tuple, std::size_t... Is>
    void print_tuple(std::basic_ostream<Ch, Tr> &out, Tuple const &tuple_arg, seq<Is...>)
    {
        using swallow = int[];
        (void)swallow{0, (void(out << (Is == 0 ? "" : ", ") << std::get<Is>(tuple_arg)), 0)...};
    }
} // namespace aux

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

template <typename T, typename S>
std::ostream &operator<<(std::ostream &out, const std::pair<T, S> &pair_arg)
{
    return out << "(" << pair_arg.first << ", " << pair_arg.second << ")";
}


template <class Ch, class Tr, class... Args>
auto operator<<(std::basic_ostream<Ch, Tr> &out, std::tuple<Args...> const &tuple_arg) -> std::basic_ostream<Ch, Tr>&
{
    out << "(";
    aux::print_tuple(out, tuple_arg, aux::gen_seq<sizeof...(Args)>());
    return out << ")";
}

#endif
#endif