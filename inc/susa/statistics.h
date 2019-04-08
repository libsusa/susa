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
 * @file statistics.h
 * @brief Basic statistical operations on STL and Susa types
 * @author Behrooz, Kamary
 * @version 1.0.0
 *
 * @defgroup Statistics Statistics
 */

#ifndef STATISTICS_H
#define STATISTICS_H

namespace susa {

/**
* @brief Mean of a STL vector
*
* @param vec_arg Input STL vector
* @ingroup Statistics
*/
template <class T> T mean(std::vector <T> vec_arg);

/**
* @brief Variance of a STL vector
*
* @param vec_arg Input STL vector
* @ingroup Statistics
*/
template <class T> double var(std::vector <T> vec_arg);

/**
* @brief Standard deviation of a STL vector
*
* @param vec_arg Input STL vector
* @ingroup Statistics
*/
template <class T> double stdv(std::vector <T> vec_arg);


template <class T> T mean(std::vector <T> vec_arg) {
    T tMean=0;
    for (int i=0; i<vec_arg.size(); i++) {
        tMean += vec_arg[i];
    }

    return (tMean/vec_arg.size());
}



template <class T> double var(std::vector <T> vec_arg) {
    T tMean = mean(vec_arg);
    T tVar = 0;
    for (int i=0; i<vec_arg.size(); i++) {
        tVar += (vec_arg[i] - tMean)*(vec_arg[i] - tMean);
    }

    return (tVar/vec_arg.size());
}




template <class T> double stdv(std::vector <T> vec_arg) {
    T tMean = mean(vec_arg);
    T tStd = 0;
    for (int i=0; i<vec_arg.size(); i++) {
        tStd += (vec_arg[i] - tMean)*(vec_arg[i] - tMean);
    }

    return (sqrt(tStd/vec_arg.size()));
}
} // NAMESPACE SUSA
#endif // STATISTICS_H
