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
 * @file utility.h
 * @brief Utility functions.
 * @author Behrooz, Aliabadi
 * @version 1.0.0
 * @defgroup Utility
 *
 */

#ifndef UTILITY_H
#define UTILITY_H

namespace susa {

/**
 * @brief Starts the internal chronometer.
 *
 *
 * @ingroup Utility
 */
void tic();

/**
 * @brief Prints out the internal chronometer.
 *
 *
 * @ingroup Utility
 */
void toc_print();

/**
 * @brief Reset the internal chronometer.
 *
 * @return The elapsed time in seconds.
 *
 * @ingroup Utility
 */
double toc();

/**
 * @brief Gives a time-stamp in string format.
 *
 * @return A string that contains year/month/day/hour/minutes/seconds.
 *
 * @ingroup Utility
 */
std::string timestamp();

} // NAMESPACE SUSA

#endif

