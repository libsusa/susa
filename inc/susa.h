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
 * @file susa.h
 * @brief Susa Package Header
 * @author Behrooz Kamary
 * @version 1.0.0
 */

/**
 * @mainpage Susa Documentation
 *
 *
 * @section intro Introduction
 * <p> Susa is a library of utility classes (types) and routines (functions) for signal processing and mathematics applications.
 * It is portable and stand-alone, therefore it does not need complicated and time consuming installation of third
 * party libraries. Susa also provides a framework for embedded devices. You may use it to develop quick and efficient
 * signal processing software for mobile platforms such as Android using the companion Native Development Kit (NDK).
 * It is published under GNU Lesser General Public License.</p>
 *
 */

#ifndef SUSA_H
#define SUSA_H

// Constants
#define PI  3.1415926535897932384626433

// STL headers
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include <climits>
#include <ctime>
#include <cstdlib>
#include <string>
#include <cstring>
#include <initializer_list>
#include <typeinfo>
#include <tuple>
#include <type_traits>

// Susa headers
#include "susa/type_traits.h"
#include "susa/memory.h"
#include "susa/sets.h"
#include "susa/matrix.h"
#include "susa/array.h"
#include "susa/base.h"
#include "susa/svd.h"
#include "susa/statistics.h"
#include "susa/rng.h"
#include "susa/mt.h"
#include "susa/fft.h"
#include "susa/rrcosine.h"
#include "susa/channel.h"
#include "susa/ccode.h"
#include "susa/modulation.h"
#include "susa/utility.h"
#include "susa/linalg.h"
#include "susa/solver.h"
#include "susa/search.h"
#include "susa/signal.h"
#include "susa/debug.h"

#endif // SUSA_H
