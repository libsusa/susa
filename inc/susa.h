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
 * @brief Package Header
 * @author Behrooz Kamary Aliabadi
 * @version 1.0.0
 */

/**
 * @mainpage Susa Documentation
 *
 *
 * @section intro Introduction
 * <p> Susa is a library of utility classes and methods (functions) for signal processing and mathematics applications.
 * It is portable and stand-alone, therefore it does not need complicated and time consuming installation of third
 * party libraries. Susa also provides a framework for embedded devices. You may use it to develop quick and efficient
 * signal processing methods for mobile platforms such as Android using the companion Native Development Kit (NDK).
 * It is published under GNU Lesser General Public License.</p>
 *
 */

#ifndef SUSA_H
#define SUSA_H

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

// Constants
#define PI  3.1415926535897932384626433

#include "memory.h"
#include "sets.h"
#include "matrix.h"
#include "array.h"
#include "base.h"
#include "svd.h"
#include "statistics.h"
#include "rng.h"
#include "mt.h"
#include "fft.h"
#include "rrcosine.h"
#include "channel.h"
#include "convolutional.h"
#include "modulation.h"
#include "utility.h"
#include "linalg.h"
#include "search.h"
#include "signal.h"
#include "debug.h"

#endif // SUSA_H
