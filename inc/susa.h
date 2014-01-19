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
 * @author Behrooz, Kamary Aliabadi
 * @version 1.0.0
 */

/**
 * @mainpage Susa Documentation
 *
 * @author Behrooz, Kamary Aliabadi
 *
 * @section intro Introduction
 * <p> Susa is a library of utility classes and functions to be used in signal processing and mathematic applications. 
 * This library package is portable and light. It does not need complicated and time consuming installation of other 
 * mathematical libraries to work with. Because of this independency at the moment it's limited basic operations. 
 * Developers of Susa have a long term plan to improve the library. It is published publicly under Lesser GPL license.</p>
 * 
 */

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <complex>
#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include <climits>
#include <ctime>
#include <cstdlib>
#include <string>
#include <cstring>

#ifndef SUSA_H
#define SUSA_H

// Constants
#define PI  3.1415926535897932384626433

// Compile Switches
#define _SUSA_TERMINATE_ON_ERROR  // If defined the program will be terminated by 'exit(1)' when error occurs.

#include "matrix.h"
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

#endif // SUSA_H




