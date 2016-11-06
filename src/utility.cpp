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
 * @file utility.cpp
 * @brief The utility functions.
 * @author Behrooz Aliabadi
 * @version 1.0.0
 */

#include <susa.h>
#include <chrono>

namespace susa {

// Elapsed time utility functions
static std::chrono::time_point<std::chrono::high_resolution_clock> __susa_start;
static std::chrono::time_point<std::chrono::high_resolution_clock> __susa_end;


void tic()
{
    __susa_start = std::chrono::high_resolution_clock::now();
}

void toc_print()
{
    __susa_end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(__susa_end - __susa_start);
    std::cout << "Elapsed time is " << elapsed.count()/1000.0f << " seconds." << std::endl;
}

double toc() {
    __susa_end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(__susa_end - __susa_start);
    return((double)elapsed.count()/1000.0f);
}

std::string timestamp()
{

    time_t rawtime;
    std::stringstream ss_itos;
    std::stringstream ss_ret;
    tm *ptm;


    time(&rawtime);
    ptm = localtime(&rawtime);
    ss_ret << std::string("20");

    ss_itos.str(std::string());
    ss_itos << ptm->tm_year - 100;
    if (ss_itos.str().length() > 1) ss_ret << ss_itos.str();
    else ss_ret << std::string("0") << ss_itos.str();

    ss_itos.str(std::string());
    ss_itos << ptm->tm_mon+1;
    if (ss_itos.str().length() > 1) ss_ret << ss_itos.str();
    else ss_ret << std::string("0") << ss_itos.str();

    ss_itos.str(std::string());
    ss_itos << ptm->tm_mday;
    if (ss_itos.str().length() > 1) ss_ret << ss_itos.str();
    else ss_ret << std::string("0") << ss_itos.str();

    ss_itos.str(std::string());
    ss_itos << ptm->tm_hour;
    if (ss_itos.str().length() > 1) ss_ret << ss_itos.str();
    else ss_ret << std::string("0") << ss_itos.str();

    ss_itos.str(std::string());
    ss_itos << ptm->tm_min;
    if (ss_itos.str().length() > 1) ss_ret << ss_itos.str();
    else ss_ret << std::string("0") << ss_itos.str();

    ss_itos.str(std::string());
    ss_itos << ptm->tm_sec;
    if (ss_itos.str().length() > 1) ss_ret << ss_itos.str();
    else ss_ret << std::string("0") << ss_itos.str();

    return ss_ret.str();
}


} // NAMESPACE SUSA
