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
 * @file debug.h
 * @brief A separate assertion method to catch the fatal errors.
 *
 * This may be useful when one needs to desactivate debugging
 * of Susa independent of STL assertion by defining SUSA_NDEBUG.
 * If NDEBUG is defined, it consequently affects this assertion.
 *
 * @author Behrooz Kamary
 */

#ifndef SUSA_DEBUG_H
#define SUSA_DEBUG_H

#include <iostream>
#include <cassert>
#include <cstdlib>

namespace susa
{

  inline void assert_log(const char* file, const char* func, int line, const char* cond, const char* msg)
  {
    std::cerr << "Assertion failed: " << "(" << cond <<")";
    std::cerr << ", in function " << func;
    std::cerr << ", file " << file;
    std::cerr << ", line " << line;
    std::cerr << ", message " << msg;
    std::cerr << std::endl;
    std::abort();
  }

} // NAMESAPCE SUSA

#ifdef SUSA_NDEBUG
#define SUSA_ABORT(MSG) ((void)0)
#define SUSA_ASSERT(EX) ((void)0)
#define SUSA_ASSERT_MESSAGE(EX,MSG) ((void)0)
#else
#define SUSA_ABORT(MSG) (susa::assert_log(__FILE__, __func__, __LINE__, "abort", #MSG))
#define SUSA_ASSERT(EX) ((EX) ? (void)0 : assert(EX))
#define SUSA_ASSERT_MESSAGE(EX,MSG) ((EX) ? (void)0 : susa::assert_log(__FILE__, __func__, __LINE__, #EX, #MSG))
#endif

#ifdef SUSA_NDEBUG
#define SUSA_LOG_DBG(MSG) ((void)0)
#define SUSA_LOG_INF(MSG) ((void)0)
#define SUSA_LOG_ERR(MSG) ((void)0)
#else
#define SUSA_LOG_DBG(MSG) (std::cout << "[DBG]"  << "[" << __func__ << "()]" << "[" << __LINE__ << "]  : " << MSG << std::endl)
#define SUSA_LOG_INF(MSG) (std::cout << "[INF]"  << "[" << __func__ << "()]" << "[" << __LINE__ << "]  : " << MSG << std::endl)
#define SUSA_LOG_ERR(MSG) (std::cerr << "[ERR]"  << "[" << __func__ << "()]" << "[" << __LINE__ << "]  : " << MSG << std::endl)
#endif

#endif // SUSA_DEBUG_H
