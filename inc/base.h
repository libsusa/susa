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
 * @file base.h
 * @brief Basic mathematical operations on STL and Susa types
 * @author Behrooz, Kamary Aliabadi
 * @version 1.0.0
 *
 * @defgroup Math Basic Mathematics
 */

#ifndef BASE_H
#define BASE_H

namespace susa {
/**
 * @brief Min
 * 
 * @param mat_arg
 * @return Returns minimum value in the input vector
 * @ingroup Math
 */
template <class T> matrix <unsigned int> min(const matrix <T> &mat_arg);



/**
 * @brief Max
 * 
 * @param mat_arg
 * @return Returns maximum value in the input vector
 * @ingroup Math
 */
template <class T> matrix <unsigned int> max(const matrix <T> &mat_arg);



/**
 * @brief Differential
 * 
 * @param vec_arg STL vector
 * @return Returns discrete differential of the input vector
 * @ingroup Math
 */
template <class T> std::vector <T> diff(std::vector <T> &vec_arg);

/**
 * @brief Sum
 * 
 * @param mat_arg
 * @return Returns sum of values in the input vector
 * @ingroup Math
 */
template <class T> matrix <T> sum(matrix <T> &mat_arg);

/**
 * @brief Mean
 * 
 * @param mat_arg
 * @return Returns mean of values in the input vector
 * @ingroup Math
 */
template <class T> matrix <double> mean(matrix <T> &mat_arg);

/**
 * @brief Magnitude
 * 
 * @param mat_arg
 * @return Returns magnitude of the complex input vector
 * @ingroup Math
 */
template <class T> matrix <double> mag(const matrix < std::complex <T> > &mat_arg);

/**
 * @brief Norm
 * 
 * @param mat_arg
 * @return Returns magnitude of the complex input vector
 * @ingroup Math
 */
template <class T> matrix <double> norm(const matrix <T> &mat_arg);

/**
 * @brief Real
 * 
 * @param vec_arg STL vector
 * @return Returns real part of the complex input vector
 * @ingroup Math
 */
template <class T> matrix <T> real(const matrix < std::complex <T> > &mat_arg);

/**
 * @brief Imaginary
 * 
 * @param vec_arg STL vector
 * @return Returns imaginary part of the complex input vector
 * @ingroup Math
 */
template <class T> matrix <T> imag(const matrix < std::complex <T> > &mat_arg);

/**
 * @brief Absolute
 * 
 * @param vec_arg STL vector
 * @return Returns absolut value of the input vector
 * @ingroup Math
 */
template <class T> std::vector <T> abs(std::vector <T> &vec_arg);

/**
 * @brief Conjugate
 * 
 * @param vec_arg STL vector
 * @return Returns conjugate of the complex input vector
 * @ingroup Math
 */
template <class T> matrix <std::complex <T> > conj(const matrix < std::complex <T> > &mat_arg);

/**
 * @brief Sign
 * 
 * @param mat_arg input argument
 * @return Returns sign of input
 * @ingroup Math
 */
template <class T> matrix <T> sign(const matrix <T> &mat_arg);

template <class T> T sign(T T_arg);

/**
 * @brief Log
 * 
 * @param mat_arg input argument
 * @return Returns natural logarithm of input
 * @ingroup Math
 */
template <class T> matrix <T> log(const matrix <T> &mat_arg);

/**
 * @brief Pow
 * 
 * @param uint_b base argument
 * @param uint_p power argument
 * @ingroup Math
 */
unsigned int pow(unsigned int uint_b, unsigned int uint_p);

/**
 * @brief Pow
 * 
 * @param uint_b base argument
 * @param uint_p power argument
 * @ingroup Math
 */
int pow(int int_b, unsigned int uint_p);

/**
 * @brief ncdf
 * 
 * @param x input argument
 * @ingroup Math
 */
double ncdf(const double x);

/**
 * @brief qfunc
 * 
 * @param x input argument
 * @ingroup Math
 */
double qfunc(const double x);


/**
 * @brief Modular operation
 *
 * @return lint_a "mod" lint_mod 
 *
 * @ingroup Math
 */
long int mod(long int lint_a,long int lint_mod);




template <class T> std::vector <T> diff(std::vector <T> &vec_arg) {
  std::vector <T> vec_diff( vec_arg.size() - 1, 0 );

  for (int i=0;i<vec_arg.size()-1;i++) {
      vec_diff[i] = vec_arg[i+1] - vec_arg[i];
  }

  return (vec_diff);
}



template <class T> matrix <T> sum(matrix <T> &mat_arg) {
  matrix <T> mat_ret;

  if (mat_arg.no_rows() == 1 || mat_arg.no_cols() == 1) {
    mat_ret = matrix <T> (1,1);
    for (unsigned int uint_i = 0; uint_i < mat_arg.size(); uint_i++) mat_ret(0) += mat_arg(uint_i);
  } else {
    mat_ret = matrix <T> (1,mat_arg.no_cols());
    for (unsigned int uint_col = 0; uint_col < mat_arg.no_cols(); uint_col++) {
      for (unsigned int uint_row = 0; uint_row < mat_arg.no_rows(); uint_row++) mat_ret(uint_col) += mat_arg(uint_row);
    }
  }

  return mat_ret;
}

template <class T> matrix <double> mean(matrix <T> &mat_arg) {
  matrix <double> mat_ret;

  if (mat_arg.no_rows() == 1 || mat_arg.no_cols() == 1) {
    mat_ret = matrix <double> (1,1);
    for (unsigned int uint_i = 0; uint_i < mat_arg.size(); uint_i++) mat_ret(0) += mat_arg(uint_i);
  } else {
    mat_ret = matrix <double> (1,mat_arg.no_cols());
    for (unsigned int uint_col = 0; uint_col < mat_arg.no_cols(); uint_col++) {
      for (unsigned int uint_row = 0; uint_row < mat_arg.no_rows(); uint_row++) mat_ret(uint_col) += mat_arg(uint_row);
      mat_ret(uint_col) /= mat_arg.no_rows();
    }
  }

  return mat_ret;
}

template <class T> std::vector <T> abs(std::vector <T> &vec_arg) {
  std::vector <T> vec_abs(vec_arg.size(),0);
  for (int i=0;i<vec_arg.size();i++) vec_abs[i] = abs(vec_arg[i]);
  return vec_abs;
}


template <class T> matrix <unsigned int> min(const matrix <T> &mat_arg) {
  matrix <unsigned int> mat_ret;
  T T_min;

  if (mat_arg.no_cols() == 1 || mat_arg.no_rows() == 1) {
    mat_ret = matrix <unsigned int> (1,1);
    T_min = mat_arg(0);
    for (unsigned int uint_i = 1; uint_i < mat_arg.size(); uint_i++) {
      if (mat_arg(uint_i) < T_min) {
        T_min = mat_arg(uint_i);
	mat_ret(0) = uint_i;
      }
    }
  } else if (mat_arg.no_rows() > 1 && mat_arg.no_cols() > 1) {
    mat_ret = matrix <unsigned int> (1,mat_arg.no_cols());
    for (unsigned int uint_col = 0; uint_col < mat_arg.no_cols(); uint_col++) {
      T_min = mat_arg(0,uint_col);
      for (unsigned int uint_row = 1; uint_row < mat_arg.no_rows(); uint_row++) {
        if (mat_arg(uint_row, uint_col) < T_min) {
	  T_min = mat_arg(uint_row,uint_col);
	  mat_ret(uint_col) = uint_row;
	}
      }
    }
  }

  return mat_ret;
}


template <class T> matrix <unsigned int> max(const matrix <T> &mat_arg) {
  matrix <unsigned int> mat_ret;
  T T_max;

  if (mat_arg.no_cols() == 1 || mat_arg.no_rows() == 1) {
    mat_ret = matrix <unsigned int> (1,1);
    T_max = mat_arg(0);
    for (unsigned int uint_i = 1; uint_i < mat_arg.size(); uint_i++) {
      if (mat_arg(uint_i) > T_max) {
        T_max = mat_arg(uint_i);
	mat_ret(0) = uint_i;
      }
    }
  } else if (mat_arg.no_rows() > 1 && mat_arg.no_cols() > 1) {
    mat_ret = matrix <unsigned int> (1,mat_arg.no_cols());
    for (unsigned int uint_col = 0; uint_col < mat_arg.no_cols(); uint_col++) {
      T_max = mat_arg(0,uint_col);
      for (unsigned int uint_row = 1; uint_row < mat_arg.no_rows(); uint_row++) {
        if (mat_arg(uint_row, uint_col) > T_max) {
	  T_max = mat_arg(uint_row,uint_col);
	  mat_ret(uint_col) = uint_row;
	}
      }
    }
  }

  return mat_ret;
}


template <class T> matrix <std::complex <T> > conj(const matrix < std::complex <T> > &mat_arg) {
  matrix < std::complex <T> > mat_ret(mat_arg.no_rows(), mat_arg.no_cols());

  for (unsigned int uint_i = 0; uint_i < mat_arg.size(); uint_i++) mat_ret(uint_i) = std::conj(mat_arg(uint_i));
  return mat_ret;
}

template <class T> matrix <double> mag(const matrix < std::complex <T> > &mat_arg) {
  matrix <double> mat_ret(mat_arg.no_rows(), mat_arg.no_cols());

  for (unsigned int uint_i = 0; uint_i < mat_arg.size(); uint_i++) mat_ret(uint_i) = (mat_arg(uint_i) * std::conj(mat_arg(uint_i))).real();
  return mat_ret;
}

template <class T> matrix <double> norm(const matrix <T> &mat_arg) {
  matrix <double> mat_ret;

  if (mat_arg.no_cols() == 1 || mat_arg.no_rows() == 1) {
    mat_ret = matrix <double> (1,1);
    for (unsigned int uint_i = 0; uint_i < mat_arg.size(); uint_i++) mat_ret(0) += mat_arg(uint_i) * mat_arg(uint_i);
    mat_ret(0) = std::sqrt(mat_ret(0));
  } else if (mat_arg.no_rows() > 1 && mat_arg.no_cols() > 1) {
    mat_ret = matrix <double> (1,mat_arg.no_cols());
    for (unsigned int uint_col = 0; uint_col < mat_arg.no_cols(); uint_col++) {
      for (unsigned int uint_row = 0; uint_row < mat_arg.no_rows(); uint_row++) {
        mat_ret(uint_col) += mat_arg(uint_row,uint_col) * mat_arg(uint_row,uint_col);
      }
      mat_ret(uint_col) = std::sqrt(mat_ret(uint_col));
    }
  }
  return mat_ret;
}

template <class T> matrix <T> real(const matrix < std::complex <T> > &mat_arg) {
  matrix <T> mat_ret(mat_arg.no_rows(), mat_arg.no_cols());

  for (unsigned int uint_i = 0; uint_i < mat_arg.size(); uint_i++) mat_ret(uint_i) = mat_arg(uint_i).real();
  
  return mat_ret;
}



template <class T> matrix <T> imag(const matrix < std::complex <T> > &mat_arg) {
  matrix < std::complex <T> > mat_ret(mat_arg.no_rows(), mat_arg.no_cols());

  for (unsigned int uint_i = 0; uint_i < mat_arg.size(); uint_i++) mat_ret(uint_i) = mat_arg(uint_i).imag();
  
  return mat_ret;
}

template <class T> matrix <T> sign(const matrix <T> &mat_arg) {
  matrix <T> mat_ret(mat_arg.no_rows(), mat_arg.no_cols());

  for (unsigned int uint_i = 0; uint_i < mat_arg.size(); uint_i++) {
    if (mat_arg(uint_i) != 0) mat_ret(uint_i) = mat_arg(uint_i) > 0 ? 1 : -1;
  }
	
  return mat_ret;
}

template <class T> T sign(T T_arg) {
  T T_ret = 0;
  if (T_arg > 0) T_arg = 1;
  else if (T_arg < 0) T_arg = -1;

  return T_arg;
}

template <class T> matrix <T> log(const matrix <T> &mat_arg) {
  matrix <T> mat_ret(mat_arg.no_rows(), mat_arg.no_cols());

  for (unsigned int uint_i = 0; uint_i < mat_arg.size(); uint_i++) {
    mat_ret(uint_i) = std::log(mat_arg(uint_i));
  }
	
  return mat_ret;
}


} // NAMESPACE SUSA
#endif


