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
 * @file svd.cpp
 * @brief Singular Value Decomposition (SVD)
 * @author Behrooz Kamary Aliabadi
 * @version 1.0.0
 */

#include "susa.h"

namespace susa {

#define MAX(x,y) ((x)>(y)?(x):(y))
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static double _susa_svd_pythag( double a, double b ) {
    /* This routine is adapted from svdecomp.c in XLISP-STAT 2.1 which is
     * a code from Numerical Recipes adapted by Luke Tierney and David Betz.
     * http://www.public.iastate.edu/~dicook/JSS/paper/code/svd.c
     * http://svn.lirec.eu/libs/magicsquares/src/SVD.cpp
     */
    double at = fabs(a);
    double bt = fabs(b);
    double bat = bt / at;
    double abt = at / bt;
    if( at > bt )
        return at * sqrt( 1.0 + bat * bat );
    else if( bt > 0.0 )
        return bt * sqrt( 1.0 + abt * abt );
    else return 0.0;
}


int svd(matrix <float> &mat_arg_a, matrix <float> &mat_arg_w, matrix <float> &mat_arg_v) {
    /* This routine is adapted from svdecomp.c in XLISP-STAT 2.1 which is
     * a code from Numerical Recipes adapted by Luke Tierney and David Betz.
     * http://www.public.iastate.edu/~dicook/JSS/paper/code/svd.c
     * http://svn.lirec.eu/libs/magicsquares/src/SVD.cpp
    */
    int m = mat_arg_a.no_rows();
    int n = mat_arg_a.no_cols();

    int flag, i, its, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    matrix <double> rv1(n,1);


    if ( mat_arg_a.no_rows() < mat_arg_a.no_cols()) {
        // The rest of the matrix should be filled with zeros.
        return(0);
    }

    /* Householder reduction to bidiagonal form */
    for (i = 0; i < n; i++) {
        /* left-hand reduction */
        l = i + 1;
        rv1(i) = scale * g;
        g = s = scale = 0.0;
        if (i < m) {
            for (k = i; k < m; k++)
                scale += fabs((double)mat_arg_a(k,i));
            if (scale) {
                for (k = i; k < m; k++) {
                    mat_arg_a(k,i) = (float)((double)mat_arg_a(k,i)/scale);
                    s += ((double)mat_arg_a(k,i) * (double)mat_arg_a(k,i));
                }
                f = (double)mat_arg_a(i,i);
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                mat_arg_a(i,i) = (float)(f - g);
                if (i != n - 1) {
                    for (j = l; j < n; j++) {
                        for (s = 0.0, k = i; k < m; k++)
                            s += ((double)mat_arg_a(k,i) * (double)mat_arg_a(k,j));
                        f = s / h;
                        for (k = i; k < m; k++)
                            mat_arg_a(k,j) += (float)(f * (double)mat_arg_a(k,i));
                    }
                }
                for (k = i; k < m; k++)
                    mat_arg_a(k,i) = (float)((double)mat_arg_a(k,i)*scale);
            }
        }
        mat_arg_w(i) = (float)(scale * g);

        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != n - 1) {
            for (k = l; k < n; k++)
                scale += fabs((double)mat_arg_a(i,k));
            if (scale) {
                for (k = l; k < n; k++) {
                    mat_arg_a(i,k) = (float)((double)mat_arg_a(i,k)/scale);
                    s += ((double)mat_arg_a(i,k) * (double)mat_arg_a(i,k));
                }
                f = (double)mat_arg_a(i,l);
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                mat_arg_a(i,l) = (float)(f - g);
                for (k = l; k < n; k++)
                    rv1(k) = (double)mat_arg_a(i,k) / h;
                if (i != m - 1) {
                    for (j = l; j < m; j++) {
                        for (s = 0.0, k = l; k < n; k++)
                            s += ((double)mat_arg_a(j,k) * (double)mat_arg_a(i,k));
                        for (k = l; k < n; k++)
                            mat_arg_a(j,k) += (float)(s * rv1(k));
                    }
                }
                for (k = l; k < n; k++)
                    mat_arg_a(i,k) = (float)((double)mat_arg_a(i,k)*scale);
            }
        }
        anorm = MAX(anorm, (fabs((double)mat_arg_w(i)) + fabs(rv1(i))));
    }

    /* accumulate the right-hand transformation */
    for (i = n - 1; i >= 0; i--) {
        if (i < n - 1) {
            if (g) {
                for (j = l; j < n; j++)
                    mat_arg_v(j,i) = (float)(((double)mat_arg_a(i,j) / (double)mat_arg_a(i,l)) / g);
                /* double division to avoid underflow */
                for (j = l; j < n; j++) {
                    for (s = 0.0, k = l; k < n; k++)
                        s += ((double)mat_arg_a(i,k) * (double)mat_arg_v(k,j));
                    for (k = l; k < n; k++)
                        mat_arg_v(k,j) += (float)(s * (double)mat_arg_v(k,i));
                }
            }
            for (j = l; j < n; j++)
                mat_arg_v(i,j) = mat_arg_v(j,i) = 0.0;
        }
        mat_arg_v(i,i) = 1.0;
        g = rv1(i);
        l = i;
    }

    /* accumulate the left-hand transformation */
    for (i = n - 1; i >= 0; i--) {
        l = i + 1;
        g = (double)mat_arg_w(i);
        if (i < n - 1)
            for (j = l; j < n; j++)
                mat_arg_a(i,j) = 0.0;
        if (g) {
            g = 1.0 / g;
            if (i != n - 1) {
                for (j = l; j < n; j++) {
                    for (s = 0.0, k = l; k < m; k++)
                        s += ((double)mat_arg_a(k,i) * (double)mat_arg_a(k,j));
                    f = (s / (double)mat_arg_a(i,i)) * g;
                    for (k = i; k < m; k++)
                        mat_arg_a(k,j) += (float)(f * (double)mat_arg_a(k,i));
                }
            }
            for (j = i; j < m; j++)
                mat_arg_a(j,i) = (float)((double)mat_arg_a(j,i)*g);
        } else {
            for (j = i; j < m; j++)
                mat_arg_a(j,i) = 0.0;
        }
        ++mat_arg_a(i,i);
    }

    /* diagonalize the bidiagonal form */
    for (k = n - 1; k >= 0; k--) {                             /* loop over singular values */
        for (its = 0; its < 30; its++) {                       /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--) {                         /* test for splitting */
                nm = l - 1;
                if (fabs(rv1(l)) + anorm == anorm) {
                    flag = 0;
                    break;
                }
                if (fabs((double)mat_arg_w(nm)) + anorm == anorm)
                    break;
            }
            if (flag) {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++) {
                    f = s * rv1(i);
                    if (fabs(f) + anorm != anorm) {
                        g = (double)mat_arg_w(i);
                        h = _susa_svd_pythag(f, g);
                        mat_arg_w(i) = (float)h;
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++) {
                            y = (double)mat_arg_a(j,nm);
                            z = (double)mat_arg_a(j,i);
                            mat_arg_a(j,nm) = (float)(y * c + z * s);
                            mat_arg_a(j,i) = (float)(z * c - y * s);
                        }
                    }
                }
            }
            z = (double)mat_arg_w(k);
            if (l == k) {                  /* convergence */
                if (z < 0.0) {             /* make singular value nonnegative */
                    mat_arg_w(k) = (float)(-z);
                    for (j = 0; j < n; j++)
                        mat_arg_v(j,k) = - mat_arg_v(j,k);
                }
                break;
            }
            if (its >= 30) {
                //free((void*) rv1);
                fprintf(stderr, "No convergence after 30,000! iterations \n");
                return(0);
            }

            /* shift from bottom 2 x 2 minor */
            x = (double)mat_arg_w(l);
            nm = k - 1;
            y = (double)mat_arg_w(nm);
            g = rv1(nm);
            h = rv1(k);
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = _susa_svd_pythag(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++) {
                i = j + 1;
                g = rv1(i);
                y = (double)mat_arg_w(i);
                h = s * g;
                g = c * g;
                z = _susa_svd_pythag(f, h);
                rv1(j) = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < n; jj++) {
                    x = (double)mat_arg_v(jj,j);
                    z = (double)mat_arg_v(jj,i);
                    mat_arg_v(jj,j) = (float)(x * c + z * s);
                    mat_arg_v(jj,i) = (float)(z * c - x * s);
                }
                z = _susa_svd_pythag(f, h);
                mat_arg_w(j) = (float)z;
                if (z) {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++) {
                    y = (double)mat_arg_a(jj,j);
                    z = (double)mat_arg_a(jj,i);
                    mat_arg_a(jj,j) = (float)(y * c + z * s);
                    mat_arg_a(jj,i) = (float)(z * c - y * s);
                }
            }
            rv1(l) = 0.0;
            rv1(k) = f;
            mat_arg_w(k) = (float)x;
        }
    }
    return 0;
}

} // NAMESPACE SUSA
