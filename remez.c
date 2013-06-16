/*
* Copyright (c) 2013, Richard P. Curnow
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of the <organization> nor the
*       names of its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
* */

/* Remez exchange algorithm */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tool.h"
#include "remez.h"

static long double poly(int n, long double *coef, long double x)
{
    long double prod = 1.0;
    int i;
    long double result;
    result = 0.0;
    for (i=0; i<n; i++) {
        result += coef[i] * prod;
        prod *= x;
    }
    return result;
}
static double poly2(int n, double *coef, double x)
{
    double prod = 1.0;
    int i;
    double result;
    result = 0.0;
    for (i=0; i<n; i++) {
        result += coef[i] * prod;
        prod *= x;
    }
    return result;
}
void remez(double a, double b, double (*fn)(double, void *), void *args, int n, double *coef, int verbose)
{
    /* Find a polynomial up to degree (n-1) with coefficients to be stored in coef,
     * that is a minimax fit to (*fn) over the range [a,b].  Assume some upper
     * bound on n */

    long double x[256]; /* The points where the error is maximal */
    long double root[256];
    long double xnew[256];
    int i, n1;
    double ab = 0.5*(a+b);
    double spread = 0.5*(b-a);
    double e;
    Matrix m;
    long double c[32], r[32];
    int iter;

    /* Degree (n-1).
     * #coefs = n
     * #roots = n
     * #extrema = (n+1) = n1
     * */
    n1 = n + 1;

    /* Compute initial choice of x's : the cheby nodes */
    for (i=0; i<n1; i++) {
        int i1 = i + 1;
        /* negate to get the nodes into ascending order */
        double node = -cos(M_PI * (double)(i1+i1-1) / (double)(n1+n1));
        x[i] = ab + node*spread;
        if (verbose) {
            printf("x[%d] = %20.12Lf\n", i, x[i]);
        }
    }

    for (iter=1; iter<10; iter++) {
        if (verbose) {
            printf("=================\n");
            printf("Iteration %d\n", iter);
        }
        /* Solve for the coefs */
        for (i=0; i<n1; i++) {
            long double prod = 1.0;
            long double z = x[i];
            int j;
            for (j=0; j<n; j++) {
                m[i][j] = prod;
                prod *= z;
            }
            m[i][n] = (i & 1) ? 1.0 : -1.0;
            r[i] = fn(z, args);
        }

        solve(m, n1, c, r);
        for (i=0; i<n; i++) {
            coef[i] = c[i];
        }
        e = c[n];
        if (verbose) {
            printf("e=%20.12e\n", e);
        }

        /* Find roots */
        for (i=0; i<n; i++) {
            long double lo, hi, mid, flo, fhi, fmid;
            long double epsilon = 1.0e-16;
            lo = x[i];
            hi = x[i+1];
            flo = fn(lo, args) - poly(n, c, lo);
            fhi = fn(hi, args) - poly(n, c, hi);
            /* bisection search */
            while ((hi - lo) > epsilon) {
                mid = 0.5 * (hi + lo);
                fmid = fn(mid, args) - poly(n, c, mid);
                if ((flo >= 0.0) ^ (fmid >= 0.0)) {
                    /* low and mid have opposite signs */
                    hi = mid;
                    fhi = fmid;
                } else {
                    lo = mid;
                    flo = fmid;
                }
            }
            root[i] = lo;
            if (verbose) {
                printf("root[%d] = %25.18Lf\n", i, lo);
            }
        }

        /* Find extrema.  Surely the endpoints always have to be? */
        xnew[0] = a;
        for (i=1; i<n; i++) {
            long double lo, hi;
            long double epsilon = 1.0e-12;
            long double gold = 0.618;
            long double mid0, mid1;
            lo = root[i-1];
            hi = root[i];
            while (hi - lo > epsilon) {
                long double f1, f0;
                mid1 = lo + gold * (hi - lo);
                mid0 = hi - gold * (hi - lo);
                f1 = fn(mid1, args) - poly(n, c, mid1);
                f0 = fn(mid0, args) - poly(n, c, mid0);
                if (fabs(f1) > fabs(f0)) {
                    /* extremum is in the upper portion */
                    lo = mid0;
                } else {
                    hi = mid1;
                }
            }
            xnew[i] = lo;
            if (verbose) {
                printf("extreme[%d] = %25.18Lf\n", i, lo);
            }

        }
        xnew[n] = b;

        for (i=0; i<n1; i++) {
            x[i] = xnew[i];
        }
    }

    return;
}

#if TEST

static double my_atan(double x, void *unused) {
    return atan2(x, 1.0);
}

static void plot(int n, double *coef, char *name)
{
    FILE *out = fopen(name, "w");
    double x, value, error;
    for(x=-1.01; x<=1.01; x+=0.01) {
        value = poly2(n, coef, x);
        error = value - my_atan(x, NULL);
        fprintf(out, "%20.12f %20.12f %20.12f %20.12f\n", x, value, error, error*180.0/M_PI);
    }
    fclose(out);
}

#define K1 (0.97179803008)
#define K2 (0.19065470515)

int main (int argc, char **argv)
{
    double coef[64];
    int n;
    int i;


    n = 5;
    remez(-1.0, 1.0, my_atan, NULL, n, coef, 1);
    for (i=0; i<n; i++) {
        printf("%2d : %20.12f\n", i, coef[i]);
    }
    plot(n, coef, "remez_cubic.dat");

    n = 7;
    remez(-1.0, 1.0, my_atan, NULL, n, coef, 1);
    for (i=0; i<n; i++) {
        printf("%2d : %20.12f\n", i, coef[i]);
    }
    plot(n, coef, "remez_quintic.dat");

    n = 9;
    remez(-1.0, 1.0, my_atan, NULL, n, coef, 1);
    for (i=0; i<n; i++) {
        printf("%2d : %20.12f\n", i, coef[i]);
    }
    plot(n, coef, "remez_deg7.dat");

    n = 13;
    remez(-1.0, 1.0, my_atan, NULL, n, coef, 1);
    for (i=0; i<n; i++) {
        printf("%2d : %20.12f\n", i, coef[i]);
    }
    plot(n, coef, "remez_deg11.dat");

    coef[0] = 0.0;
    coef[2] = 0.0;
    coef[4] = 0.0;
    coef[1] = K1;
    coef[3] = -K2;
    plot(n, coef, "web_model.dat");

    return 0;

}

#endif


/* vim:et:sw=4:sts=4:ht=4
 * */
