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
#include "simplex.h"

static int verbose = 1;

void simplex(int n, double *xx, double scale, double (*fn)(int n, double *, void *), void *args, int max_iter)
{
    double **x;
    double *f;
    double *xr;
    double *xe;
    double *xc;
    double *x0;
    int n1, i, j;
    double value;
    double fr, fe, fc;
    double alpha = 1.0;
    double gamma = 2.0;
    double rho = -0.5;
    double sigma = 0.5;
    int iter;

    n1 = n + 1;
    x = malloc(n1 * sizeof(double *));
    for (i=0; i<n1; i++) {
        x[i] = malloc(n * sizeof(double));
    }

    for (i=0; i<n; i++) {
        memcpy(x[i], xx, n*sizeof(double));
        x[i][i] += scale;
    }
    memcpy(x[n], xx, n*sizeof(double));

    f = malloc(n1 * sizeof(double));
    xr = malloc(n * sizeof(double));
    xe = malloc(n * sizeof(double));
    xc = malloc(n * sizeof(double));
    x0 = malloc(n * sizeof(double));

    for (i=0; i<n1; i++) {
        value = fn(n, x[i], args);
        f[i] = value;
    }

    iter = 0;

    while (iter < max_iter) {

        /* Now, crap bubble sort.  We'll only have small problems to deal with. */
        for (i=0; i<n; i++) {
            for (j=i+1; j<n1; j++) {
                if (f[i] > f[j]) {
                    double t;
                    double *tt;
                    t = f[i];
                    f[i] = f[j];
                    f[j] = t;
                    tt = x[i];
                    x[i] = x[j];
                    x[j] = tt;
                }
            }
        }

        if (f[0] == f[n]) {
            /* obviously it has converged */
            break;
        }

        iter++;
        if (verbose) {
            printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
            printf("Iteration %d\n", iter);
            for (j=0; j<n1; j++) {
                printf("%2d %20.12f :", j, f[j]);
                for (i=0; i<n; i++) {
                    printf(" %20.12f", x[j][i]);
                }
                printf("\n");
            }
        }

        /* COG */
        for (i=0; i<n; i++) { /* coordinate */
            x0[i] = 0.0;
            for (j=0; j<n1; j++) { /* contributing point */
                x0[i] += x[j][i];
            }
            x0[i] /= (double) n1;
        }
        if (verbose) {
            printf("COG %20.12f :", f[0]);
            for (i=0; i<n; i++) {
                printf(" %20.12f", x0[i]);
            }
            printf("\n");
        }

        /* Reflect */
        for (i=0; i<n; i++) {
            xr[i] = x0[i] + alpha * (x0[i] - x[n][i]);
        }
        fr = fn(n, xr, args);
        if ((fr < f[n-1]) && (fr >= f[0])) {
            memcpy(x[n], xr, n*sizeof(double));
            f[n] = fr;
            if (verbose) printf("Reflected\n");
            continue;
        }

        if (fr < f[0]) {
            /* Otherwise, expand */
            for (i=0; i<n; i++) {
                xe[i] = x0[i] + gamma * (x0[i] - x[n][i]);
            }
            fe = fn(n, xe, args);
            if (fe < fr) {
                memcpy(x[n], xe, n*sizeof(double));
                f[n] = fe;
                if (verbose) printf("Expanded\n");
            } else {
                memcpy(x[n], xr, n*sizeof(double));
                f[n] = fr;
                if (verbose) printf("Reflected after expand\n");
            }
            continue;
        } else {
            /* Otherwise, contract */
            for (i=0; i<n; i++) {
                xc[i] = x0[i] + rho * (x0[i] - x[n][i]);
            }
            fc = fn(n, xc, args);
            if (fc < f[n]) {
                memcpy(x[n], xc, n*sizeof(double));
                f[n] = fc;
                if (verbose) printf("Contracted\n");
                continue;
            } else {
                /* Reduction */
                for (i=1; i<n1; i++) {
                    for (j=0; j<n; j++) {
                        x[i][j] = x[0][j] + sigma * (x[i][j] - x[0][j]);
                    }
                    f[i] = fn(n, x[i], args);
                }
                if (verbose) printf("Reduced\n");
            }
        }
    }

    for (i=0; i<n; i++) {
        /* Keep the best point */
        xx[i] = x[0][i];
    }

    for (i=0; i<n; i++) {
        free(x[i]);
    }
    free(x);
    free(f);
    free(xr);
    free(xe);
    free(xc);
    free(x0);
    return;
}

#if TEST

static double quadra(int n, double *x, void *args)
{
    double a, b, c;
    double ab, ac, bc1;
    a = x[0];
    b = x[1];
    c = x[2];

    ac = a - c + 2;
    ab = a + b;
    bc1 = b - c - 4;

    return sqrt(a*a+2.0*b*b+1.1*ac*ac) + ab*ab + bc1*bc1 + 3.0*ac*ac + 1.5*a*b*c;

}

int main (int argc, char **argv)
{
    double x[3];
    x[0] = 5.0;
    x[1] = 5.0;
    x[2] = 5.0;
    simplex(3, x, 0.01, quadra, NULL);
}

#endif

/* vim:et:sw=4:sts=4:ht=4
 * */

