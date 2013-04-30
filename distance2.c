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

/* Find approximation to distance between two points.  Actually, the scaling
 * per unit web-mercator space, as a function of the Y-coordinate
 * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tool.h"

static const double a = 6378137.0;
static const double b = 6356752.3141;

static double latofy(double y) {
    double t = (1.0 - 2.0*y) * M_PI;
    double M = (2.0 * atan(exp(t))) - 0.5*M_PI;
    return M;
}

static double latdegofy(double y) {
    double M = latofy(y);
    return M * (180.0/M_PI);
}

static double fy(double y)
{
    /* given y in [0,1] as the Y-coord of the web mercator position,
     * compute the 'radius', equivalent to Re * cos(lat) */

    double M = latofy(y);
    double e2 = 1.0 - (b*b)/(a*a);
    double SM = sin(M);
    double nu = a / sqrt(1.0 - e2*SM*SM);
    return nu * cos(M);
}

static void do_fit(int n, double *coef)
{
    Matrix m;
    long double c[32], r[32];
    double nodes[32];
    int i, j;

    for (i=1; i<=n; i++) {
        /* In range [-1,1] */
        double node = cos(M_PI*(double)(i+i-1)/(double)(n+n));
        printf("%2d : %12.5f\n", i, node);

#if 1
        /* Only cover the range 71N-71S */
        node = (0.5-0.215)*node + 0.5;
#endif

#if 0
        /* Only cover the range 75N-75S */
        node = (0.5-0.18)*node + 0.5;
#endif

#if 0
        /* Scale/shift to range [0,1] */
        node = 0.5*node + 0.5;
#endif
        printf("%2d : %10.5f\n", i, node);
        nodes[i-1] = node;
    }

    for (i=0; i<n; i++) {
        double nn = nodes[i];
        double prod = 1.0;
        double z = nn - 0.5;
        for (j=0; j<n; j++) {
            m[i][j] = prod;
            prod *= z;
        }
        r[i] = a/fy(nn);
    }

    solve(m, n, c, r);
    for (i=0; i<n; i++) {
        coef[i] = c[i];
    }
    return;

}

static double dmxy_to_metres(int n, double *coef, int ny, double *ycoef, struct mxy *p1, struct mxy *p2)
{
    double dx = p1->X - p2->X;
    double dy = p1->Y - p2->Y;
    double ay = 0.5 * (p1->Y + p2->Y);
    double az = ay - 0.5;
    double az2 = az * az;
    double total, prod;
    double scale_factor;
    double y_scale_factor;
    double d;
    int i;
    total = 1.0;
    prod = az2;
    for (i=2; i<n; i+=2) {
        total += prod * coef[i];
        prod *= az2;
    }
    scale_factor = (2.0 * M_PI * a) / total;
    if (ny > 0) {
        total = 0.0;
        prod = 1.0;
        for (i=0; i<ny; i+=2) {
            total += prod * ycoef[i];
            prod *= az2;
        }
        y_scale_factor = 1.0 / total;
    } else {
        y_scale_factor = 1.0;
    }
    dy /= y_scale_factor;
    d = sqrt(dx*dx + dy*dy);
    return d * scale_factor;
}

double dmxy_to_metres_C(struct mxy *p1, struct mxy *p2)
{
    const double KK = (2.0*M_PI*a);
    const double K0 = 1.0/KK;
    const double K2 = 19.42975297/KK;
    const double K4 = 74.22319781/KK;
    const double J0 = 0.99330562;
    const double J2 = 0.18663111;
    const double J4 = -1.45510549;
    double dx = p1->X - p2->X;
    double dy = p1->Y - p2->Y;
    double ay = 0.5 * (p1->Y + p2->Y);
    double az = ay - 0.5;
    double az2 = az * az;
    double az4 = az2 * az2;
    double sfy = J0 + J2*az2 + J4*az4;
    double sft = K0 + K2*az2 + K4*az4;
    double Dy = dy * sfy;
    double d = sqrt(dx*dx + Dy*Dy);
    return d / sft;
}

static void run_test(int n, double *coef, int ny, double *ycoef)
{
    struct mxy base, o1, o2;
    base.X = 0.45;
    base.Y = 0.33;

    struct llh base_llh, o1_llh, o2_llh;
    mxy_to_wgs84(&base, &base_llh);


    double actual;
    double predicted;
    int octa;
    int iteration;
    int octa_n = 6;
    double scaling = pow(10.0, 1.0/4.0);
    double step = 0.000001;
    for (iteration=0; iteration<4*5; iteration++) {
        for (octa=0; octa<octa_n; octa++) {
            double ang = M_PI * (double) octa / (double) octa_n;

            o1.X = base.X + step*sin(ang);
            o1.Y = base.Y + step*cos(ang);
            o2.X = base.X - step*sin(ang);
            o2.Y = base.Y - step*cos(ang);
            mxy_to_wgs84(&o1, &o1_llh);
            mxy_to_wgs84(&o2, &o2_llh);
            actual = llh_deg_to_metres(&o1_llh, &o2_llh);
#if 0
            predicted = dmxy_to_metres(n, coef, ny, ycoef, &o1, &o2);
#endif
            predicted = dmxy_to_metres_C(&o1, &o2);
            printf("%15.8f %6.1f %12.3f %12.3f %12.3f\n",
                step, (180.0/M_PI)*ang,
                actual, predicted, 1000.0*fabs(actual-predicted)/actual);
        }
        step *= scaling;
    }

    return;
}

static void run_data(const char *filename, double offset)
{
    struct mxy base, o1, o2, tmp;
    FILE *out;
    base.X = 0.50;

    struct llh o1_llh, o2_llh, tmp_llh;

    double actual;
    double predicted;
    int octa;
    int octa_n = 36;
    double y_incr = 0.005;
    double bound0 = 0.19 + 0.5*y_incr;
    double bound1 = 1.0 - bound0;
    double y;
    y = bound0;

    out = fopen(filename, "w");
    /* Avoid using 0.5 (the Equator) : we get NaN from somewhere then. */
    while (y < bound1) {
        for (octa=0; octa<octa_n; octa++) {
            double ang = M_PI * (double) octa / (double) octa_n;
            double ymod = y + (double)octa * 0.5 * y_incr / (double)octa_n;
            tmp.X = 0.5;
            tmp.Y = ymod;
            mxy_to_wgs84(&tmp, &tmp_llh);

            o1.X = base.X + offset*sin(ang);
            o1.Y = y + offset*cos(ang);
            o2.X = base.X - offset*sin(ang);
            o2.Y = y - offset*cos(ang);
            mxy_to_wgs84(&o1, &o1_llh);
            mxy_to_wgs84(&o2, &o2_llh);
            actual = llh_deg_to_metres(&o1_llh, &o2_llh);
            predicted = dmxy_to_metres_C(&o1, &o2);
            fprintf(out, "%12.5f %12.5f %12.3f %12.3f %12.3f\n",
                ymod, tmp_llh.lat,
                actual, predicted, 1000.0*fabs(actual-predicted)/actual);
        }
        y += y_incr;
    }
    fclose(out);

    return;
}

static double calc_y_scaling(int n, double *coef, int ny, double *ycoef, double y, double *Actual, double *Predicted)
{
    struct mxy base, o1, o2;
    base.X = 0.5;
    double step = 0.000001;
    struct llh o1_llh, o2_llh;
    double actual;
    double predicted;

    base.Y = y;
    o1.X = base.X;
    o1.Y = base.Y + step;
    o2.X = base.X;
    o2.Y = base.Y - step;
    mxy_to_wgs84(&o1, &o1_llh);
    mxy_to_wgs84(&o2, &o2_llh);
    actual = llh_deg_to_metres(&o1_llh, &o2_llh);
    predicted = dmxy_to_metres(n, coef, ny, ycoef, &o1, &o2);
    if (Actual) *Actual = actual;
    if (Predicted) *Predicted = predicted;

    return predicted / actual;
}

static void y_scaling(int n, double *coef, int ny, double *ycoef)
{
    double y;
    double actual;
    double predicted;
    y = 0.20;
    while (y < 0.80) {
        double total, prod;
        double z = y - 0.5;
        double z2 = z*z;
        int i;
        total = 0.0;
        prod = 1.0;
        for (i=0; i<ny; i+=2) {
            total += prod * ycoef[i];
            prod *= z2;
        }
        double ratio = calc_y_scaling(n, coef, 0, NULL, y, &actual, &predicted);
#if 1
        double predicted = 1.0 / total;
#else
        double predicted = total;
#endif
        printf("%8.3f %15.8f %15.8f %15.8f\n",
                y, ratio, predicted, predicted/ratio);

        y += 0.01;
    }
    return;
}

static void fit_y_scaling(int n, double *coef, int ny, double *ycoef)
{
    Matrix m;
    long double c[32], r[32];
    double nodes[32];
    int i, j;

    for (i=1; i<=ny; i++) {
        /* In range [-1,1] */
        double node = cos(M_PI*(double)(i+i-1)/(double)(ny+ny));

#if 1
        /* Only cover the range 71N-71S */
        node = (0.5-0.215)*node + 0.5;
#endif
        nodes[i-1] = node;
    }

    for (i=0; i<ny; i++) {
        double nn = nodes[i];
        double prod = 1.0;
        double z = nn - 0.5;
        for (j=0; j<ny; j++) {
            m[i][j] = prod;
            prod *= z;
        }
#if 1
        r[i] = 1.0 / calc_y_scaling(n, coef, 0, NULL, nn, NULL, NULL);
#else
        r[i] = calc_y_scaling(n, coef, nn, 0, NULL, NULL, NULL);
#endif
    }

    solve(m, ny, c, r);
    for (i=0; i<ny; i++) {
        ycoef[i] = c[i];
    }
    return;

}

int main (int argc, char **argv)
{
#if 1
    double coef[32];
    double ycoef[32];
    int i;
    int n = 5;
    int ny = 5;

    do_fit(n, coef);
    for (i=0; i<n; i++) {
        printf("%2d : %15.8f\n", i, coef[i]);
    }

    fit_y_scaling(n, coef, ny, ycoef);
    for (i=0; i<ny; i++) {
        printf("%2d : %15.8f\n", i, ycoef[i]);
    }
#endif

#if 0
    double y;
    double z, z2, fz;
    double fyy;
    double err;
    y = 0.0;
    while (y < 1.0) {
        z = (y - 0.5);
        z2 = z*z;
        double total, prod;
        total = 1.0;
        prod = z2;
        for (i=2; i<n; i+=2) {
            total += prod * coef[i];
            prod *= z2;
        }
        fz = a / total;
        fyy = fy(y);
        err = (fyy - fz)/fyy;

        printf("%15.8f %10.2f %15.8f %15.8f %15.8f %8.2f\n",
                y, latdegofy(y), fyy, fz, err, err*1000);
        y += 0.01;
    }
#endif

#if 0
    y_scaling(n, coef, ny, ycoef);
#endif

#if 0
    run_test(n, coef, ny, ycoef);
#endif

    run_data("geodesic_200m.dat", 0.0000025); /* +/- 100m around centre point at Equator */
    run_data("geodesic_2km.dat", 0.000025); /* +/- 1km around centre point at Equator */
    run_data("geodesic_20km.dat", 0.00025); /* +/- 10km around centre point at Equator */
    run_data("geodesic_200km.dat", 0.0025); /* +/- 100km around centre point at Equator */
    run_data("geodesic_2000km.dat", 0.025); /* +/- 1000km around centre point at Equator */

    return 0;
}

/* vim:et:sw=4:sts=4:ht=4
 * */

