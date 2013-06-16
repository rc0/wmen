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
#include "remez.h"

struct params/*{{{*/
{
    int n;
    int ny;
    double coef[32];
    double ycoef[32];
};
/*}}}*/
static const double a = 6378137.0;
static const double b = 6356752.3141;

static double latofy(double y) {/*{{{*/
    double t = (1.0 - 2.0*y) * M_PI;
    double M = (2.0 * atan(exp(t))) - 0.5*M_PI;
    return M;
}
/*}}}*/
static double latdegofy(double y) {/*{{{*/
    double M = latofy(y);
    return M * (180.0/M_PI);
}
/*}}}*/
static double poly2(int n, double *coef, double x)/*{{{*/
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
/*}}}*/
static double fy(double y)/*{{{*/
{
    /* given y in [0,1] as the Y-coord of the web mercator position,
     * compute the 'radius', equivalent to Re * cos(lat) */

    double M = latofy(y);
    double e2 = 1.0 - (b*b)/(a*a);
    double SM = sin(M);
    double nu = a / sqrt(1.0 - e2*SM*SM);
    return nu * cos(M);
}
/*}}}*/
static void do_fit(int n, double *coef)/*{{{*/
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
        node = (0.5-0.211)*node + 0.5;
#endif
#if 0
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
/*}}}*/
static double dmxy_to_metres(int n, double *coef, int ny, double *ycoef, struct mxy *p1, struct mxy *p2)/*{{{*/
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
/*}}}*/
double dmxy_to_metres_C(struct mxy *p1, struct mxy *p2)/*{{{*/
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
/*}}}*/
double dmxy_to_metres_model(struct params *p, struct mxy *p1, struct mxy *p2)/*{{{*/
{
    int i;
    double prod;
    double dx = p1->X - p2->X;
    double dy = p1->Y - p2->Y;
    double ay = 0.5 * (p1->Y + p2->Y);
    double az = ay - 0.5;
    prod = 1.0;
    double sfy = 0.0;
    for (i=0; i<p->ny; i++) {
        sfy += p->ycoef[i] * prod;
        prod *= az;
    }
    prod = 1.0;
    double sft = 0.0;
    for (i=0; i<p->n; i++) {
        sft += p->coef[i] * prod;
        prod *= az;
    }
    double Dy = dy * sfy;
    double d = sqrt(dx*dx + Dy*Dy);
    return d * 2.0 * M_PI * a / sft;
}
/*}}}*/
static void run_test(int n, double *coef, int ny, double *ycoef)/*{{{*/
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
/*}}}*/
static void run_data(const char *filename, double offset)/*{{{*/
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
            fprintf(out, "%12.5f %12.5f %12.5f %12.5f %12.3f\n",
                ymod, tmp_llh.lat,
                actual, predicted, 1000.0*fabs(actual-predicted)/actual);
        }
        y += y_incr;
    }
    fclose(out);

    return;
}
/*}}}*/
static void run_data_with_model(const char *filename, struct params *p, double offset)/*{{{*/
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
    double bound0 = 0.20 + 0.5*y_incr;
    double bound1 = 1.0 - bound0;
    double y;
    y = bound0;

    out = fopen(filename, "w");
    /* Avoid using 0.5 (the Equator) : we get NaN from somewhere then. */
    while (y < bound1) {
        double min_err = 99999999.0;
        double max_err = -99999999.0;
        tmp.X = 0.5;
        tmp.Y = y;
        mxy_to_wgs84(&tmp, &tmp_llh);
        for (octa=0; octa<octa_n; octa++) {
            double ang = M_PI * (double) octa / (double) octa_n;
            double ymod = y + (double)octa * 0.5 * y_incr / (double)octa_n;
            double err;

            o1.X = base.X + offset*sin(ang);
            o1.Y = y + offset*cos(ang);
            o2.X = base.X - offset*sin(ang);
            o2.Y = y - offset*cos(ang);
            mxy_to_wgs84(&o1, &o1_llh);
            mxy_to_wgs84(&o2, &o2_llh);
            actual = llh_deg_to_metres(&o1_llh, &o2_llh);
            predicted = dmxy_to_metres_model(p, &o1, &o2);
            err = 1000.0 * (actual - predicted) / actual;
            if (err < min_err) min_err = err;
            if (err > max_err) max_err = err;
#if 0
            fprintf(out, "%12.5f %12.5f %12.3f %12.3f %12.3f\n",
                ymod, tmp_llh.lat,
                actual, predicted, 1000.0*fabs(actual-predicted)/actual);
#endif
        }
        fprintf(out, "%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n",
                y, tmp_llh.lat,
                0.5 * (max_err + min_err),
                0.5 * (max_err - min_err),
                min_err, max_err);
        y += y_incr;
    }
    fclose(out);

    return;
}
/*}}}*/
static double calc_y_scaling_old(int n, double *coef, int ny, double *ycoef, double y, double *Actual, double *Predicted)/*{{{*/
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
/*}}}*/
static double calc_y_scaling(int n, double *coef, int ny, double *ycoef, double y, double *Actual, double *Predicted)/*{{{*/
{
    struct mxy base, o1, o2;
    double offset = 0.000001;
    struct llh o1_llh, o2_llh;
    double actual;
    double predicted;
    double min_scaling = 9999999999.0;
    double max_scaling = -9999999999.0;
    int octa;
    int octa_n = 36;

    base.X = 0.5;
    base.Y = y;
    for (octa=0; octa<octa_n; octa++) {
        double ang = M_PI * (0.5 + (double) octa) / (double) octa_n;
        double scaling;
        o1.X = base.X + offset*sin(ang);
        o1.Y = y + offset*cos(ang);
        o2.X = base.X - offset*sin(ang);
        o2.Y = y - offset*cos(ang);
        mxy_to_wgs84(&o1, &o1_llh);
        mxy_to_wgs84(&o2, &o2_llh);
        actual = llh_deg_to_metres(&o1_llh, &o2_llh);
        predicted = dmxy_to_metres(n, coef, ny, ycoef, &o1, &o2);

        scaling = predicted / actual;
        if (scaling < min_scaling) min_scaling = scaling;
        if (scaling > max_scaling) max_scaling = scaling;

    }
    /* Mid point of scaling range */
    return 0.5 * (min_scaling + max_scaling);
}
/*}}}*/
static int compare_double(const void *a, const void *b)
{
    a = (double *) a;
    b = (double *) b;
    if (a < b) return -1;
    else if (a > b) return +1;
    else return 0;
}

static double calc_y_scaling_remez(int n, double *coef, double y, double *lo, double *hi)/*{{{*/
{
    struct mxy base, o1, o2;
    double offset = 0.000001;
    struct llh o1_llh, o2_llh;
    double actual;
    double predicted;
    double min_scaling = 9999999999.0;
    double max_scaling = -9999999999.0;
    int step;
    int step_n = 12;
    double product;
    int count;
    double result;
    double storage[1024];
    double total;
    double total_weight;
    double weight;

    base.X = 0.5;
    base.Y = y;
    product = 1.0;
    count = 0;
    total = 0.0;
    total_weight = 0.0;
    for (step=0; step<step_n; step++) {
        /* Only need to cover a quadrant */
        double ang = (0.5 * M_PI) * (0.5 + (double) step) / (double) step_n;
        double scaling;
        double q;
        double qq, qq2, qq3, qq4;
        double dx, dy;
        double sft;
        dx = offset * sin(ang);
        dy = offset * cos(ang);
        o1.X = base.X + dx;
        o1.Y = y + dy;
        o2.X = base.X - dx;
        o2.Y = y - dy;
        mxy_to_wgs84(&o1, &o1_llh);
        mxy_to_wgs84(&o2, &o2_llh);
        actual = llh_deg_to_metres(&o1_llh, &o2_llh);
        sft = poly2(n, coef, y-0.5);
        q = 2.0*M_PI*a / sft;
        qq = actual / q;
        qq2 = qq*qq;
        qq3 = qq2 - 4.0*dx*dx;
        qq4 = qq3 / (4.0*dy*dy);
        scaling = sqrt(qq4);
#if 0
        printf("y=%20.8f step=%2d sft=%20.10f actual=%20.10f q=%20.10f qq=%20.10f qq2=%20.10f qq3=%20.10f qq4=%20.10f scaling=%20.10f\n",
                y, step, sft, actual, q, qq, qq2, qq3, qq4, scaling);
#endif
        if (!isnan(scaling)) {
            if (scaling < min_scaling) min_scaling = scaling;
            if (scaling > max_scaling) max_scaling = scaling;
            product *= scaling;
            storage[count] = scaling;
            count++;

            if (fabs(dx) > 0.0) {
                weight = fabs(dy/dx);
                weight *= weight;
                total_weight += weight;
                total += scaling * weight;
            }
        }
    }
    /* Mid point of scaling range */
    if (lo) *lo = min_scaling;
    if (hi) *hi = max_scaling;
    return total / total_weight;
#if 0
    qsort(storage, count, sizeof(double), compare_double);
    return storage[count>>1]; /* median */
#endif
#if 0
    result = pow(product, 1.0 / (double) count);
    return result;
#endif
}
/*}}}*/
static void y_scaling(int n, double *coef, int ny, double *ycoef)/*{{{*/
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
/*}}}*/
static void fit_y_scaling(int n, double *coef, int ny, double *ycoef)/*{{{*/
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
        r[i] = 1.0 / calc_y_scaling_old(n, coef, 0, NULL, nn, NULL, NULL);
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
/*}}}*/
static double fn_for_remez(double x, void *args) {/*{{{*/
    return a/fy(x+0.5);
}
/*}}}*/
static double fnyscale_for_remez(double x, void *args) {/*{{{*/
    struct params *p = args;
    double nn = x + 0.5;
    return calc_y_scaling_remez(p->n, p->coef, nn, NULL, NULL);
}
/*}}}*/
static void fny2scale_for_remez(double x, void *args, double *lo, double *hi) {/*{{{*/
    struct params *p = args;
    double nn = x + 0.5;
    calc_y_scaling_remez(p->n, p->coef, nn, lo, hi);
    *lo = 1.0 / *lo;
    *hi = 1.0 / *hi;
}
/*}}}*/
int main (int argc, char **argv)/*{{{*/
{
#if 1
    int i;
    struct params cheby_params;
    struct params remez_params;
    double remez_bound;

    cheby_params.n = 5;
    cheby_params.ny = 5;

    do_fit(cheby_params.n, cheby_params.coef);
    for (i=0; i<cheby_params.n; i++) {
        printf("%2d : %15.8f\n", i, cheby_params.coef[i]);
    }
    fit_y_scaling(cheby_params.n, cheby_params.coef, cheby_params.ny, cheby_params.ycoef);
    for (i=0; i<cheby_params.ny; i++) {
        printf("%2d : %15.8f\n", i, cheby_params.ycoef[i]);
    }

    remez_params.n = 6;
    remez_params.ny = 6;
    remez_bound = 0.265;
    remez(-remez_bound, remez_bound, fn_for_remez, &remez_params, remez_params.n, remez_params.coef, 0);
#if 1
    remez(-remez_bound, remez_bound, fnyscale_for_remez, &remez_params, remez_params.ny, remez_params.ycoef, 0);
#endif
#if 0
    remez2(-0.25, 0.25, fny2scale_for_remez, &remez_params, remez_params.ny, remez_params.ycoef, 0);
#endif

    printf("Remez coefficients\n");
    for (i=0; i<remez_params.n; i++) {
        printf("%2d : %15.8f\n", i, remez_params.coef[i]);
    }
    for (i=0; i<remez_params.ny; i++) {
        printf("%2d : %15.8f\n", i, remez_params.ycoef[i]);
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
    run_data_with_model("geodesic_remez.dat", &remez_params, 0.0000025);
    run_data_with_model("geodesic_cheby.dat", &cheby_params, 0.0000025);

    return 0;
}
/*}}}*/
/* vim:et:sw=4:sts=4:ht=4
 * */

