/* Caister water tower example from OS paper */

#include <stdio.h>
#include <stdlib.h>
#include "tool.h"

int main (int argc, char **argv)
{
  struct llh etrs;
  double lon, lat;
  double P, Q;
  double odn;
  struct mxy mxy;
  double E, N;
  double X, Y;



  etrs.lon = lon = +1.0 + (42.0/60.0) + (57.8663/3600.0);
  etrs.lat = lat = +52.0 + (39.0/60.0) + (28.8282/3600.0);
  etrs.h = 108.05;

  P = 2.0*(etrs.lat - 54.5)/9.0;
  Q = (etrs.lon + 2.0)/4.0;

  odn = etrs.h + (-50.1 + 6.1*Q - 1.1*P - 1.5*P*Q);

  wgs84_to_mxy(&etrs, &mxy);

  do {
    double x = 63.0*(mxy.X - 0.4944401);
    double y = 37.0*(mxy.Y - 0.3126638);
    double x2 = x*x;
    double x3 = x2*x;
    double x4 = x3*x;
    double x5 = x4*x;
    double y2 = y*y;
    double y3 = y2*y;
    double y4 = y3*y;
    double y5 = y4*y;
    double E0 = +400000.86 -16.61*y -1.65*y2;
    double E1 = +358760.59 +50238.41*y +1856.79*y2 -215.05*y3 -36.24*y4;
    double E2 = +3.08 -17.33*y +13.16*y2 +8.99*y3;
    double E3 = -215.68 +68.76*y +35.39*y2 +5.66*y3;
    E = E0 +E1*x +E2*x2 +E3*x3;
    double N0 = +650001.19 -609558.71*y -42504.18*y2 -1018.04*y3 +96.32*y4 +12.16*y5;
    double N1 = -11.18 +4.35*y2 -3.84*y3;
    double N2 = +14783.67 +1106.27*y -188.21*y2 -46.05*y3 -2.88*y4;
    double N3 = -3.27 +3.48*y +9.38*y2 -1.56*y3;
    double N4 = +10.94 +6.83*y +1.96*y2 +2.85*y3 -3.21*y4;
    double N5 = -3.67*y4;
    N = N0 +N1*x +N2*x2 +N3*x3 +N4*x4 +N5*x5;
  } while (0);

  printf("ODN height : %f\n", odn);
  printf("Eastings   : %f\n", E);
  printf("Northings  : %f\n", N);
  printf("=========================\n");

  E = 651409.792;
  N = 313177.448;

  do {
    double e = (E - 400000.0) / 400000.0;
    double n = (N - 650000.0) / 650000.0;
    double e2 = e*e;
    double e3 = e2*e;
    double e4 = e3*e;
    double e5 = e4*e;
    double n2 = n*n;
    double n3 = n2*n;
    double n4 = n3*n;
    double n5 = n4*n;
    double n6 = n5*n;

    double X0 = +12361001.62 -19.66*n -2.24*n2;
    double X1 = +442441.04 +66066.30*n +12171.89*n2 +2057.38*n3 +364.45*n4 +62.11*n5;
    double X2 = -2.52 -23.19*n -23.76*n2 +3.72*n3;
    double X3 = -1540.23 -778.77*n -273.35*n2 -63.49*n3 -25.26*n4 -12.21*n5;
    double X4 = +7.76*n3 +6.93*n4;
    double X5 = +10.01 +12.09*n4 +10.61*n5;
    X = X0 +X1*e +X2*e2 +X3*e3 +X4*e4 +X5*e5;
    double Y0 = +7816596.43 -720504.09*n -53575.89*n2 -6596.44*n3 -830.91*n4 -117.01*n5 -15.94*n6;
    double Y1 = -13.88 -5.83*n +4.10*n2 +4.32*n3;
    double Y2 = +20368.69 +7485.47*n +1909.20*n2 +454.99*n3 +98.62*n4 +23.58*n5;
    double Y3 = -9.34 -12.10*n +26.34*n2 +23.54*n3 -5.75*n4;
    double Y4 = -119.10 -84.19*n -41.68*n2 -7.53*n3 +20.95*n4 +15.10*n5;
    double Y5 = +14.47*n4 +30.69*n5 +14.37*n6;
    Y = Y0 +Y1*e +Y2*e2 +Y3*e3 +Y4*e4 +Y5*e5;

  } while (0);

  mxy.X = X/25000000.0;
  mxy.Y = Y/25000000.0;
  mxy_to_wgs84(&mxy, &etrs);
  printf("lat : %.12f (orig=%.12f, diff=%.6e)\n", etrs.lat, lat, etrs.lat-lat);
  printf("lon : %.12f (orig=%.12f, diff=%.6e)\n", etrs.lon, lon, etrs.lon-lon);

  return 0;
}

