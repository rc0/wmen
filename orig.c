#include <stdio.h>

int main (int argc, char **argv) {
  double E = 400000.0;
  double N = 650000.0;
  double X, Y;
  do {
    double x = (E - 400000.0) / 400000.0;
    double y = (N - 600000.0) / 600000.0;
    double x2 = x*x;
    double x3 = x2*x;
    double x4 = x3*x;
    double x5 = x4*x;
    double y2 = y*y;
    double y3 = y2*y;
    double y4 = y3*y;
    double y5 = y4*y;

    double X0 = +12361003.1 -17.9*y -1.9*y2;
    double X1 = +437429.8 +59288.7*y +9978.1*y2 +1536.3*y3 +247.0*y4 +38.5*y5;
    double X2 = -1.0 -17.7*y -21.4*y2 +3.0*y3;
    double X3 = -1481.1 -681.3*y -220.6*y2 -53.8*y3 -17.2*y4;
    double X4 = +4.1*y3 +4.1*y4;
    double X5 = +9.3 +5.8*y3 +7.5*y4;
    X = X0 +X1*x +X2*x2 +X3*x3 +X4*x4 +X5*x5;
    double Y0 = +7871705.5 -657579.0*y -44372.5*y2 -4992.5*y3 -586.4*y4 -73.5*y5;
    double Y1 = -13.4 -5.9*y +1.6*y2 +2.0*y3;
    double Y2 = +19803.9 +6645.0*y +1541.1*y2 +338.8*y3 +64.0*y4 +12.6*y5;
    double Y3 = -7.9 -15.1*y +17.8*y2 +26.1*y3;
    double Y4 = -112.9 -72.1*y -33.3*y2 -9.7*y3 +10.8*y4 +9.6*y5;
    double Y5 = -5.0*y2 -6.5*y3 +13.3*y4 +16.1*y5;
    Y = Y0 +Y1*x +Y2*x2 +Y3*x3 +Y4*x4 +Y5*x5;

  } while (0);

  X = X/25000000.0;
  Y = Y/25000000.0;
  printf("X=%.10f Y=%.10f\n", X, Y);
  return 0;
}

