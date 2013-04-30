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

/* Use Vincenty's formulae to compute geodesic distance between 2 points
 * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tool.h"

static const double a = 6378137.0;
static const double b = 6356752.3141;
static const double f = 1.0 / 298.257223563;

/* Radians positions -> metres apart along geodesic */
double llh_rad_to_metres(const struct llh *p1, const struct llh *p2)
{
  double U1 = atan((1.0 - f) * tan(p1->lat));
  double U2 = atan((1.0 - f) * tan(p2->lat));
  double L, lambda, lambda0;
  double cos_2sigma_m;
  double cos2_alpha;
  double sigma;
  L = p2->lon - p1->lon;
  lambda = L;
  do {
    lambda0 = lambda;
    double t1 = cos(U2)*sin(lambda);
    double t2 = cos(U1)*sin(U2) - sin(U1)*cos(U2)*cos(lambda);
    double sin_sigma = sqrt(t1*t1 + t2*t2);
    double cos_sigma = sin(U1)*sin(U2) + cos(U1)*cos(U2)*cos(lambda);
    sigma = atan2(sin_sigma, cos_sigma);
    /* is this meant to use sin_sigma or a freshly computed sin(sigma) ? */
    double sin_alpha = cos(U1)*cos(U2)*sin(lambda)/sin(sigma);
    double sin2_alpha = sin_alpha * sin_alpha;
    cos2_alpha = 1.0 - sin2_alpha;
    cos_2sigma_m = cos(sigma) - 2.0*sin(U1)*sin(U2)/cos2_alpha;
    double C = (1.0/16.0) * f * cos2_alpha * (4.0 + f * (4.0 - 3.0*cos2_alpha));
    double squares = cos_2sigma_m + C * cos(sigma) * (-1.0 + 2.0 * cos_2sigma_m * cos_2sigma_m);
    double braces = sigma + C * sin(sigma) * squares;
    lambda = L + (1.0 - C) * f * sin_alpha * braces;
  } while (fabs(lambda0 - lambda) > 1.0e-12);

  double u2 = cos2_alpha * (a - b) * (a + b) / (b*b);
  double A = 1.0 + (u2 / 16384.0) * (
      4096.0 + u2 * (-768.0 + u2*(320.0 - 175.0*u2))
      );
  double B = (u2 / 1024.0) * (
      256.0 + u2*(-128.0 + u2*(74.0 - 47.0*u2))
      );
  double z1 = cos(sigma)*(-1.0 + 2.0*cos_2sigma_m*cos_2sigma_m);
  double z2 = -3.0 + 4.0 * sin(sigma)*sin(sigma);
  double z3 = -3.0 + 4.0 * cos_2sigma_m * cos_2sigma_m;
  double z4 = (B/6.0) * cos_2sigma_m * z2 * z3;
  double DELTA_sigma = B * sin(sigma) * (
      cos_2sigma_m + (B/4.0) * (z1 - z4)
      );
  double s = b * A * (sigma - DELTA_sigma);
  return s;

}

double llh_deg_to_metres(const struct llh *p1, const struct llh *p2)
{
  struct llh q1;
  struct llh q2;
  q1.lat = (M_PI/180.0) * p1->lat;
  q1.lon = (M_PI/180.0) * p1->lon;
  q2.lat = (M_PI/180.0) * p2->lat;
  q2.lon = (M_PI/180.0) * p2->lon;
  return llh_rad_to_metres(&q1, &q2);
}

