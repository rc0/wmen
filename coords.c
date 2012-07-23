/*
* Copyright (c) 2012, Richard P. Curnow
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

/* Convert WGS84 to OSGB36 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tool.h"

struct param
{
  double a0, b0;
  double tx, ty, tz;
  double rx, ry, rz;
  double s;
  double a1, b1;
};

#define WGS84_A 6378137.0
#define WGS84_B 6356752.3142

#define AIRY1830_A 6377563.396
#define AIRY1830_B 6356256.910

static struct param wgs84_to_osgb36_param =
{
  WGS84_A,
  WGS84_B,
  -446.448,
  125.157,
  -542.060,
  -0.1502,
  -0.2470,
  -0.8421,
  20.4894,
  AIRY1830_A,
  AIRY1830_B
};

static struct param osgb36_to_wgs84_param =
{
  AIRY1830_A,
  AIRY1830_B,
  446.448,
  -125.157,
  542.060,
  0.1502,
  0.2470,
  0.8421,
  -20.4894,
  WGS84_A,
  WGS84_B
};

static struct param itrs2005_to_etrs89_param =
{
  WGS84_A,
  WGS84_B,
  0.056,
  0.048,
  -0.037,
  /* THE FOLLOWING MUST BE SCALED BY THE NUMBER OF YEARS!! */
  0.000054,
  0.000518,
  -0.000781,

  0.0, /* s */
  WGS84_A,
  WGS84_B
};

static inline double radof(double deg)
{
  return deg * (M_PI/180.0);
}

static inline double degof(double rad)
{
  return rad * (180.0/M_PI);
}

static void convert(const struct llh *in, struct llh *out, const struct param *pp)
{
  double lat, lon;
  double a, b;
  double sinlat, coslat, sinlon, coslon;
  double esq;
  double nu;
  double x1, y1, z1;
  double tx, ty, tz;
  double rx, ry, rz;
  double s1;
  double x2, y2, z2;
  double precision;
  double p, phi, lambda, phiP, H;

  lat = radof(in->lat);
  lon = radof(in->lon);
  sinlat = sin(lat);
  coslat = cos(lat);
  sinlon = sin(lon);
  coslon = cos(lon);

  a = pp->a0;
  b = pp->b0;
  esq = (a*a - b*b)/(a*a);

  nu = a / sqrt(1.0 - esq*sinlat*sinlat);
  x1 = (nu + in->h)  * coslat * coslon;
  y1 = (nu + in->h)  * coslat * sinlon;
  z1 = ((1.0-esq)*nu + in->h) * sinlat;

  tx = pp->tx;
  ty = pp->ty;
  tz = pp->tz;
  rx = radof(pp->rx/3600.0);
  ry = radof(pp->ry/3600.0);
  rz = radof(pp->rz/3600.0);
  s1 = 1.0 + (pp->s / 1.0e6);

  x2 = tx + x1*s1 - y1*rz + z1*ry;
  y2 = ty + x1*rz + y1*s1 - z1*rx;
  z2 = tz - x1*ry + y1*rx + z1*s1;

  a = pp->a1;
  b = pp->b1;
  precision = 4.0 / a;
  esq = (a*a - b*b)/(a*a);
  p = sqrt(x2*x2 + y2*y2);
  phi = atan2(z2, p*(1.0-esq));
  phiP = 2.0*M_PI;
  while(abs(phi-phiP) > precision) {
    double sinphi;
    sinphi = sin(phi);
    nu = a/sqrt(1.0 - esq*sinphi*sinphi);
    phiP = phi;
    phi = atan2(z2 + esq*nu*sinphi, p);
  }
  lambda = atan2(y2, x2);
  H = p/cos(phi) - nu;
  out->lat = degof(phi);
  out->lon = degof(lambda);
  out->h = H;
}

void wgs84_to_osgb36(const struct llh *in, struct llh *out)
{
  convert(in, out, &wgs84_to_osgb36_param);
}

void osgb36_to_wgs84(const struct llh *in, struct llh *out)
{
  convert(in, out, &osgb36_to_wgs84_param);
}

void itrs2005_to_etrs89(const struct llh *in, double year, struct llh *out)
{
  struct param temp;
  double dy = year - 1989.0;
  temp = itrs2005_to_etrs89_param;
  temp.rx *= dy;
  temp.ry *= dy;
  temp.rz *= dy;
  convert(in, out, &temp);
}

void foo_to_grid(const struct llh *in,
    double a, double b,
    struct en *out)
{
  double lat, lon;
  double F0;
  double lat0, lon0;
  double N0, E0, e2, n, n2, n3;
  double coslat, sinlat, nu, rho, eta2;
  double Ma, Mb, Mc, Md, M;
  double cos3lat, cos5lat, tanlat, tan2lat, tan4lat;
  double I, II, III, IIIA, IV, V, VI;
  double dLon, dLon2, dLon3, dLon4, dLon5, dLon6;

  lat = radof(in->lat);
  lon = radof(in->lon);
  F0 = 0.9996012717;
  lat0 = radof(49.0);
  lon0 = radof(-2.0);
  N0 = -100000.0;
  E0 = 400000.0;
  e2 = 1.0 - (b*b)/(a*a);
  n = (a-b)/(a+b);
  n2 = n*n;
  n3 = n*n2;

  coslat = cos(lat);
  sinlat = sin(lat);
  nu = a*F0/sqrt(1.0 - e2*sinlat*sinlat);
  rho = a*F0*(1.0-e2)/pow(1.0-e2*sinlat*sinlat, 1.5);
  eta2 = nu/rho - 1.0;

  Ma = (1.0 + n + (5.0/4.0)*n2 + (5.0/4.0)*n3) * (lat - lat0);
  Mb = (3.0*n + 3.0*n2 + (21.0/8.0)*n3) * sin(lat-lat0) * cos(lat+lat0);
  Mc = ((15.0/8.0)*n2 + (15.0/8.0)*n3) * sin(2.0*(lat-lat0)) * cos(2.0*(lat+lat0));
  Md = (35.0/24.0)*n3 * sin(3.0*(lat-lat0)) * cos(3.0*(lat+lat0));
  M = b * F0 * (Ma - Mb + Mc - Md);

  cos3lat = coslat*coslat*coslat;
  cos5lat = cos3lat*coslat*coslat;
  tanlat = tan(lat);
  tan2lat = tanlat * tanlat;
  tan4lat = tan2lat * tan2lat;

  I = M + N0;
  II = (nu/2.0)*sinlat*coslat;
  III = (nu/24.0)*sinlat*cos3lat*(5.0 - tan2lat + 9*eta2);
  IIIA = (nu/720.0)*sinlat*cos5lat*(61.0-58.0*tan2lat+tan4lat);
  IV = nu*coslat;
  V = (nu/6.0)*cos3lat*(nu/rho - tan2lat);
  VI = (nu/120.0) * cos5lat * (5.0 - 18.0*tan2lat + tan4lat + 14.0*eta2 - 58.0*tan2lat*eta2);

  dLon = lon-lon0;
  dLon2 = dLon*dLon;
  dLon3 = dLon*dLon2;
  dLon4 = dLon*dLon3;
  dLon5 = dLon*dLon4;
  dLon6 = dLon*dLon5;

  out->N = I + II*dLon2 + III*dLon4 + IIIA*dLon6;
  out->E = E0 + IV*dLon + V*dLon3 + VI*dLon5;
}

void osgb36_to_grid(const struct llh *in, struct en *out)
{
  foo_to_grid(in, AIRY1830_A, AIRY1830_B, out);
}

void wgs84_to_grid(const struct llh *in, struct en *out)
{
  foo_to_grid(in, WGS84_A, WGS84_B, out);
}

void wgs84_to_mxy(const struct llh *in, struct mxy *out)
{
  double x, yy, y;
  x = radof(in->lon);
  yy = radof(in->lat);
  y = log(tan(yy) + 1.0/cos(yy));
  out->X = 0.5 * (1.0 + x/M_PI);
  out->Y = 0.5 * (1.0 - y/M_PI);
}

void mxy_to_wgs84(const struct mxy *in, struct llh *out)
{
  double t;
  out->lon = degof(M_PI * ((2.0 * in->X)-1.0));
  t = M_PI * (1.0 - (2.0 * in->Y));
  t = 2.0 * atan(exp(t)) - (0.5*M_PI);
  out->lat = degof(t);
  out->h = 0.0;
}

#if TEST

static void do_test(double lat, double lon, double h)
{
  struct llh wgs, osgb;
  struct en en;
  wgs.lat = lat;
  wgs.lon = lon;
  wgs.h = h;
  wgs84_to_osgb36(&wgs, &osgb);
  osgb36_to_grid(&osgb, &en);
  printf("lat=%f, lon=%f, h=%f, E=%f, N=%f\n",
      osgb.lat,
      osgb.lon,
      osgb.h, en.E, en.N);
  return;
}


int main(int argc, char **argv)
{
  struct llh wgs, osgb;
  struct en en;
  do_test(53.0, 1.0, 24.7);
  do_test(51.36935, -2.92031, 60.0);
  do_test(52.6579766111, 1.716038416666, 50.0);
  return 0;

}
#endif

#if MAP

struct model
{


};

static void fit(struct model *result)
{
  struct llh osgb_orig, wgs84_orig;
  struct mxy mxy_orig;
  double lat, lon;
  double Sx2, Sx, n;
  double Sy2, Sy, SNy, SN;
  double R;
  double C, H;
  osgb_orig.lat = 52.0 /* 49.0 */;
  osgb_orig.lon = -2.0;
  osgb_orig.h = 0.0;
  osgb36_to_wgs84(&osgb_orig, &wgs84_orig);
  wgs84_to_mxy(&wgs84_orig, &mxy_orig);

  Sy2 = 0.0;
  SNy = 0.0;
  SN = 0.0;
  Sy = 0.0;
  n = 0.0;

  /* TODO : make the upper bounds deterministic */
  for (lat=50.0; lat<=59.0; lat+=0.5) {
    for (lon=-6.0; lon<=2.0; lon+=0.5) {
      struct llh osgb, wgs;
      struct mxy xy;
      struct en en;
      double xi, yi;
      osgb.lat = lat;
      osgb.lon = lon;
      osgb36_to_wgs84(&osgb, &wgs);
      wgs84_to_mxy(&wgs, &xy);
      osgb36_to_grid(&osgb, &en);
      xi = xy.X - mxy_orig.X;
      yi = xy.Y - mxy_orig.Y;

      Sy2 += yi*yi;
      Sy += yi;
      SNy += en.N * yi;
      SN += en.N;
      n += 1.0;
    }
  }

  printf("%f points\n", n);

  R = Sy2 / Sy;
  H = (SNy - R*SN) / (Sy - R*n);

  R = Sy / n;
  C = (SNy - R*SN) / (Sy2 - R*Sy);

  for (lat=50.0; lat<=59.0; lat+=0.5) {
    for (lon=-6.0; lon<=2.0; lon+=0.5) {
      struct llh osgb, wgs;
      struct mxy xy;
      struct en en;
      double ee, ne;
      double y0;
      double err;
      osgb.lat = lat;
      osgb.lon = lon;
      osgb36_to_wgs84(&osgb, &wgs);
      wgs84_to_mxy(&wgs, &xy);
      osgb36_to_grid(&osgb, &en);

      y0 = xy.Y - mxy_orig.Y;
      ne = H + C*y0;
      err = ne - en.N;
      printf("%f %f %f\n", lat, lon, err);
    }
  }

}

int main(int argc, char **argv)
{
  struct llh osgb_orig, wgs84_orig;
  struct mxy mxy_orig;
  struct model model;
  osgb_orig.lat = 49.0;
  osgb_orig.lon = -2.0;
  osgb_orig.h = 0.0;
  osgb36_to_wgs84(&osgb_orig, &wgs84_orig);
  wgs84_to_mxy(&wgs84_orig, &mxy_orig);
  printf("%15.10f %15.10f\n", mxy_orig.X, mxy_orig.Y);

  fit(&model);

  return 0;
}
#endif


