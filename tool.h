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

#ifndef TOOL_H
#define TOOL_H

struct llh
{
  double lat;
  double lon;
  double h;
};

struct en
{
  double E;
  double N;
};

struct mxy
{
  /* As fraction */
  double X;
  double Y;
};

extern void wgs84_to_osgb36(const struct llh *in, struct llh *out);
extern void osgb36_to_wgs84(const struct llh *in, struct llh *out);
extern void itrs2005_to_etrs89(const struct llh *in, double year, struct llh *out);
extern void osgb36_to_grid(const struct llh *in, struct en *out);
extern void wgs84_to_grid(const struct llh *in, struct en *out);

/* Degrees lat/lon */
extern void wgs84_to_mxy(const struct llh *in, struct mxy *out);
extern void mxy_to_wgs84(const struct mxy *in, struct llh *out);

extern void uk_range(double lat, double *lon0, double *lon1);

#define SZ 100
typedef long double Matrix[SZ][SZ];
typedef char * SMatrix[SZ][SZ];

extern void solve(Matrix m, int n, long double *l, long double *r);
extern void solve2(Matrix m, int n, long double *l0, long double *l1, long double *r0, long double *r1);

extern void estrin(Matrix m, int n, const char *vr, const char *vc, char *res);
extern void estrin2(SMatrix m, int nx, int ny, const char *vr, const char *vc, const char *res);

extern double llh_rad_to_metres(const struct llh *p1, const struct llh *p2);
extern double llh_deg_to_metres(const struct llh *p1, const struct llh *p2);

#endif
