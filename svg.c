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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "contour.h"

static FILE *out;
/* pngs are loaded into inkscape at 1 pixel per point and each tile is 256
 * pixels in the png */
static double scale = 256.0 ;
static double stroke_width = 1.0;

void start_svg(char *name)
{
  out = fopen(name, "w");

  fprintf(out, "<?xml version=\"1.0\"?>\n");
  fprintf(out, "<svg\n");

  fprintf(out, "xmlns:svg=\"http://www.w3.org/2000/svg\"\n"
         "xmlns=\"http://www.w3.org/2000/svg\"\n"
         "xmlns:xlink=\"http://www.w3.org/1999/xlink\"\n"
         "id=\"svg2\"\n"
         "height=\"1052.3622\"\n"
         "width=\"744.09448\"\n"
         "y=\"0.0000000\"\n"
         "x=\"0.0000000\"\n"
         "version=\"1.0\">\n");
  fprintf(out, "<defs\n"
         "id=\"defs3\" />\n");

#if 0
  fprintf(out, "<image x=\"0\" y=\"0\" xlink:href=\"file:///mnt/flash/optima/wgs/osm/uk12.png\" height=\"1024\" width=\"768\" />\n");
#endif
  fprintf(out, "<image x=\"0\" y=\"0\" xlink:href=\"osm/uk12.png\" height=\"1024\" width=\"768\" />\n");
  fprintf(out, "<g id=\"layer1\">\n");
}

static inline double mapx(double x)
{
  return (x - (30.0/64.0)) * 64.0;
}

static inline double mapy(double y)
{
  return (y - (18.0/64.0)) * 64.0;
}

void emit_svg(double level, const char *colour, double width_scale, struct cline *lines)
{
  struct cline *line;

  for (line=lines; line; line=line->next) {
    struct cpoint *p;
    struct cpoint *fp, *lp;
    double minx, miny, maxx, maxy;
    double mhd;
    int first;
    fprintf(out, "<path style=\"fill:none;stroke:%s;stroke-width:%f;stroke-linecap:square;stroke-linejoin:miter;stroke-miterlimit:4.0;stroke-opacity:1.0\"\n",
        colour,
        stroke_width * width_scale);

    first = 1;
    fp = line->points;
    for (p=line->points; p; p=p->next) {
      double x = mapx(p->x);
      double y = mapy(p->y);
      if (first) {
        fprintf(out, "d=\"M %f,%f",
            scale * x, scale * y);
        first = 0;
        minx = maxx = x;
        miny = maxy = y;
      } else {
        fprintf(out, " L %f,%f",
            scale * x,
            scale * y);
        if (x < minx) minx = x;
        if (x > maxx) maxx = x;
        if (y < miny) miny = y;
        if (y > maxy) maxy = y;
      }
      lp = p;
    }
    fprintf(out, "\" />\n");

    mhd = (maxx - minx) + (maxy - miny);
    if (mhd > 0.05) {
      double xc, yc;
      xc = mapx(fp->x);
      yc = mapy(fp->y);
      /* Draw caption on the line */
      fprintf(out, "<text style=\"font-size:10;font-style:normal;font-variant:normal;"
               "font-weight:bold;fill:%s;fill-opacity:1.0;stroke:none;"
               "font-family:Luxi Sans;text-anchor:left;writing-mode:lr-tb\"\n",
               colour);
      if (level < 1.0) {
        fprintf(out, "x=\"%f\" y=\"%f\">%.2f</text>\n",
            scale*xc,
            scale*yc, level);
      } else if (level < 10.0) {
        fprintf(out, "x=\"%f\" y=\"%f\">%.1f</text>\n",
            scale*xc,
            scale*yc, level);
      } else {
        fprintf(out, "x=\"%f\" y=\"%f\">%.0f</text>\n",
            scale*xc,
            scale*yc, level);
      }
    }

  }


  return;
}

void emit_reg_marks(double xoff)
{
  int i, j;
  double x, y;
  double half;
  half = scale / 32.0;
  for (i=0; i<=3; i++) {
    for (j=0; j<=4; j++) {
      x = xoff + scale * (double) i;
      y = scale * (double) j;
      fprintf(out, "<path style=\"fill:none;stroke:%s;stroke-width:%f;stroke-linecap:square;stroke-linejoin:miter;stroke-miterlimit:4.0;stroke-opacity:1.0\"\n",
        "#000",
        0.5 * stroke_width);
      fprintf(out, "d=\"M %f,%f L %f,%f M %f,%f L %f,%f\" />\n",
          x - half, y,
          x + half, y,
          x, y - half,
          x, y + half
          );
    }
  }
}

void finish_svg(void)
{
  fprintf(out, "</g>\n");
  fprintf(out, "</svg>\n");

}
