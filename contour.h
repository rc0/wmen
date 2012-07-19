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

#ifndef CONTOUR_H
#define CONTOUR_H

struct crow;
struct cdata;

struct cpoint
{
  struct cpoint *next;
  double x, y;
};

struct cline
{
  struct cline *next;
  struct cpoint *points;
};

extern struct cdata *new_cdata(void);
extern void free_cdata(struct cdata *c);
extern struct crow *next_crow(struct cdata *c);
extern void add_cnode(struct crow *row, double data);
extern void add_empty_cnode(struct crow *row);

extern struct cline *generate_isolines(struct cdata *cd, double level);
extern void free_clines(struct cline *x);

extern void start_svg(char *name);
extern void emit_svg(double level, const char *colour, double width_scale, struct cline *lines);
extern void emit_reg_marks(double xoff);
extern void finish_svg(void);

#define REF6 64.0

extern void en_to_wm(double E, double N, double *X, double *Y);
extern void XY_to_xy(double X, double Y, double *x, double *y);
extern void start_tikz(char *name);
extern void emit_tikz(double level, const char *colour, const char *thickness, struct cline *lines);
extern void draw_map_tikz(FILE *out, const struct cline *map);
extern void draw_map_mxy_tikz(FILE *out, const struct cline *map);
extern struct cline *load_map(const char *filename);
extern void finish_tikz(void);
extern void generate_contour_key_tikz(FILE *out);

struct level {/*{{{*/
  char *colour;
  char *svg_colour;
  double level;
  char *thickness;
  double width_scale;
};
/*}}}*/
#define N_LEVELS_1 26
struct level levels_1side[N_LEVELS_1];
#define N_LEVELS_2 43
struct level levels_2side[N_LEVELS_2];

#endif
