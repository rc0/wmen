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

#include "contour.h"

struct node/*{{{*/
{
  struct node *next;
  struct node *prev;
  int valid;
  double value;
};
/*}}}*/
struct crow/*{{{*/
{
  struct crow *next;
  struct crow *prev;
  struct node nodes;
};
/*}}}*/
struct cdata/*{{{*/
{
  struct crow rows;
};
/*}}}*/
struct cdata *new_cdata(void)/*{{{*/
{
  struct cdata *result = malloc(sizeof(struct cdata));
  result->rows.next = result->rows.prev = &result->rows;
  return result;
}
/*}}}*/
static void free_row(struct crow *r)/*{{{*/
{
  struct node *n, *nn;
  for (n=r->nodes.next; n!=&r->nodes; n=nn) {
    nn = n->next;
    free(n);
  }
  free(r);
}
/*}}}*/
void free_cdata(struct cdata *c)/*{{{*/
{
  struct crow *r, *nr;
  for (r=c->rows.next; r!=&c->rows; r=nr) {
    nr = r->next;
    free_row(r);
  }
  free(c);
}
/*}}}*/
struct crow *next_crow(struct cdata *c)/*{{{*/
{
  struct crow *result;
  result = malloc(sizeof(struct crow));
  result->nodes.next = result->nodes.prev = &result->nodes;
  result->next = &c->rows;
  result->prev = c->rows.prev;
  c->rows.prev->next = result;
  c->rows.prev = result;
  return result;
}
/*}}}*/
static void insert_node(struct crow *row, struct node *nn)/*{{{*/
{
  nn->prev = row->nodes.prev;
  nn->next = &row->nodes;
  row->nodes.prev->next = nn;
  row->nodes.prev = nn;
}
/*}}}*/
void add_cnode(struct crow *row, double data)/*{{{*/
{
  struct node *nn;
  nn = malloc(sizeof(struct node));
  nn->valid = 1;
  nn->value = data;
  insert_node(row, nn);
}
/*}}}*/
void add_empty_cnode(struct crow *row)/*{{{*/
{
  struct node *nn;
  nn = malloc(sizeof(struct node));
  nn->valid = 0;
  nn->value = 0.0;
  insert_node(row, nn);
}
/*}}}*/

/* Do the endpoints exist? */
#define E_HAS0    (1<<0)
#define E_HAS1    (1<<1)
#define E_HASBOTH (E_HAS0 | E_HAS1)

#define E_HORIZ   (1<<2)
#define E_VOID    (1<<3)
#define E_ACTIVE  (1<<4)
#define E_START   (1<<5)

struct edge/*{{{*/
{
  int flags;
  double v0, v1; /* values at the 0, 1 ends */
  double frac; /* if E_CROSSED, how far along the crossing is */
};
/*}}}*/
static void generate_edges(const struct cdata *cd, double level, int *nx, int *ny, struct edge ***h, struct edge ***v)/*{{{*/
{
  struct crow *row;
  struct node *node;
  int tx, mx, my;
  int i, j;
  mx = my = 0;
  struct edge **H, **V;

  for (row=cd->rows.next; row!=&cd->rows; row=row->next) {
    my++;
    tx = 0;
    for (node=row->nodes.next; node!=&row->nodes; node=node->next) {
      tx++;
    }
    if (tx > mx) {
      mx = tx;
    }
  }
  /* Number of GRID SQUARES in each direction, not number of grid lines */
  mx--;
  my--;
  *nx = mx;
  *ny = my;

  /* Allocate edge arrays */
  H = malloc(mx * sizeof(struct edge *));
  for (i=0; i<mx; i++) {
    H[i] = malloc((my+1) * sizeof(struct edge));
    for (j=0; j<(my+1); j++) {
      H[i][j].flags = E_HORIZ;
    }
  }
  V = malloc((mx+1) * sizeof(struct edge *));
  for (i=0; i<(mx+1); i++) {
    V[i] = malloc(my * sizeof(struct edge));
    for (j=0; j<my; j++) {
      V[i][j].flags = 0;
    }
  }

  /* Build edges */
  for (i=0, row=cd->rows.next; row!=&cd->rows; i++, row=row->next) {
    for (j=0, node=row->nodes.next; node!=&row->nodes; j++, node=node->next) {
      if (node->valid) {
        if (j>0) {
          H[j-1][i].flags |= E_HAS1;
          H[j-1][i].v1 = node->value;
        }
        if (j<mx) {
          H[j][i].flags   |= E_HAS0;
          H[j][i].v0 = node->value;
        }
        if (i>0) {
          V[j][i-1].flags |= E_HAS1;
          V[j][i-1].v1 = node->value;
        }
        if (i<my) {
          V[j][i].flags |= E_HAS0;
          V[j][i].v0 = node->value;
        }
      }
    }
  }

  /* Find live edges and cut points */
  for (i=0; i<mx; i++) {
    for (j=0; j<(my+1); j++) {
      if ((H[i][j].flags & E_HASBOTH) == E_HASBOTH) {
        if ((H[i][j].v0 <= level) && (level < H[i][j].v1)) {
          H[i][j].flags |= E_ACTIVE;
          H[i][j].frac = (level - H[i][j].v0) / (H[i][j].v1 - H[i][j].v0);
        } else if ((H[i][j].v0 >= level) && (level > H[i][j].v1)) {
          H[i][j].flags |= E_ACTIVE;
          H[i][j].frac = (level - H[i][j].v0) / (H[i][j].v1 - H[i][j].v0);
        }
      } else {
        H[i][j].flags |= E_VOID;
      }
    }
  }
  for (i=0; i<(mx+1); i++) {
    for (j=0; j<my; j++) {
      if ((V[i][j].flags & E_HASBOTH) == E_HASBOTH) {
        if ((V[i][j].v0 <= level) && (level < V[i][j].v1)) {
          V[i][j].flags |= E_ACTIVE;
          V[i][j].frac = (level - V[i][j].v0) / (V[i][j].v1 - V[i][j].v0);
        } else if ((V[i][j].v0 >= level) && (level > V[i][j].v1)) {
          V[i][j].flags |= E_ACTIVE;
          V[i][j].frac = (level - V[i][j].v0) / (V[i][j].v1 - V[i][j].v0);
        }
      } else {
        V[i][j].flags |= E_VOID;
      }
    }
  }

  *h = H;
  *v = V;

}
/*}}}*/
enum direction/*{{{*/
{
  D_LEFT,
  D_RIGHT,
  D_UP,
  D_DOWN
};
/*}}}*/
static struct cpoint *new_cpoint(struct cpoint *last,/*{{{*/
    int i, int j, int is_horiz, double frac)
{
  struct cpoint *pt;
  pt = malloc(sizeof(struct cpoint));
  pt->next = last;
  pt->x = (double) i;
  pt->y = (double) j;
  if (is_horiz) pt->x += frac;
  else          pt->y += frac;
  return pt;
}
/*}}}*/

#define CHECK_V(x,y) (((x)>=0) && ((x)<=mx) && ((y)>=0) && ((y)<my) && ((v[x][y].flags & E_VOID) == 0) && ((v[x][y].flags & (E_ACTIVE|E_START)) != 0))
#define CHECK_H(x,y) (((x)>=0) && ((x)<mx) && ((y)>=0) && ((y)<=my) && ((h[x][y].flags & E_VOID) == 0) && ((h[x][y].flags & (E_ACTIVE|E_START)) != 0))
#define CLEAR_V(x,y) do { v[x][y].flags &= ~(E_ACTIVE|E_START);} while(0)
#define CLEAR_H(x,y) do { h[x][y].flags &= ~(E_ACTIVE|E_START);} while(0)

static struct cpoint *follow_path(int i, int j, int mx, int my,
    struct edge **h, struct edge **v,
    struct cpoint *pt,
    enum direction dir)/*{{{*/
{
  /* This always takes the first exit edge in a clockwise direction, i.e.
   * doesn't try to handle saddles intelligently. */
  do {
    switch (dir) {
      case D_UP: /* entering box from below */
        if (CHECK_V(i,j)) {
          pt = new_cpoint(pt, i, j, 0, v[i][j].frac);
          CLEAR_V(i,j);
          dir = D_LEFT;
        } else if (CHECK_H(i,j+1)) {
          j++;
          pt = new_cpoint(pt, i, j, 1, h[i][j].frac);
          CLEAR_H(i,j);
        } else if (CHECK_V(i+1,j)) {
          i++;
          pt = new_cpoint(pt, i, j, 0, v[i][j].frac);
          CLEAR_V(i,j);
          dir = D_RIGHT;
        } else {
          return pt;
        }
        break;
      case D_DOWN:
        if (CHECK_V(i+1,j-1)) {
          i++;
          j--;
          pt = new_cpoint(pt, i, j, 0, v[i][j].frac);
          CLEAR_V(i,j);
          dir = D_RIGHT;
        } else if (CHECK_H(i,j-1)) {
          j--;
          pt = new_cpoint(pt, i, j, 1, h[i][j].frac);
          CLEAR_H(i,j);
        } else if (CHECK_V(i,j-1)) {
          j--;
          pt = new_cpoint(pt, i, j, 0, v[i][j].frac);
          CLEAR_V(i,j);
          dir = D_LEFT;
        } else {
          return pt;
        }
        break;
      case D_LEFT:
        if (CHECK_H(i-1,j)) {
          i--;
          pt = new_cpoint(pt, i, j, 1, h[i][j].frac);
          CLEAR_H(i,j);
          dir = D_DOWN;
        } else if (CHECK_V(i-1,j)) {
          i--;
          pt = new_cpoint(pt, i, j, 0, v[i][j].frac);
          CLEAR_V(i,j);
        } else if (CHECK_H(i-1,j+1)) {
          i--;
          j++;
          pt = new_cpoint(pt, i, j, 1, h[i][j].frac);
          CLEAR_H(i,j);
          dir = D_UP;
        } else {
          return pt;
        }
        break;
      case D_RIGHT:
        if (CHECK_H(i,j+1)) {
          j++;
          pt = new_cpoint(pt, i, j, 1, h[i][j].frac);
          CLEAR_H(i,j);
          dir = D_UP;
        } else if (CHECK_V(i+1,j)) {
          i++;
          pt = new_cpoint(pt, i, j, 0, v[i][j].frac);
          CLEAR_V(i,j);
        } else if (CHECK_H(i,j)) {
          pt = new_cpoint(pt, i, j, 1, h[i][j].frac);
          CLEAR_H(i,j);
          dir = D_DOWN;
        } else {
          return pt;
        }
        break;
    }
  } while (1);
}
/*}}}*/
static struct cpoint *reverse(struct cpoint *list)/*{{{*/
{
  struct cpoint *p, *np, *result;
  for (p=list, result=NULL; p; p=np) {
    np = p->next;
    p->next = result;
    result = p;
    p = np;
  }
  return result;
}
/*}}}*/
static struct cline *add_cline(struct cline *last, struct cpoint *points)/*{{{*/
{
  struct cline *result;
  result = malloc(sizeof(struct cline));
  result->next = last;
  result->points = points;
  return result;
}
/*}}}*/
static struct cline *do_scan(int mx, int my, struct edge **h, struct edge **v)/*{{{*/
{
  struct cline *result = NULL;
  int i, j;
  for (i=0; i<mx; i++) {
    for (j=0; j<(my+1); j++) {
      if (h[i][j].flags & E_ACTIVE) {
        struct cpoint *pt;
        h[i][j].flags |= E_START;
        h[i][j].flags &= ~E_ACTIVE;
        pt = new_cpoint(NULL, i, j, 1, h[i][j].frac);
        pt = follow_path(i, j, mx, my, h, v, pt, D_UP);
        if (h[i][j].flags & E_START) {
          pt = reverse(pt);
          pt = follow_path(i, j, mx, my, h, v, pt, D_DOWN);
        }
        h[i][j].flags &= ~E_START;
        result = add_cline(result, pt);
      }
    }
  }
  for (i=0; i<(mx+1); i++) {
    for (j=0; j<my; j++) {
      if (v[i][j].flags & E_ACTIVE) {
        struct cpoint *pt;
        v[i][j].flags |= E_START;
        v[i][j].flags &= ~E_ACTIVE;
        pt = new_cpoint(NULL, i, j, 0, v[i][j].frac);
        pt = follow_path(i, j, mx, my, h, v, pt, D_RIGHT);
        if (v[i][j].flags & E_START) {
          pt = reverse(pt);
          pt = follow_path(i, j, mx, my, h, v, pt, D_LEFT);
        }
        v[i][j].flags &= ~E_START;
        result = add_cline(result, pt);
      }
    }
  }
  return result;
}
/*}}}*/
static void free_edges(int mx, int my, struct edge **h, struct edge **v)/*{{{*/
{
  int i;
  for (i=0; i<mx; i++) {
    free(h[i]);
  }
  free(h);
  for (i=0; i<(mx+1); i++) {
    free(v[i]);
  }
  free(v);
}
/*}}}*/
struct cline *generate_isolines(struct cdata *cd, double level)/*{{{*/
{
  struct cline *result;
  struct edge **h, **v;
  int mx, my;
  generate_edges(cd, level, &mx, &my, &h, &v);
  result = do_scan(mx, my, h, v);
  free_edges(mx, my, h, v);
  return result;
}
/*}}}*/
void free_clines(struct cline *x)/*{{{*/
{
  struct cline *nx;
  while (x!=NULL) {
    struct cpoint *p, *np;
    nx = x->next;
    p = x->points;
    while (p!=NULL) {
      np = p->next;
      free(p);
      p = np; /* :-) */
    }
    free(x);
    x = nx;
  }
}
/*}}}*/

#if 0
#define C0 "#000000"

#define C1 "#0000ff"
#define C2 "#00b0b0"
#define C3 "#00b000"
#define C4 "#b0b000"
#define C5 "#d07000"
#define C6 "#ff0000"
#define C7 "#b000b0"
#endif
#define C0 "black", "#000000"

#define C1 "purple!70!blue!80!white", "#8000c0"
#define C2 "blue", "#0000ff"
#define C3 "cyan!90!black", "#00c0c0"
#define C4 "green!70!black", "#00c000"
#define C5 "yellow!95!red!70!black", "#c0b000"
#define C6 "orange!80!black", "#e0a000"
#define C7 "red", "#ff0000"

#define THICK1 "thin", 1.0
#define THICK2 "thick", 1.4
#define THICK3 "very thick", 2.0
#define THICK4 "ultra thick", 2.8

/* Line data */
struct level levels_1side[N_LEVELS_1] =
{
  {C0,   0.00, THICK2},
  {C1,   0.10, THICK1},
  {C2,   0.13, THICK1},
  {C3,   0.18, THICK1},
  {C4,   0.25, THICK1},
  {C5,   0.35, THICK1},
  {C6,   0.50, THICK1},
  {C7,   0.70, THICK1},
  {C1,   1.00, THICK2},
  {C2,   1.30, THICK2},
  {C3,   1.80, THICK2},
  {C4,   2.50, THICK2},
  {C5,   3.50, THICK2},
  {C6,   5.00, THICK2},
  {C7,   7.00, THICK2},
  {C1,  10.00, THICK3},
  {C2,  13.00, THICK3},
  {C3,  18.00, THICK3},
  {C4,  25.00, THICK3},
  {C5,  35.00, THICK3},
  {C6,  50.00, THICK3},
  {C7,  70.00, THICK3},
  {C1, 100.00, THICK4},
  {C2, 130.00, THICK4},
  {C3, 180.00, THICK4},
  {C4, 250.00, THICK4},
};

struct level levels_2side[N_LEVELS_2] =
{
  {C7, -70.00, THICK3},
  {C6, -50.00, THICK3},
  {C5, -35.00, THICK3},
  {C4, -25.00, THICK3},
  {C3, -18.00, THICK3},
  {C2, -13.00, THICK3},
  {C1, -10.00, THICK3},
  {C7,  -7.00, THICK2},
  {C6,  -5.00, THICK2},
  {C5,  -3.50, THICK2},
  {C4,  -2.50, THICK2},
  {C3,  -1.80, THICK2},
  {C2,  -1.30, THICK2},
  {C1,  -1.00, THICK2},
  {C7,  -0.70, THICK1},
  {C6,  -0.50, THICK1},
  {C5,  -0.35, THICK1},
  {C4,  -0.25, THICK1},
  {C3,  -0.18, THICK1},
  {C2,  -0.13, THICK1},
  {C1,  -0.10, THICK1},
  {C0,   0.00, THICK2},
  {C1,   0.10, THICK1},
  {C2,   0.13, THICK1},
  {C3,   0.18, THICK1},
  {C4,   0.25, THICK1},
  {C5,   0.35, THICK1},
  {C6,   0.50, THICK1},
  {C7,   0.70, THICK1},
  {C1,   1.00, THICK2},
  {C2,   1.30, THICK2},
  {C3,   1.80, THICK2},
  {C4,   2.50, THICK2},
  {C5,   3.50, THICK2},
  {C6,   5.00, THICK2},
  {C7,   7.00, THICK2},
  {C1,  10.00, THICK3},
  {C2,  13.00, THICK3},
  {C3,  18.00, THICK3},
  {C4,  25.00, THICK3},
  {C5,  35.00, THICK3},
  {C6,  50.00, THICK3},
  {C7,  70.00, THICK3}
};

