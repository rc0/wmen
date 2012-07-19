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

/* Implement Estrin's scheme for computing a 2D polynomial (in this case, with
 * some zero coefficients) */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tool.h"

enum type/*{{{*/
{
  X_CONST,
  X_MACMUL,
};
/*}}}*/
struct tree/*{{{*/
{
  enum type type;

  int number;

  char *value;

  struct tree *t0;

  int horiz;
  int pow;
  struct tree *t1;

};
/*}}}*/
static int get_step(int n)/*{{{*/
{
  int count=0;
  while (n > 0) {
    count++;
    n >>= 1;
  }
  count--;
  return 1<<count;
}
/*}}}*/
static struct tree *new_const(const char *value)/*{{{*/
{
  struct tree *result;
  result = malloc(sizeof(struct tree));
  result->type = X_CONST;
  result->value = strdup(value);
  return result;
}
/*}}}*/
static struct tree *new_mult(struct tree *t, int horiz, int pow)/*{{{*/
{
  struct tree *result;
  result = malloc(sizeof(struct tree));
  result->type = X_MACMUL;
  result->horiz = horiz;
  result->pow = pow;
  result->t1 = t;
  result->t0 = NULL;
  return result;
}
/*}}}*/
static struct tree *new_mac(struct tree *t0, struct tree *t1, int horiz, int pow)/*{{{*/
{
  struct tree *result;
  result = malloc(sizeof(struct tree));
  result->type = X_MACMUL;
  result->horiz = horiz;
  result->pow = pow;
  result->t0 = t0;
  result->t1 = t1;
  return result;
}
/*}}}*/
static void inner_number_temps(struct tree *cursor, int *next)/*{{{*/
{
  if (cursor->type == X_MACMUL) {
    if (cursor->t0) {
      inner_number_temps(cursor->t0, next);
    }
    if (cursor->t1) {
      inner_number_temps(cursor->t1, next);
    }

    cursor->number = (*next)++;
  }
}
/*}}}*/
static void write_decls(int start, int finish)/*{{{*/
{
  int i;
  int j;
  int perline = 8;
  for (i=0; i+start < finish; i += perline) {
    int limit = i + perline;
    if (limit + start > finish) {
      limit = finish - start;
    }
    printf("double ");
    for (j=i; j<limit; j++) {
      if (j>i) printf(", ");
      printf("t%d", j+start);
    }
    printf(";\n");
  }
}
/*}}}*/
static void number_temps(struct tree *top)/*{{{*/
{
  static int next = 0;
  int start, finish;
  start = next;
  if (top) {
    if (top->type == X_MACMUL) {
      if (top->t0) inner_number_temps(top->t0, &next);
      inner_number_temps(top->t1, &next);
    }
  }
  finish = next;
  write_decls(start, finish);
}
/*}}}*/
static void inner_print_tree(struct tree *cursor, const char *res, const char *vr, const char *vc)/*{{{*/
{
  if (cursor->t0 && cursor->t0->type == X_MACMUL) {
    inner_print_tree(cursor->t0, NULL, vr, vc);
  }
  if (cursor->t1 && cursor->t1->type == X_MACMUL) {
    inner_print_tree(cursor->t1, NULL, vr, vc);
  }
  if (res) {
    printf("%s = ", res);
  } else {
    printf("t%d = ", cursor->number);
  }
  if (cursor->t0) {
    switch(cursor->t0->type) {
      case X_CONST:
        printf("(%s)", cursor->t0->value);
        break;
      case X_MACMUL:
        printf("t%d", cursor->t0->number);
        break;
    }
    if (cursor->t1) printf(" + ");
  }
  if (cursor->t1) {
    switch(cursor->t1->type) {
      case X_CONST:
        printf("(%s)", cursor->t1->value);
        break;
      case X_MACMUL:
        printf("t%d", cursor->t1->number);
        break;
    }
    printf("*%s", cursor->horiz ? vr : vc);
    if (cursor->pow > 1) {
      printf("%d", cursor->pow);
    }
  }
  printf(";\n");
}
/*}}}*/
static void print_tree(struct tree *top, const char *res, const char *vr, const char *vc)/*{{{*/
{
  inner_print_tree(top, res, vr, vc);
}
/*}}}*/

static int count_mul(const struct tree *t)/*{{{*/
{
  int total;
  total = 0;
  if (t) {
    if (t->type == X_MACMUL) {
      total += 1;
      if (t->t0) {
        total += count_mul(t->t0);
      }
      total += count_mul(t->t1); /* MUST BE NON-NULL else error! */
    }
  }
  return total;
}
/*}}}*/
static void free_tree(struct tree *t)/*{{{*/
{
  if (t) {
    if (t->type == X_MACMUL) {
      if (t->t0) {
        free_tree(t->t0);
      }
      free_tree(t->t1);
    } else if (t->type == X_CONST) {
      if (t->value) free(t->value);
    }
    free(t);
  }
}
/*}}}*/
static void find_max_pow(const struct tree *cursor, int *max_r, int *max_c)/*{{{*/
{
  if (cursor) {
    if (cursor->type == X_MACMUL) {
      if (cursor->horiz) {
        if (cursor->pow > *max_r) *max_r = cursor->pow;
      } else {
        if (cursor->pow > *max_c) *max_c = cursor->pow;
      }
      if (cursor->t0) {
        find_max_pow(cursor->t0, max_r, max_c);
      }
      find_max_pow(cursor->t1, max_r, max_c);
    }
  }
}
/*}}}*/
static void emit_powers(const struct tree *top, const char *vr, const char *vc)/*{{{*/
{
  int max_r, max_c;
  int i;
  max_r = 0;
  max_c = 0;
  find_max_pow(top, &max_r, &max_c);
  for (i=1; i<max_r; i<<=1) {
    printf("double %s%d = ", vr, i<<1);
    if (i > 1) {
      printf("%s%d * %s%d;\n", vr, i, vr, i);
    } else {
      printf("%s * %s;\n", vr, vr);
    }
  }
  for (i=1; i<max_c; i<<=1) {
    printf("double %s%d = ", vc, i<<1);
    if (i > 1) {
      printf("%s%d * %s%d;\n", vc, i, vc, i);
    } else {
      printf("%s * %s;\n", vc, vc);
    }
  }
}
/*}}}*/
static struct tree *generate(SMatrix m,/*{{{*/
    int nx, int ny,
    int bx, int by,
    int sx, int sy)
{
  if ((bx >= nx) || (by >= ny)) return NULL;
  if ((sx >= 1) && (sy >= 1)) {

    struct tree *t0x, *t1x, *t0y, *t1y;
    int count_x, count_y;
    t0x = generate(m, nx, ny, bx, by, sx>>1, sy);
    t1x = generate(m, nx, ny, bx+sx, by, sx>>1, sy);
    t0y = generate(m, nx, ny, bx, by, sx, sy>>1);
    t1y = generate(m, nx, ny, bx, by+sy, sx, sy>>1);

    count_x = count_mul(t0x) + count_mul(t1x);
    count_y = count_mul(t0y) + count_mul(t1y);

#if 0
    /* Force xy-first */
    count_x = -1;
#endif

    if (count_x <= count_y) {
      free_tree(t0y);
      free_tree(t1y);
      if (t0x && t1x) {
        return new_mac(t0x, t1x, 1, sx);
      } else if (t0x) {
        return t0x;
      } else if (t1x) {
        return new_mult(t1x, 1, sx);
      }
    } else {
      free_tree(t0x);
      free_tree(t1x);
      if (t0y && t1y) {
        return new_mac(t0y, t1y, 0, sy);
      } else if (t0y) {
        return t0y;
      } else if (t1y) {
        return new_mult(t1y, 0, sy);
      }
    }

  } else if (sy >= 1) {
    struct tree *t0, *t1;
    t0 = generate(m, nx, ny, bx, by, sx, sy>>1);
    t1 = generate(m, nx, ny, bx, by+sy, sx, sy>>1);
    if (t0 && t1) {
      return new_mac(t0, t1, 0, sy);
    } else if (t0) {
      return t0;
    } else if (t1) {
      return new_mult(t1, 0, sy);
    }

  } else if (sx >= 1) {
    /* do both */
    struct tree *t0, *t1;
    t0 = generate(m, nx, ny, bx, by, sx>>1, sy);
    t1 = generate(m, nx, ny, bx+sx, by, sx>>1, sy);
    if (t0 && t1) {
      return new_mac(t0, t1, 1, sx);
    } else if (t0) {
      return t0;
    } else if (t1) {
      return new_mult(t1, 1, sx);
    }
  } else {
    /* Just looking at a single square */
    if (m[bx][by]) {
      return new_const(m[bx][by]);
    } else {
      return NULL;
    }
  }
}
/*}}}*/
void estrin2(SMatrix m, int nx, int ny, const char *vr, const char *vc, const char *res)/*{{{*/
{
  struct tree *top;
  int sx, sy;

  sx = get_step(nx);
  sy = get_step(ny);

  top = generate(m, nx, ny, 0, 0, sx, sy);
  number_temps(top);
  emit_powers(top, vr, vc);
  print_tree(top, res, vr, vc);
  free_tree(top);
  return;
}
/*}}}*/
#if defined (TEST)
int main (int argc, char **argv)/*{{{*/
{
  SMatrix foo;
  int nx, ny;
  int i, j;

  srand48(12345);

  memset(foo, 0, sizeof(foo));
  nx = 7;
  ny = 17;
  for (i=0; i<nx; i++) {
    for (j=0; j<ny; j++) {
      if (drand48() < 0.2) {
        char buffer[32];
        sprintf (buffer, "%.10f", drand48() * 16.0);
        foo[i][j] = strdup(buffer);
        printf("foo[%d][%d] = %s\n", i, j, foo[i][j]);
      } else {
        foo[i][j] = NULL;
      }
    }
  }

  estrin2(foo, nx, ny, "P", "Q", "XYZ");
  return 0;
}
/*}}}*/
#endif

