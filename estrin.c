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

static int get_step(int n)
{
  int count=0;
  while (n > 0) {
    count++;
    n >>= 1;
  }
  count--;
  return 1<<count;
}

static void gather(Matrix m, int n, int i, int j,
    int horiz, int step,
    int *next,
    const char *vr, const char *vc, char **res)
{
  if (step > 0) {
    if ((horiz ? (i + step) : (j + step)) < n) {
      char *x1, *x2;
      x1 = x2 = NULL;
      if (horiz) {
        gather(m, n, i,      j, 1-horiz, step>>1, next, vr, vc, &x1);
        gather(m, n, i+step, j, 1-horiz, step>>1, next, vr, vc, &x2);
      } else {
        gather(m, n, i,      j, 1-horiz, step   , next, vr, vc, &x1);
        gather(m, n, i, j+step, 1-horiz, step   , next, vr, vc, &x2);
      }
      if (x1 && x2) {
        if (!*res) {
          char buffer[16];
          sprintf(buffer, "t%d", *next);
          ++*next;
          *res = strdup(buffer);
        }
        if (step > 1) {
          printf("%s = %s + (%s * %s%d);\n", *res, x1, x2,
              horiz ? vc : vr,
              step);
        } else {
          printf("%s = %s + (%s * %s);\n", *res, x1, x2,
              horiz ? vc : vr);
        }
        free(x1);
        free(x2);
      } else if (x1) {
        if (*res) {
          printf("%s = %s;\n", *res, x1);
          free(x1);
        } else {
          *res = x1;
        }
      } else if (x2) {
        if (!*res) {
          char buffer[16];
          sprintf(buffer, "t%d", *next);
          ++*next;
          *res = strdup(buffer);
        }
        if (step > 1) {
          printf("%s = %s * %s%d;\n", *res, x2,
              horiz ? vc : vr,
              step);
        } else {
          printf("%s = %s * %s;\n", *res, x2,
              horiz ? vc : vr);
        }
        free(x2);
      } else {
        if (*res) {
          printf("%s = 0.0;\n", *res);
        } else {
          /* leave null */
        }
      }
    } else {
      if (horiz) {
        gather(m, n, i, j, 1-horiz, step>>1, next, vr, vc, res);
      } else {
        gather(m, n, i, j, 1-horiz, step, next, vr, vc, res);
      }
    }
  } else {
    if (*res) {
      printf("%s = %.8Le;\n", *res, m[i][j]);
    } else {
      char buf[32];
      if (fabs(m[i][j]) > 0.0) {
        sprintf(buf, "%.8Le", m[i][j]);
        *res = strdup(buf);
      } else {
        /* leave *res == NULL : nothing to accumulate here */
      }
    }
  }
}

/* Emit polynomial expansions using Estrin's scheme.  'n' is the max number of
 * rows/cols.  */
void estrin(Matrix m, int n, const char *vr, const char *vc, char *res)
{
  static int next = 0;

  gather(m, n, 0, 0, 0, get_step(n), &next, vr, vc, &res);

}

#if defined (TEST)

int main (int argc, char **argv)
{
  Matrix foo;
  int n;
  int i, j;

  srand48(12345);

  n = 43;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      if (drand48() < 0.2) {
        foo[i][j] = drand48() * 16.0;
        printf("foo[%d][%d] = %.10Le\n", i, j, foo[i][j]);
      } else {
        foo[i][j] = 0.0;
      }
    }
  }

  estrin(foo, n, "P", "Q", "XYZ");
  return 0;
}
#endif
