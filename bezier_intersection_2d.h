//
// Copyright © 2026 Barry Schwartz
// 
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use, copy,
// modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
// 

#ifndef BEZIER_INTERSECTION_2D_H__INCLUDED__
#define BEZIER_INTERSECTION_2D_H__INCLUDED__

//
// The algorithm here is iterated "Bézier clipping" [1] of quadratic
// Bézier curves in a plane. However, the curves will be represented
// in the s-power basis [2,3], and there will be these other
// differences from the algorithms suggested by [1]:
//
// * Before clipping, the problem will be moved to a coordinate system
//   in which the "skeleton" of the clipper curve—the line segment
//   between its endpoints—extends from (0,0) and (1,0). In this
//   coordinate system, the "distance function" is the y-coordinate.
//
// * The fat line will be the best possible parallel to the "skeleton"
//   of the clipper curve. This is easily done in the s-power basis,
//   once one has transformed coordinate systems as described above:
//   the center coefficient of the (y,t)-polynomial is exactly four
//   times the optimum signed fat line width.
//
// * Intersections with fat line boundaries will be found not by
//   extending the sides of control polygons, but instead by the
//   quadratic formula.
//
// * Mechanisms to ensure termination of the algorithm, and to improve
//   the chances of finding intersections, may differ.
//
// Inputs and outputs to the algorithms will be in single
// precision. However, the internals will be in double precision. This
// decision helps resolve some minor numerical quandaries.
//
// References:
//
// [1] Sederberg, Thomas W., "Computer Aided Geometric Design"
//     (2012). Faculty
//     Publications. 1. https://scholarsarchive.byu.edu/facpub/1
//     http://hdl.lib.byu.edu/1877/2822
//
// [2] J. Sánchez-Reyes, ‘The symmetric analogue of the polynomial
//     power basis’, ACM Transactions on Graphics, vol 16 no 3, July
//     1997, page 319.
//
// [3] J. Sánchez-Reyes, ‘Applications of the polynomial s-power basis
//     in geometry processing’, ACM Transactions on Graphics, vol 19
//     no 1, January 2000, page 35.
//

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>

#ifndef BZX2D_REAL
#define BZX2D_REAL double
#endif

typedef BZX2D_REAL bzx2d_real;

#ifndef BZX2D_ABS
#define BZX2D_ABS fabs
#endif

#ifndef BZX2D_SQRT
#define BZX2D_SQRT sqrt
#endif

#define BZX2D_CHOOSE(pred,T,F) ((pred)?(T):(F))

#ifdef NAN
#define BZX2D_ARBITRARY_REAL NAN
#else
#define BZX2D_ARBITRARY_REAL (nan (""))
#endif

#ifndef BZX2D_PRIVATE
#  ifdef __GNUC__
#    if defined __STDC_VERSION__ && 201710L < __STDC_VERSION__
#      define BZX2D_PRIVATE [[gnu::const]] static inline
#    else
#      define BZX2D_PRIVATE __attribute__((__const__)) static inline
#    endif
#  else
#    define BZX2D_PRIVATE static inline
#  endif
#endif

BZX2D_PRIVATE bool
bzx2d_between_0_and_1 (const bzx2d_real x)
{
  return (0.0 <= x) * (x <= 1.0);
}

typedef struct
{
  const bool has_smaller;
  const bool has_larger;
  const bzx2d_real smaller;
  const bzx2d_real larger;
} bzx2d_quadratic_formula_result;

BZX2D_PRIVATE bzx2d_quadratic_formula_result
bzx2d_quadratic_formula (bzx2d_real a,
                         bzx2d_real b,
                         bzx2d_real c)
{
  //
  // This routine is based on QUADPL from the Naval Surface Warfare
  // Center mathematical library nswc.f90 https://archive.ph/2IpWb
  //
  // We are interested only in real roots. Be aware that, if there are
  // two roots, they might be a double root.
  //

  if (a == 0.0)
    {
      if (b == 0.0)
        {
          // The equation is a constant equation. I do not know what
          // to do with that. There may be no solution, or there may
          // be infinitely many solutions. Let us treat it as having
          // no USEFUL solutions.
          const bzx2d_quadratic_formula_result result = {
            .has_smaller = false,
            .has_larger = false,
            .smaller = BZX2D_ARBITRARY_REAL,
            .larger = BZX2D_ARBITRARY_REAL
          };
          return result;
        }
      else
        {
          // The equation is linear. There truly is only one root.
          const bzx2d_quadratic_formula_result result = {
            .has_smaller = true,
            .has_larger = false,
            .smaller = -c / b,
            .larger = BZX2D_ARBITRARY_REAL
          };
          return result;
        }
    }
  else if (c == 0.0)
    {
      // The equation is quadratic, but the constant term vanishes.
      // One of the roots is trivially zero.
      const bzx2d_quadratic_formula_result result = {
        .has_smaller = true,
        .has_larger = true,
        .smaller = 0.0,
        .larger = -b / a
      };
      return result;
    }
  else
    {
      // Compute the discriminant, avoiding overflow.
      const BZX2D_REAL bh = b / 2.0;
      const BZX2D_REAL e =
        BZX2D_CHOOSE (BZX2D_ABS (c) <= BZX2D_ABS (bh),
                      1.0 - ((a / bh) * (c / bh)),
                      (bh * (bh / BZX2D_ABS (c)))
                      - ((1 - (2 * (c < 0.0))) * a));
      const BZX2D_REAL d =
        BZX2D_CHOOSE (BZX2D_ABS (c) <= BZX2D_ABS (bh),
                      BZX2D_SQRT (BZX2D_ABS (e)) * BZX2D_ABS (bh),
                      BZX2D_SQRT (BZX2D_ABS (e))
                      * BZX2D_SQRT (BZX2D_ABS (c)));
      if (e < 0.0)
        {
          // The solutions are complex conjugates, but we are
          // interested only in real solutions.
          const bzx2d_quadratic_formula_result result = {
            .has_smaller = false,
            .has_larger = false,
            .smaller = BZX2D_ARBITRARY_REAL,
            .larger = BZX2D_ARBITRARY_REAL
          };
          return result;
        }
      else
        {
          const BZX2D_REAL dd = (1 - (2 * (0.0 <= bh))) * d;
          const BZX2D_REAL larger = (-bh + dd) / a;
          // If larger == 0.0, the root is obviously double.
          const BZX2D_REAL smaller =
            BZX2D_CHOOSE (larger == 0.0, 0.0, (c / larger) / a);
          const bzx2d_quadratic_formula_result result = {
            .has_smaller = true,
            .has_larger = true,
            .smaller = smaller,
            .larger = larger
          };
          return result;
        }
    }
}

#endif // BEZIER_INTERSECTION_2D_H__INCLUDED__
