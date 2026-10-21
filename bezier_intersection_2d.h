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
#include <float.h>

//----------------------------------------------------------------------

#define BZX2D_CHOOSE(pred,T,F) ((pred)?(T):(F))

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

#ifdef NAN
#  define BZX2D_NAN NAN
#else
#  define BZX2D_NAN (nan (""))
#endif

//----------------------------------------------------------------------

typedef struct
{
  const double min;
  const double max;
} bzx2d_double_interval;

typedef struct
{
  const bool has_smaller;
  const bool has_larger;
  const double smaller;
  const double larger;
} bzx2d_quadratic_formula_result;

typedef struct
{
  //
  // A point is single precision because the algorithm is designed for
  // this. DO NOT NAÏVELY INCREASE PRECISION AND BELIEVE YOU ARE
  // IMPROVING THINGS.
  //
  const float x;
  const float y;
} bzx2d_point;

typedef struct
{
  // A quadratic parametric spline given by its Bernstein-basis
  // control points.
  const bzx2d_point pt0;
  const bzx2d_point pt1;
  const bzx2d_point pt2;
} bzx2d_control_point_curve2;

typedef struct
{
  // A quadratic polynomial in Sánchez-Reyes s-power basis.
  //
  // The coefficients are single precision because the algorithm is
  // designed for this. DO NOT NAÏVELY INCREASE PRECISION AND BELIEVE
  // YOU ARE IMPROVING THINGS.
  //
  const float c0;
  const float c1;
  const float c2;
} bzx2d_spower_poly2;

typedef struct
{
  // A quadratic parametric spline in Sánchez-Reyes s-power basis.
  const bzx2d_spower_poly2 x;
  const bzx2d_spower_poly2 y;
} bzx2d_spower_curve2;

typedef struct
{
  // A point relative to the origin.
  const double x;
  const double y;
} bzx2d_opoint;

typedef struct
{
  // A point relative to some other point.
  const double dx;
  const double dy;
} bzx2d_dpoint;

typedef struct
{
  //
  // A plane curve in s-power basis, represented geometrically.
  //
  // pt0 and pt2 form the curve’s ‘skeleton’. Their convex
  // combination, ((1.0-t)*pt0)+(t*pt1), form the linear backbone of
  // the curve.
  //
  // (1.0-t)*t*pt1 forms (x,t) and (y,t) parabolas that are
  // convex-up/down, and which (because this is an s-power curve)
  // vanish as the curve becomes linear.
  //
  // The actual geocurve is the locus of (1.0-t)*t*pt1 relative to the
  // skeleton, respectively for each value of t in [0,1].
  //
  // If [t0,t1] is not [0.0,1.0], then this is a portion of some
  // parent curve that is determined by context.
  //
  const bzx2d_opoint pt0;
  const bzx2d_dpoint pt1;
  const bzx2d_opoint pt2;
  const double t0;
  const double t1;
} bzx2d_geocurve2;

typedef struct
{
  //
  // An affine transformation used to project objects into the
  // coordinate system of what we are calling a "skeleton".
  //
  // If it is an opoint, first add transpose[dx,dy]. If it is a
  // dpoint, ignore dx,dy. Then, in either case, premultiply by
  // [[a,b],[c,d]]. (One can view dx,dy as a dpoint.)
  //
  const double dx, dy;
  const double a, b;
  const double c, d;
} bzx2d_affine_projection;

typedef enum {
  bzx2d_clip_everywhere,
  bzx2d_clip_nowhere,
  bzx2d_clip_left_and_right
} bzx2d_clipping_places;

typedef struct
{
  bzx2d_clipping_places where;
  double t_lft;
  double t_rgt;
} bzx2d_clipping_plan;

typedef struct
{
  const size_t ncrossings;
  const double tcrossings[4];
  const double dcrossings[4];
} bzx2d_crossings;

//----------------------------------------------------------------------

BZX2D_PRIVATE bool
bzx2d_between_0_and_1 (const double x)
{
  return (0.0 <= x) * (x <= 1.0);
}

BZX2D_PRIVATE bzx2d_quadratic_formula_result
bzx2d_quadratic_formula (const double a,
                         const double b,
                         const double c)
{
  //
  // The quadratic formula to solve 0 = ax² + bx + c.
  //
  // Naïve coding of the formula is numerically unstable. One would be
  // much better off using a numerical rootfinder. This routine,
  // however, is based on QUADPL from the Naval Surface Warfare Center
  // mathematical library nswc.f90 https://archive.ph/2IpWb
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
          const double nan_ = (double) BZX2D_NAN;
          const bzx2d_quadratic_formula_result result = {
            .has_smaller = false,
            .has_larger = false,
            .smaller = nan_,
            .larger = nan_
          };
          return result;
        }
      else
        {
          // The equation is linear. There truly is only one root.
          const double nan_ = (double) BZX2D_NAN;
          const bzx2d_quadratic_formula_result result = {
            .has_smaller = true,
            .has_larger = false,
            .smaller = -c / b,
            .larger = nan_
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
      const double bh = b / 2.0;
      const double e =
        BZX2D_CHOOSE (fabs (c) <= fabs (bh),
                      1.0 - ((a / bh) * (c / bh)),
                      (bh * (bh / fabs (c)))
                      - ((1 - (2 * (c < 0.0))) * a));
      const double d =
        BZX2D_CHOOSE (fabs (c) <= fabs (bh),
                      sqrt (fabs (e)) * fabs (bh),
                      sqrt (fabs (e))
                      * sqrt (fabs (c)));
      if (e < 0.0)
        {
          // The solutions are complex conjugates, but we are
          // interested only in real solutions.
          const double nan_ = (double) BZX2D_NAN;
          const bzx2d_quadratic_formula_result result = {
            .has_smaller = false,
            .has_larger = false,
            .smaller = nan_,
            .larger = nan_
          };
          return result;
        }
      else
        {
          const double larger =
            (-bh + (1 - (2 * (0.0 <= bh))) * d) / a;
          // If larger == 0.0, the root is obviously double.
          const double smaller =
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

BZX2D_PRIVATE double
bzx2d_eval_spower_poly2 (const bzx2d_spower_poly2 poly,
                         const double t)
{
  return
    (BZX2D_CHOOSE
     (t <= 0.5,
      poly.c0 + (t * ((poly.c2 - poly.c0) + ((1.0 - t) * poly.c1))),
      poly.c2 + ((1.0 - t) * ((poly.c0 - poly.c2) + (t * poly.c1)))));
}

BZX2D_PRIVATE bzx2d_point
bzx2d_eval_spower_curve2 (const bzx2d_spower_curve2 curve,
                          const double t)
{
  const bzx2d_point result = {
    .x = bzx2d_eval_spower_poly2 (curve.x, t),
    .y = bzx2d_eval_spower_poly2 (curve.y, t)
  };
  return result;
}

BZX2D_PRIVATE bzx2d_spower_curve2
bzx2d_to_spower_curve (const bzx2d_control_point_curve2 ctrl)
{
  const double pt0_x = (double) ctrl.pt0.x;
  const double pt1_x = (double) ctrl.pt1.x;
  const double pt2_x = (double) ctrl.pt2.x;
  const double pt0_y = (double) ctrl.pt0.y;
  const double pt1_y = (double) ctrl.pt1.y;
  const double pt2_y = (double) ctrl.pt2.y;
  const bzx2d_spower_curve2 result = {
    .x = {
      .c0 = ctrl.pt0.x,
      .c1 = (float) ((2.0 * pt1_x) - pt0_x - pt2_x),
      .c2 = ctrl.pt2.x
    },
    .y = {
      .c0 = ctrl.pt0.y,
      .c1 = (float) ((2.0 * pt1_y) - pt0_y - pt2_y),
      .c2 = ctrl.pt2.y
    }
  };
  return result;
}

BZX2D_PRIVATE bzx2d_dpoint
bzx2d_sub_opoint_opoint (const bzx2d_opoint pt1,
                         const bzx2d_opoint pt2)
{
  const bzx2d_dpoint result = {
    .dx = pt1.x - pt2.x,
    .dy = pt1.y - pt2.y
  };
  return result;
}

BZX2D_PRIVATE bzx2d_opoint
bzx2d_add_opoint_dpoint (const bzx2d_opoint opt,
                         const bzx2d_dpoint dpt)
{
  const bzx2d_opoint result = {
    .x = opt.x + dpt.dx,
    .y = opt.y + dpt.dy
  };
  return result;
}

BZX2D_PRIVATE bzx2d_dpoint
bzx2d_add_dpoint_dpoint (const bzx2d_dpoint pt1,
                         const bzx2d_dpoint pt2)
{
  const bzx2d_dpoint result = {
    .dx = pt1.dx + pt2.dx,
    .dy = pt1.dy + pt2.dy
  };
  return result;
}

BZX2D_PRIVATE bzx2d_dpoint
bzx2d_mul_scalar_dpoint (const double scalar,
                         const bzx2d_dpoint pt)
{
  const bzx2d_dpoint result = {
    .dx = scalar * pt.dx,
    .dy = scalar * pt.dy
  };
  return result;
}

BZX2D_PRIVATE bzx2d_geocurve2
bzx2d_to_geocurve (const bzx2d_spower_curve2 scurve)
{
  const bzx2d_geocurve2 result = {
    .pt0 = {
      .x = (double) scurve.x.c0,
      .y = (double) scurve.y.c0
    },
    .pt1 = {
      .dx = (double) scurve.x.c1,
      .dy = (double) scurve.y.c1
    },
    .pt2 = {
      .x = (double) scurve.x.c2,
      .y = (double) scurve.y.c2
    },
    .t0 = 0.0,
    .t1 = 1.0
  };
  return result;
}

BZX2D_PRIVATE bzx2d_opoint
bzx2d_eval_geocurve2 (const bzx2d_geocurve2 gcurve,
                      const double t)
{
  //
  // if t <= 0.5
  // then pt0 + t * ((pt2 - pt0) + ((1.0 - t) * pt1))
  // else pt2 + (1.0 - t) * ((pt0 - pt2) + (t * pt1))
  //
  return
    (BZX2D_CHOOSE
     (t <= 0.5,
      bzx2d_add_opoint_dpoint
      (gcurve.pt0,
       bzx2d_mul_scalar_dpoint
       (t, (bzx2d_add_dpoint_dpoint
            (bzx2d_sub_opoint_opoint (gcurve.pt2, gcurve.pt0),
             bzx2d_mul_scalar_dpoint (1.0 - t, gcurve.pt1))))),
      bzx2d_add_opoint_dpoint
      (gcurve.pt2,
       bzx2d_mul_scalar_dpoint
       (1.0 - t,
        (bzx2d_add_dpoint_dpoint
         (bzx2d_sub_opoint_opoint (gcurve.pt0, gcurve.pt2),
          bzx2d_mul_scalar_dpoint (t, gcurve.pt1)))))));
}

BZX2D_PRIVATE bzx2d_opoint
bzx2d_project_opoint (const bzx2d_opoint pt,
                      const bzx2d_affine_projection trans)
{
  //
  // Apply an affine transformation to an opoint, projecting it from
  // one coordinate system to another.
  //

  // Change of origin.
  const double x = pt.x + trans.dx;
  const double y = pt.y + trans.dy;

  // Rotation, scaling, etc.
  const bzx2d_opoint result = {
    .x = (trans.a * x) + (trans.b * y),
    .y = (trans.c * x) + (trans.d * y)
  };
  return result;
}

BZX2D_PRIVATE bzx2d_dpoint
bzx2d_project_dpoint (const bzx2d_dpoint pt,
                      const bzx2d_affine_projection trans)
{
  //
  // Apply an affine transformation to a dpoint, projecting it from
  // one coordinate system to another.
  //

  // For a dpoint, there is no change of origin, because a dpoint is
  // independent of the origin.
  const double dx = trans.dx;
  const double dy = trans.dy;

  // Rotation, scaling, etc.
  const bzx2d_dpoint result = {
    .dx = (trans.a * dx) + (trans.b * dy),
    .dy = (trans.c * dx) + (trans.d * dy)
  };
  return result;
}

BZX2D_PRIVATE bzx2d_geocurve2
bzx2d_project_geocurve2 (const bzx2d_geocurve2 gcurve,
                         const bzx2d_affine_projection trans)
{
  //
  // Apply an affine transformation to a geocurve, projecting it from
  // one coordinate system to another.
  //
  const bzx2d_geocurve2 result = {
    .pt0 = bzx2d_project_opoint (gcurve.pt0, trans),
    .pt1 = bzx2d_project_dpoint (gcurve.pt1, trans),
    .pt2 = bzx2d_project_opoint (gcurve.pt2, trans),
    .t0 = 0.0,
    .t1 = 1.0
  };
  return result;
}

BZX2D_PRIVATE bzx2d_affine_projection
bzx2d_segment_projection_to_0_1 (const bzx2d_opoint pt0,
                                 const bzx2d_opoint pt1)
{
  //
  // Return the transformation that projects pt0 to (0,0) and pt1 to
  // (1,0) without any distortions except translation, rotation, and
  // scaling.
  //
  const double a1 = pt1.x - pt0.x;
  const double b1 = pt1.y - pt0.y;
  const double denom = (a1 * a1) + (b1 * b1);
  const double a = a1 / denom;
  const double b = b1 / denom;
  const bzx2d_affine_projection result = {
    .dx = -pt0.x,
    .dy = -pt0.y,
    .a = a,
    .b = b,
    .c = -b,
    .d = a
  };
  return result;
}

BZX2D_PRIVATE bzx2d_affine_projection
bzx2d_skeleton_projection (const bzx2d_geocurve2 clipper)
{
  //
  // Return the transformation that projects the “skeleton” of the
  // clipper into a coordinate system where the skeleton goes from the
  // origin (0,0) to the point (1,0).
  //
  return bzx2d_segment_projection_to_0_1 (clipper.pt0, clipper.pt2);
}

BZX2D_PRIVATE bzx2d_double_interval
bzx2d_fat_line_bounds (const bzx2d_geocurve2 projected_clipper)
{
  //
  // (1-t)*t is bounded on [0,1] by 0 and 1/4. Thus we have our fat
  // line boundaries. We leave some padding, however, for numerical
  // safety. The padding is much smaller than our tolerance. (We are
  // calculating in double precision, but to obtain results in single
  // precision.)
  //
  // The fat line bounds will be y=result.min and y=result.max. If you
  // ignore the padding, the fat line is the best possible among those
  // parallel to the skeleton. It is likely to be much tighter than
  // what you would infer from the Bernstein-basis control polygon.
  //
  const double padding = ((double) FLT_EPSILON) / 1000.0;
  const double d = projected_clipper.pt1.dy / 4.0;
  const double d1 = -padding * d;
  const double d2 = (1.0 + padding) * d;
  const bzx2d_double_interval result = {
    .min = fmin (d1, d2),
    .max = fmax (d1, d2)
  };
  return result;
}

BZX2D_PRIVATE bool
bzx2d_outside_of_bounds (const bzx2d_geocurve2 projected_clippehend,
                         const bzx2d_double_interval d)
{
  //
  // (1-t)*t is bounded on [0,1] by 0 and 1/4. Use that fact to make
  // padding for a simple (although not especially tight) bounds on y.
  //
  const double dy = fabs (projected_clippehend.pt1.dy) / 4.0;
  const double ymin =
    -dy + fmin (projected_clippehend.pt0.y, projected_clippehend.pt2.y);
  const double ymax =
    dy + fmax (projected_clippehend.pt0.y, projected_clippehend.pt2.y);
  return (ymax <= d.min || d.max <= ymin);
}

BZX2D_PRIVATE bzx2d_crossings
bzx2d_make_crossings (void)
{
  const double nan_ = (double) BZX2D_NAN;
  const bzx2d_crossings result = {
    .ncrossings = 0,
    .tcrossings = {
      [0] = nan_,
      [1] = nan_,
      [2] = nan_,
      [3] = nan_
    },
    .dcrossings = {
      [0] = nan_,
      [1] = nan_,
      [2] = nan_,
      [3] = nan_
    }
  };
  return result;
}

BZX2D_PRIVATE bzx2d_crossings
bzx2d_add_crossing (const double t,
                    const double yvalue,
                    const bzx2d_crossings crossings0)
{
  switch (crossings0.ncrossings)
    {
    case 0:
      {
        const bzx2d_crossings crossings1 = {
          .ncrossings = 1,
          .tcrossings = {
            [0] = t,
            [1] = crossings0.tcrossings[1],
            [2] = crossings0.tcrossings[2],
            [3] = crossings0.tcrossings[3]
          },
          .dcrossings = {
            [0] = yvalue,
            [1] = crossings0.dcrossings[1],
            [2] = crossings0.dcrossings[2],
            [3] = crossings0.dcrossings[3]
          }
        };
        return crossings1;
      }
    case 1:
      {
        const bzx2d_crossings crossings1 = {
          .ncrossings = 2,
          .tcrossings = {
            [0] = crossings0.tcrossings[0],
            [1] = t,
            [2] = crossings0.tcrossings[2],
            [3] = crossings0.tcrossings[3]
          },
          .dcrossings = {
            [0] = crossings0.dcrossings[0],
            [1] = yvalue,
            [2] = crossings0.dcrossings[2],
            [3] = crossings0.dcrossings[3]
          }
        };
        return crossings1;
      }
    case 2:
      {
        const bzx2d_crossings crossings1 = {
          .ncrossings = 3,
          .tcrossings = {
            [0] = crossings0.tcrossings[0],
            [1] = crossings0.tcrossings[1],
            [2] = t,
            [3] = crossings0.tcrossings[3]
          },
          .dcrossings = {
            [0] = crossings0.dcrossings[0],
            [1] = crossings0.dcrossings[1],
            [2] = yvalue,
            [3] = crossings0.dcrossings[3]
          }
        };
        return crossings1;
      }
    case 3:
      {
        const bzx2d_crossings crossings1 = {
          .ncrossings = 4,
          .tcrossings = {
            [0] = crossings0.tcrossings[0],
            [1] = crossings0.tcrossings[1],
            [2] = crossings0.tcrossings[2],
            [3] = t
          },
          .dcrossings = {
            [0] = crossings0.dcrossings[0],
            [1] = crossings0.dcrossings[1],
            [2] = crossings0.dcrossings[2],
            [3] = yvalue
          }
        };
        return crossings1;
      }
    default:
      // Programmer error.
      abort ();
    }
}

BZX2D_PRIVATE bzx2d_crossings
bzx2d_find_crossings (const double c0,
                      const double c1,
                      const double c2,
                      const double yvalue,
                      const bzx2d_crossings crossings0)
{
  const bzx2d_quadratic_formula_result roots =
    bzx2d_quadratic_formula (-c1, c1 + c2 - c0, c0 - yvalue);
  const bzx2d_crossings crossings1 =
    BZX2D_CHOOSE ((roots.has_smaller *
                   bzx2d_between_0_and_1 (roots.smaller)),
                  bzx2d_add_crossing (roots.smaller, yvalue,
                                      crossings0),
                  crossings0);
  const bzx2d_crossings crossings2 =    
    BZX2D_CHOOSE ((roots.has_larger *
                   bzx2d_between_0_and_1 (roots.larger)),
                  bzx2d_add_crossing (roots.larger, yvalue,
                                      crossings1),
                  crossings1);
  return crossings2;
}

BZX2D_PRIVATE bzx2d_clipping_plan
bzx2d_create_a_clipping_plan (const bzx2d_crossings crossings)
{
  //
  // For simplicity, there is just one situation in which we will
  // clip: when there is one crossing of d.min and one crossing of
  // d.max.
  //
  if (crossings.ncrossings == 2
      && crossings.dcrossings[0] != crossings.dcrossings[1])
    {
      const bzx2d_clipping_plan result = {
        .where = bzx2d_clip_left_and_right,
        .t_lft = fmin (crossings.tcrossings[0],
                       crossings.tcrossings[1]),
        .t_rgt = fmax (crossings.tcrossings[0],
                       crossings.tcrossings[1])
      };
      return result;
    }
  else
    {
      const double nan_ = (double) BZX2D_NAN;
      const bzx2d_clipping_plan result = {
        .where = bzx2d_clip_nowhere,
        .t_lft = nan_,
        .t_rgt = nan_
      };
      return result;
    }
}

BZX2D_PRIVATE bzx2d_clipping_plan
bzx2d_analyze_and_plan (const bzx2d_geocurve2 projected_clippehend,
                        const bzx2d_double_interval d)
{
  const double c0 = projected_clippehend.pt0.y;
  const double c1 = projected_clippehend.pt1.dy;
  const double c2 = projected_clippehend.pt2.y;
  const bzx2d_crossings crossings0 = bzx2d_make_crossings ();
  const bzx2d_crossings crossings1 =
    bzx2d_find_crossings (c0, c1, c2, d.min, crossings0);
  const bzx2d_crossings crossings2 =
    bzx2d_find_crossings (c0, c1, c2, d.max, crossings1);
  return bzx2d_create_a_clipping_plan (crossings2);
}

BZX2D_PRIVATE bzx2d_clipping_plan
bzx2d_make_clipping_plan (const bzx2d_geocurve2 projected_clippehend,
                          const bzx2d_double_interval d)
{
  if (bzx2d_outside_of_bounds (projected_clippehend, d))
    {
      const double nan_ = (double) BZX2D_NAN;
      const bzx2d_clipping_plan result = {
        .where = bzx2d_clip_everywhere,
        .t_lft = nan_,
        .t_rgt = nan_
      };
      return result;
    }
  else
    return bzx2d_analyze_and_plan (projected_clippehend, d);
}

#endif // BEZIER_INTERSECTION_2D_H__INCLUDED__
