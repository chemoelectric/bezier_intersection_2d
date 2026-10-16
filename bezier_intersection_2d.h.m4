//----------------------------------------------------------------------
//
// Intersection of quadratic and cubic parametric splines in the
// plane, using the symmetric power (s-power) basis and fat-line
// Bézier clipping.
//
// Written in C99 as pure static inline routines in tail-recursive
// style.
//
// References
// ----------
//
//   Javier Sánchez-Reyes. 1997. The symmetric analogue of the
//   polynomial power basis. ACM Trans. Graph. 16, 3 (July 1997),
//   319–357. DOI: https://doi.org/10.1145/256157.256162
//
//   Javier Sánchez-Reyes. 2000. Applications of the polynomial
//   s-power basis in geometry processing. ACM Trans. Graph. 19, 1
//   (Jan. 2000), 27–55. DOI: https://doi.org/10.1145/343002.343018
//
//   Thomas W. Sederberg and Tomoyuki Nishita. 1990. Curve
//   intersection using Bézier clipping. Comput. Aided Des. 22, 9
//   (Nov. 1990), 538–549.
//   https://doi.org/10.1016/0010-4485(90)90039-F
//
//----------------------------------------------------------------------
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
//----------------------------------------------------------------------

m4_divert(-1)

m4_changequote(`[:',`:]')
m4_changecom

m4_define([:m4_comma_newline_if_under_eight:],
  [:m4_ifelse(m4_eval([:($1) <= 7:]),[:1:],
[:,
:],[::]):])

m4_dnl
m4_dnl  m4_forloop(variable, from, to, statement)
m4_dnl
m4_define([:m4_forloop:],[:m4_ifelse(m4_eval([:($2) <= ($3):]),[:1:],
  [:m4_pushdef([:$1:])_$0([:$1:],m4_eval([:$2:]),
    m4_eval([:$3:]),[:$4:])m4_popdef([:$1:]):]):])
m4_define([:_m4_forloop:],
  [:m4_define([:$1:],[:$2:])$4[::]m4_ifelse([:$2:],[:$3:],[::],
    [:$0([:$1:],m4_incr([:$2:]),[:$3:],[:$4:]):]):])

m4_divert[::]m4_dnl

#ifndef BEZIER_INTERSECTION_2D_H__INCLUDED__
#define BEZIER_INTERSECTION_2D_H__INCLUDED__

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>

#ifndef INTERSECTION_REAL
#define INTERSECTION_REAL double
#endif

typedef INTERSECTION_REAL intersection_real;
#define INTERSECTION_REAL_ZERO ((intersection_real) 0.0)
#define INTERSECTION_REAL_HALF ((intersection_real) 0.5)
#define INTERSECTION_REAL_0NE ((intersection_real) 1.0)
#define INTERSECTION_REAL_TWO ((intersection_real) 2.0)
#define INTERSECTION_REAL_THREE ((intersection_real) 3.0)
#define INTERSECTION_REAL_FOUR ((intersection_real) 4.0)

#ifndef INTERSECTION_DEFAULT_RELTOL
#define INTERSECTION_DEFAULT_RELTOL 1e-7
#endif

#ifndef INTERSECTION_DEFAULT_VERYSMALL
#define INTERSECTION_DEFAULT_VERYSMALL 1e-14
#endif

#ifndef INTERSECTION_DEFAULT_PADDING
#define INTERSECTION_DEFAULT_PADDING 1e-12
#endif

#ifndef INTERSECTION_REAL_SQRT
#define INTERSECTION_REAL_SQRT sqrt
#endif

#ifndef INTERSECTION_MAXIMUM_DEPTH
#define INTERSECTION_MAXIMUM_DEPTH 60
#endif

#ifndef BEZIER_INTERSECTION_2D_CONST
#  if defined __GNUC__
#    if defined __STDC_VERSION__ && 201710L < __STDC_VERSION__
#      define BEZIER_INTERSECTION_2D_CONST [[gnu::const]]
#    else
#      define BEZIER_INTERSECTION_2D_CONST __attribute__((__const__))
#    endif
#  else
#    define BEZIER_INTERSECTION_2D_CONST
#  endif
#endif

#ifndef BEZIER_INTERSECTION_2D_INLINE
#define BEZIER_INTERSECTION_2D_INLINE inline
#endif

#ifndef BEZIER_INTERSECTION_2D_STATIC
#define BEZIER_INTERSECTION_2D_STATIC static
#endif

#ifndef BEZIER_INTERSECTION_2D_STATICINLINE
#define BEZIER_INTERSECTION_2D_STATICINLINE \
  BEZIER_INTERSECTION_2D_CONST \
  BEZIER_INTERSECTION_2D_INLINE \
  BEZIER_INTERSECTION_2D_STATIC
#endif

typedef struct
{
  bool degenerate;
  size_t n;
  intersection_real tp[9];  // p parameters
  intersection_real tq[9];  // q parameters
} intersection_record;

// A point in the plane.
typedef struct
{
  intersection_real x;
  intersection_real y;
} intersection_real_point;

//
// Quadratic curve in s-power basis:
//
//     C(t) = (1 - t)*P0 + t*P2 + (1 - t)*t*P1_delta
//
//       for t_min <= t <= t_max
//
typedef struct
{
  intersection_real_point p0;
  intersection_real_point p2;
  intersection_real_point p1_delta;
  intersection_real t_min;
  intersection_real t_max;
} intersection_real_quadratic;

BEZIER_INTERSECTION_2D_STATICINLINE intersection_real
intersection_real_nan (void)
{
#if defined(NAN)
  return (intersection_real) NAN;
#else
  return (intersection_real) nan ("");
#endif
}

BEZIER_INTERSECTION_2D_STATICINLINE intersection_real
intersection_real_abs (const intersection_real x)
{
  return (x < INTERSECTION_REAL_ZERO) ? (-x) : x;
}

BEZIER_INTERSECTION_2D_STATICINLINE intersection_real
intersection_real_sgn (const intersection_real x)
{
  return ((x < INTERSECTION_REAL_ZERO)
          ? (intersection_real) (-1.0)
          : ((x == INTERSECTION_REAL_ZERO)
             ? INTERSECTION_REAL_ZERO
             : INTERSECTION_REAL_0NE));
}

BEZIER_INTERSECTION_2D_STATICINLINE intersection_real
intersection_real_square (const intersection_real x)
{
  return (x * x);
}

BEZIER_INTERSECTION_2D_STATICINLINE intersection_real
intersection_real_sqrt (const intersection_real x)
{
  return INTERSECTION_REAL_SQRT (x);
}

BEZIER_INTERSECTION_2D_STATICINLINE intersection_real
intersection_real_min (const intersection_real x,
                       const intersection_real y)
{
  return (x < y) ? x : y;
}

BEZIER_INTERSECTION_2D_STATICINLINE intersection_real
intersection_real_max (const intersection_real x,
                       const intersection_real y)
{
  return (x < y) ? y : x;
}

BEZIER_INTERSECTION_2D_STATICINLINE intersection_real
intersection_real_bound_0_1 (const intersection_real x)
{
  return ((x < INTERSECTION_REAL_ZERO)
          ? (INTERSECTION_REAL_ZERO)
          : ((INTERSECTION_REAL_0NE < x)
             ? (INTERSECTION_REAL_0NE)
             : x));
}

BEZIER_INTERSECTION_2D_STATICINLINE intersection_record
make_intersection_record (const bool degenerate,
                          const size_t n,
                          //
m4_forloop([:i:],0,8,
[:                          const intersection_real tp[::]i[::],
:])m4_dnl
                          //
m4_forloop([:i:],0,8,
[:                          const intersection_real tq[::]i[::]m4_dnl
m4_comma_newline_if_under_eight(i):]))
{
  const intersection_record result = {
    .degenerate = degenerate,
    .n = n,

m4_forloop([:i:],0,8,
[:    .tp[i] = tp[::]i[::],
:])m4_dnl

m4_forloop([:i:],0,8,
[:    .tq[i] = tq[::]i[::]m4_dnl
m4_comma_newline_if_under_eight(i):])
  };
  return result;
}

BEZIER_INTERSECTION_2D_STATICINLINE intersection_record
empty_intersection_record (void)
{
  const intersection_real nan = intersection_real_nan ();
  return make_intersection_record (false, 0,
                                   m4_forloop([:i:],0,16,[:nan, :])nan);
}

BEZIER_INTERSECTION_2D_STATICINLINE bool
intersection_record_contains (const intersection_record rec,
                              const intersection_real tp,
                              const intersection_real tq,
                              const intersection_real tol,
                              const size_t i)
{
  if (i == rec.n)
    return false;
  else
    {
      const intersection_real dp =
        (tp < rec.tp[i]) ? (rec.tp[i] - tp) : (tp - rec.tp[i]);
      const intersection_real dq =
        (tq < rec.tq[i]) ? (rec.tq[i] - tq) : (tq - rec.tq[i]);
      if (dp <= tol && dq <= tol)
        return true;
      else
        return intersection_record_contains (rec, tp, tq, tol, i + 1);
    }
}

BEZIER_INTERSECTION_2D_STATICINLINE intersection_record
intersection_record_add (const intersection_record rec,
                         const intersection_real tp,
                         const intersection_real tq,
                         const intersection_real tol)
{
  const intersection_real tp_bounded = intersection_real_bound_0_1 (tp);
  const intersection_real tq_bounded = intersection_real_bound_0_1 (tq);

  if (intersection_record_contains (rec, tp_bounded, tq_bounded,
                                    tol, 0))
    return rec;
  else
    switch (rec.n)
      {
m4_forloop([:i:],0,8,
[:      case i:
          return (make_intersection_record
                  (rec.degenerate, rec.n + 1,
                   m4_ifelse(i,0,tp_bounded,rec.tp[0]),
                   m4_ifelse(i,1,tp_bounded,rec.tp[1]),
                   m4_ifelse(i,2,tp_bounded,rec.tp[2]),
                   m4_ifelse(i,3,tp_bounded,rec.tp[3]),
                   m4_ifelse(i,4,tp_bounded,rec.tp[4]),
                   m4_ifelse(i,5,tp_bounded,rec.tp[5]),
                   m4_ifelse(i,6,tp_bounded,rec.tp[6]),
                   m4_ifelse(i,7,tp_bounded,rec.tp[7]),
                   m4_ifelse(i,8,tp_bounded,rec.tp[8]),
                   m4_ifelse(i,0,tq_bounded,rec.tq[0]),
                   m4_ifelse(i,1,tq_bounded,rec.tq[1]),
                   m4_ifelse(i,2,tq_bounded,rec.tq[2]),
                   m4_ifelse(i,3,tq_bounded,rec.tq[3]),
                   m4_ifelse(i,4,tq_bounded,rec.tq[4]),
                   m4_ifelse(i,5,tq_bounded,rec.tq[5]),
                   m4_ifelse(i,6,tq_bounded,rec.tq[6]),
                   m4_ifelse(i,7,tq_bounded,rec.tq[7]),
                   m4_ifelse(i,8,tq_bounded,rec.tq[8])));
:])m4_dnl
        default:
          return rec;
    }
}

BEZIER_INTERSECTION_2D_STATICINLINE intersection_real
intersection_real_distance_sq (const intersection_real_point a,
                               const intersection_real_point b)
{
  return (intersection_real_square (b.x - a.x)
          + intersection_real_square (b.y - a.y));
}

BEZIER_INTERSECTION_2D_STATICINLINE intersection_real_point
intersection_real_eval_quadratic (const intersection_real_quadratic q,
                                  const intersection_real u)
{
  const intersection_real s = INTERSECTION_REAL_0NE - u;
  const intersection_real blend = s * u;
  const intersection_real_point pt = {
    .x = (s * q.p0.x) + (u * q.p2.x) + (blend * q.p1_delta.x),
    .y = (s * q.p0.y) + (u * q.p2.y) + (blend * q.p1_delta.y)
  };
  return pt;
}

BEZIER_INTERSECTION_2D_STATICINLINE intersection_real_quadratic
intersection_real_subdivide_quadratic (const intersection_real_quadratic q,
                                       const intersection_real u0,
                                       const intersection_real u1)
{
  const intersection_real du = u1 - u0;
  const intersection_real orig_span = q.t_max - q.t_min;

  const intersection_real_point p0 =
    intersection_real_eval_quadratic (q, u0);
  const intersection_real_point p2 =
    intersection_real_eval_quadratic (q, u1);

  //
  // Compute the s-power deviation vector:
  //
  //   P1_delta = 4 * ( C((u0+u1)/2) - 0.5*(P0 + P2) )
  //
  const intersection_real_point mid =
    intersection_real_eval_quadratic
    (q, INTERSECTION_REAL_HALF * (u0 + u1));
  const intersection_real_point lin_mid = {
    .x = INTERSECTION_REAL_HALF * (p0.x + p2.x),
    .y = INTERSECTION_REAL_HALF * (p0.y + p2.y)
  };
  const intersection_real_quadratic result = {
    .t_min = q.t_min + u0 * orig_span,
    .t_max = q.t_min + u1 * orig_span,
    .p0 = p0,
    .p2 = p2,
    .p1_delta = {
      .x = INTERSECTION_REAL_FOUR * (mid.x - lin_mid.x),
      .y = INTERSECTION_REAL_FOUR * (mid.y - lin_mid.y)
    }
  };

  return result;
}

typedef struct
{
  const size_t nroots;
  const intersection_real root0;
  const intersection_real root1;
  const intersection_real root2;
  const intersection_real root3;
} intersection_real_quadratic_roots;

BEZIER_INTERSECTION_2D_STATICINLINE const intersection_real_quadratic_roots
make_intersection_real_quadratic_roots (void)
{
  const intersection_real nan = intersection_real_nan ();
  const intersection_real_quadratic_roots result = {
    .nroots = 0,
    .root0 = nan,
    .root1 = nan,
    .root2 = nan,
    .root3 = nan
  };
  return result;
}

BEZIER_INTERSECTION_2D_STATICINLINE const intersection_real_quadratic_roots
intersection_real_quadratic_insert_root (const intersection_real_quadratic_roots qroots,
                                         const intersection_real new_root)
{
  switch (qroots.nroots)
    {
    case 0:
      {
        const intersection_real_quadratic_roots result = {
          .nroots = 1,
          .root0 = new_root,
          .root1 = qroots.root1,
          .root2 = qroots.root2,
          .root3 = qroots.root3
        };
        return result;
      }
    case 1:
      {
        const intersection_real_quadratic_roots result = {
          .nroots = 2,
          .root0 = qroots.root0,
          .root1 = new_root,
          .root2 = qroots.root2,
          .root3 = qroots.root3
        };
        return result;
      }
    case 2:
      {
        const intersection_real_quadratic_roots result = {
          .nroots = 3,
          .root0 = qroots.root0,
          .root1 = qroots.root1,
          .root2 = new_root,
          .root3 = qroots.root3
        };
        return result;
      }
    case 3:
      {
        const intersection_real_quadratic_roots result = {
          .nroots = 4,
          .root0 = qroots.root0,
          .root1 = qroots.root1,
          .root2 = qroots.root2,
          .root3 = new_root
        };
        return result;
      }
    default:
      // Programmer error.
      abort ();
    }
}

BEZIER_INTERSECTION_2D_STATIC const intersection_real_quadratic_roots
intersection_real_quadratic_roots_solve
(const intersection_real_quadratic_roots qroots,
 const intersection_real a,
 const intersection_real b,
 const intersection_real c)
{
  //
  // Solve av² + bv + c = 0.
  //

  if (intersection_real_abs (a) <
      ((intersection_real) INTERSECTION_DEFAULT_VERYSMALL))
    {
      if (((intersection_real) INTERSECTION_DEFAULT_VERYSMALL) <
          intersection_real_abs (b))
        return (intersection_real_quadratic_insert_root
                (qroots, -c / b));
      else
        return qroots;
    }
  else
    {
      const intersection_real discrim =
        (b * b) - (INTERSECTION_REAL_FOUR * a * c);
      if (INTERSECTION_REAL_ZERO <= discrim)
        {
          const intersection_real s_discrim =
            intersection_real_sqrt (discrim);
          return (intersection_real_quadratic_insert_root
                  ((intersection_real_quadratic_insert_root
                    (qroots, (-b - s_discrim) / (a + a))),
                   ((-b + s_discrim) / (a + a))));
        }
      else
        return qroots;
    }
}

BEZIER_INTERSECTION_2D_STATICINLINE intersection_record
intersection_of_two_quadratics (const intersection_real_quadratic p,
                                const intersection_real_quadratic q,
                                const intersection_real rel_tol,
                                const intersection_real abs_tol,
                                const int depth,
                                const intersection_record accum0)
{
  switch (accum0.n)
    {
m4_forloop([:i:],0,8,[:case i:
                      :])
 const intersection_real p_span = p.t_max - p.t_min;
 const intersection_real q_span = q.t_max - q.t_min;

 if ((p_span <= rel_tol) * (q_span <= rel_tol))
   // Both parameter spans are within relative tolerance. Record a
   // root.
   return (intersection_record_add
           (accum0, INTERSECTION_REAL_HALF * (p.t_min + p.t_max),
            INTERSECTION_REAL_HALF * (q.t_min + q.t_max),
            rel_tol + rel_tol));
 else if (INTERSECTION_MAXIMUM_DEPTH < depth)
   // Limit recursion depth to prevent infinite loops.
   return (intersection_record_add
           (accum0, INTERSECTION_REAL_HALF * (p.t_min + p.t_max),
            INTERSECTION_REAL_HALF * (q.t_min + q.t_max),
            rel_tol + rel_tol));
 else
   {
     // Skeleton baseline vector of curve P.
     const intersection_real skel_dx = p.p2.x - p.p0.x;
     const intersection_real skel_dy = p.p2.y - p.p0.y;
     const intersection_real skel_len_sq =
       (skel_dx * skel_dx) + (skel_dy * skel_dy);

     // If curve P's skeleton has collapsed to a point, swap roles or
     // split.
     if (skel_len_sq <= abs_tol * abs_tol)
       {
         const intersection_real q_skel_dx = q.p2.x - q.p0.x;
         const intersection_real q_skel_dy = q.p2.y - q.p0.y;
         const intersection_real q_skel_len_sq =
           (q_skel_dx * q_skel_dx) + (q_skel_dy * q_skel_dy);

         if (q_skel_len_sq <= abs_tol * abs_tol)
           {
             // Both curves have collapsed to points.
             if (intersection_real_distance_sq (p.p0, q.p0)
                 <= abs_tol * abs_tol)
               return (intersection_record_add
                       (accum0, p.t_min, q.t_min, rel_tol + rel_tol));
             else
               return accum0;
           }
         else
           // Swap P and Q.
           return intersection_of_two_quadratics (q, p, rel_tol, abs_tol,
                                                  depth + 1, accum0);
       }
     else
       {
         const intersection_real skel_len =
           intersection_real_sqrt (skel_len_sq);

         // Transverse deviation of P from its skeleton line.
         const intersection_real d_y =
           ((-skel_dy * p.p1_delta.x) + (skel_dx * p.p1_delta.y))
           / skel_len;

         //
         // Extreme bounds of P along perpendicular axis. The
         // polynomial is
         //
         //   (1-t)*0 + t*0 + (1-t)*t*d_y = (1-t)*t*d_y
         //
         const intersection_real d_min =
           (INTERSECTION_REAL_ZERO <= d_y)
           ? -(intersection_real) INTERSECTION_DEFAULT_PADDING
           : ((intersection_real) 0.25 * d_y
              - (intersection_real) INTERSECTION_DEFAULT_PADDING);
         const intersection_real d_max =
           (INTERSECTION_REAL_ZERO <= d_y)
           ? ((intersection_real) 0.25 * d_y
              + (intersection_real) INTERSECTION_DEFAULT_PADDING)
           : (intersection_real) INTERSECTION_DEFAULT_PADDING;

         // Project Q onto perpendicular distance to P's skeleton
         // line.
         const intersection_real q0_perp =
           (-skel_dy * (q.p0.x - p.p0.x) + skel_dx * (q.p0.y - p.p0.y))
           / skel_len;
         const intersection_real q2_perp =
           (-skel_dy * (q.p2.x - p.p0.x) + skel_dx * (q.p2.y - p.p0.y))
           / skel_len;
         const intersection_real qd_perp =
           (-skel_dy * q.p1_delta.x + skel_dx * q.p1_delta.y)
           / skel_len;

         //
         // Q's perpendicular coordinate polynomial:
         //
         //   y(v) = (1 - v)*q0_perp + v*q2_perp + (1 - v)*v*qd_perp
         //        = a*v^2 + b*v + c
         //
         const intersection_real a = -qd_perp;
         const intersection_real b = q2_perp - q0_perp + qd_perp;
         const intersection_real c = q0_perp;

         // Compute range of y(v) for 0 <= v <= 1
         const bool a_is_very_small =
           (intersection_real_abs (a)
            <= (intersection_real) INTERSECTION_DEFAULT_VERYSMALL);
         const intersection_real v_vertex =
           ((a_is_very_small)
            ? (intersection_real_sgn (b) * intersection_real_sgn (a)
               * HUGE_VAL)
            : (-b / (a + a)));
         const bool v_vertex_in_interval =
           ((INTERSECTION_REAL_ZERO < v_vertex) *
            (v_vertex < INTERSECTION_REAL_0NE));
         const intersection_real v_min_0 =
           intersection_real_min (q0_perp, q2_perp);
         const intersection_real v_max_0 =
           intersection_real_max (q0_perp, q2_perp);
         const intersection_real y_vertex =
           (v_vertex_in_interval)
           ? ((((a * v_vertex) + b) * v_vertex) + c)
           : intersection_real_nan ();       
         const intersection_real v_min =
           ((!v_vertex_in_interval)
            ? v_min_0
            : ((y_vertex < v_min_0)
               ? y_vertex
               : v_min_0));
         const intersection_real v_max =
           ((!v_vertex_in_interval)
            ? v_max_0
            : ((v_max_0 < y_vertex)
               ? y_vertex
               : v_max_0));

         if (!((v_min <= d_max) * (d_min <= v_max)))
           // Q's perpendicular band does not overlap P's fat line.
           return accum0;
         else
           {
             // Find clipped parameter interval [u_min, u_max] on Q.
             intersection_real u_min = INTERSECTION_REAL_ZERO;
             intersection_real u_max = INTERSECTION_REAL_0NE;

             // Intersect quadratic polynomial with d_min and d_max.
             {
               const intersection_real_quadratic_roots qroots =
                 intersection_real_quadratic_roots_solve
                 (intersection_real_quadratic_roots_solve
                  (make_intersection_real_quadratic_roots (),
                   a, b, c - d_min),
                  a, b, c - d_max);

               intersection_real roots[4] = {
                 [0] = qroots.root0,
                 [1] = qroots.root1,
                 [2] = qroots.root2,
                 [3] = qroots.root3
               };
               int nroots = qroots.nroots;
               int k;

               // Candidate parameter bounds from valid roots in
               // [0,1].
               {
                 intersection_real inside_min;
                 intersection_real inside_max;
                 bool has_inside;

                 inside_min = INTERSECTION_REAL_0NE;
                 inside_max = INTERSECTION_REAL_ZERO;
                 has_inside = false;

                 if (q0_perp >= d_min && q0_perp <= d_max)
                   {
                     inside_min = INTERSECTION_REAL_ZERO;
                     inside_max = INTERSECTION_REAL_ZERO;
                     has_inside = true;
                   }
                 if (q2_perp >= d_min && q2_perp <= d_max)
                   {
                     if (!has_inside || INTERSECTION_REAL_0NE < inside_min)
                       inside_min = INTERSECTION_REAL_0NE;
                     if (!has_inside || INTERSECTION_REAL_0NE > inside_max)
                       inside_max = INTERSECTION_REAL_0NE;
                     has_inside = true;
                   }

                 for (k = 0; k != nroots; k += 1)
                   {
                     intersection_real r;
                     r = roots[k];
                     if (r >= INTERSECTION_REAL_ZERO &&
                         r <= INTERSECTION_REAL_0NE)
                       {
                         if (!has_inside || r < inside_min)
                           inside_min = r;
                         if (!has_inside || r > inside_max)
                           inside_max = r;
                         has_inside = true;
                       }
                   }

                 if (has_inside)
                   {
                     u_min = inside_min;
                     u_max = inside_max;
                   }
                 )
                 }

             // If clipping achieved poor progress (reduction < 20%),
             // subdivide the larger curve.
             if ((u_max - u_min) > (intersection_real) 0.8)
               {
                 if (p_span >= q_span)
                   {
                     const intersection_real_quadratic p_left =
                       intersection_real_subdivide_quadratic
                       (p, INTERSECTION_REAL_ZERO,
                        INTERSECTION_REAL_HALF);
                     const intersection_real_quadratic p_right =
                       intersection_real_subdivide_quadratic
                       (p, INTERSECTION_REAL_HALF,
                        INTERSECTION_REAL_0NE);
                     return (intersection_of_two_quadratics
                             (p_right, q, rel_tol, abs_tol, depth + 1,
                              (intersection_of_two_quadratics
                               (p_left, q, rel_tol, abs_tol, depth + 1,
                                accum0))));
                   }
                 else
                   {
                     intersection_real_quadratic q_left =
                       intersection_real_subdivide_quadratic
                       (q, INTERSECTION_REAL_ZERO,
                        INTERSECTION_REAL_HALF);
                     intersection_real_quadratic q_right =
                       intersection_real_subdivide_quadratic
                       (q, INTERSECTION_REAL_HALF,
                        INTERSECTION_REAL_0NE);
                     return (intersection_of_two_quadratics
                             (p, q_right, rel_tol, abs_tol, depth + 1,
                              (intersection_of_two_quadratics
                               (p, q_left, rel_tol, abs_tol, depth + 1,
                                accum0))));
                   }
               }

             // Apply clipped parameter sub-interval to Q
             // and swap roles for rapid convergence.
             {
               intersection_real_quadratic q_clipped;
               // Add slight margin.
               intersection_real m0;
               intersection_real m1;
               m0 = u_min - (intersection_real) 0.01;
               m1 = u_max + (intersection_real) 0.01;
               if (m0 < INTERSECTION_REAL_ZERO)
                 m0 = INTERSECTION_REAL_ZERO;
               if (m1 > INTERSECTION_REAL_0NE)
                 m1 = INTERSECTION_REAL_0NE;
               q_clipped = intersection_real_subdivide_quadratic (q, m0, m1);
               // Swap roles (tail call).
               return intersection_of_two_quadratics (q_clipped, p, rel_tol, abs_tol,
                                                      depth + 1, accum0);
             }
           }
       }
   }

    case 9:
      return accum0;
    default:
      // Programmer error.
      abort ();
    }
}

/* -------------------------------------------------------------------------
 * Cubic Decomposition into Quadratic Arcs in s-power Basis
 * ------------------------------------------------------------------------- */

/* Number of sub-intervals required for cubic degree reduction:
 * E_K = ||Q1_delta - Q2_delta|| / (12 * sqrt(3) * K^2) <= tol */
BEZIER_INTERSECTION_2D_STATICINLINE int
sp_cubic_subdivisions (intersection_real_point q1_delta,
                       intersection_real_point q2_delta,
                       intersection_real tol)
{
  intersection_real diff_x;
  intersection_real diff_y;
  intersection_real diff_norm;
  intersection_real k_sq;
  int k;

  diff_x = q1_delta.x - q2_delta.x;
  diff_y = q1_delta.y - q2_delta.y;
  diff_norm =
    intersection_real_sqrt ((diff_x * diff_x) + (diff_y * diff_y));

  if (diff_norm <= tol || tol <= INTERSECTION_REAL_ZERO)
    {
      return 1;
    }

  /* 12 * sqrt(3) = 20.784609690826528 */
  k_sq = diff_norm / ((intersection_real) 20.78460969 * tol);
  if (k_sq <= INTERSECTION_REAL_0NE)
    {
      return 1;
    }
  k = (int)ceil (intersection_real_sqrt (k_sq));
  if (k > 16)
    {
      k = 16;
    }
  return k;
}

/* Evaluate cubic curve in s-power basis:
 * C(t) = (1-t)*P0 + t*P3 + (1-t)^2*t*Q1_delta + (1-t)*t^2*Q2_delta */
BEZIER_INTERSECTION_2D_STATICINLINE intersection_real_point
sp_eval_cubic (intersection_real_point p0,
               intersection_real_point p3,
               intersection_real_point q1_delta,
               intersection_real_point q2_delta,
               intersection_real t)
{
  intersection_real s;
  intersection_real b1;
  intersection_real b2;
  intersection_real_point pt;

  s = INTERSECTION_REAL_0NE - t;
  b1 = s * s * t;
  b2 = s * t * t;

  pt.x = s * p0.x + t * p3.x + b1 * q1_delta.x + b2 * q2_delta.x;
  pt.y = s * p0.y + t * p3.y + b1 * q1_delta.y + b2 * q2_delta.y;
  return pt;
}

/* Extract optimal approximating quadratic arc for sub-interval [t0, t1] */
BEZIER_INTERSECTION_2D_STATICINLINE intersection_real_quadratic
sp_cubic_to_quad_sub (intersection_real_point p0,
                      intersection_real_point p3,
                      intersection_real_point q1_delta,
                      intersection_real_point q2_delta,
                      intersection_real t0,
                      intersection_real t1)
{
  intersection_real_quadratic q;
  intersection_real_point mid;
  intersection_real_point lin_mid;

  q.t_min = t0;
  q.t_max = t1;
  q.p0 = sp_eval_cubic (p0, p3, q1_delta, q2_delta, t0);
  q.p2 = sp_eval_cubic (p0, p3, q1_delta, q2_delta, t1);

  mid = sp_eval_cubic (p0, p3, q1_delta, q2_delta,
                       INTERSECTION_REAL_HALF * (t0 + t1));
  lin_mid.x = INTERSECTION_REAL_HALF * (q.p0.x + q.p2.x);
  lin_mid.y = INTERSECTION_REAL_HALF * (q.p0.y + q.p2.y);

  q.p1_delta.x = INTERSECTION_REAL_FOUR * (mid.x - lin_mid.x);
  q.p1_delta.y = INTERSECTION_REAL_FOUR * (mid.y - lin_mid.y);

  return q;
}

/* -------------------------------------------------------------------------
 * Primary Public Entry Point (Static Inline, Pure)
 * ------------------------------------------------------------------------- */

BEZIER_INTERSECTION_2D_STATICINLINE intersection_record
bezier_intersection_2d (bool p_quadratic, bool p_spower,
                        intersection_real px0, intersection_real py0,
                        intersection_real px1, intersection_real py1,
                        intersection_real px2, intersection_real py2,
                        intersection_real px3, intersection_real py3,
                        bool q_quadratic, bool q_spower,
                        intersection_real qx0, intersection_real qy0,
                        intersection_real qx1, intersection_real qy1,
                        intersection_real qx2, intersection_real qy2,
                        intersection_real qx3, intersection_real qy3,
                        intersection_real rel_tol)
{
  if (rel_tol <= INTERSECTION_REAL_ZERO)
    return bezier_intersection_2d (p_quadratic, p_spower,
                                   px0, py0, px1, py1,
                                   px2, py2, px3, py3,
                                   q_quadratic, q_spower,
                                   qx0, qy0, qx1, qy1,
                                   qx2, qy2, qx3, qy3,
                                   (intersection_real)
                                   (INTERSECTION_DEFAULT_RELTOL));
  else
    {
      intersection_record accum;
      intersection_real abs_tol;
      intersection_real coord_scale;
      int kp;
      int kq;
      int ip;
      int iq;

      /* S-power canonical vectors */
      intersection_real_point p_pt0;
      intersection_real_point p_pt3;
      intersection_real_point p_q1_delta;
      intersection_real_point p_q2_delta;

      intersection_real_point q_pt0;
      intersection_real_point q_pt3;
      intersection_real_point q_q1_delta;
      intersection_real_point q_q2_delta;

      accum = empty_intersection_record ();

      /* Compute coordinate scale for relative tolerance conversion */
      coord_scale = intersection_real_abs (px0);
      if (intersection_real_abs (py0) > coord_scale)
        coord_scale = intersection_real_abs (py0);
      if (intersection_real_abs (px2) > coord_scale)
        coord_scale = intersection_real_abs (px2);
      if (intersection_real_abs (py2) > coord_scale)
        coord_scale = intersection_real_abs (py2);
      if (intersection_real_abs (qx0) > coord_scale)
        coord_scale = intersection_real_abs (qx0);
      if (intersection_real_abs (qy0) > coord_scale)
        coord_scale = intersection_real_abs (qy0);
      if (intersection_real_abs (qx2) > coord_scale)
        coord_scale = intersection_real_abs (qx2);
      if (intersection_real_abs (qy2) > coord_scale)
        coord_scale = intersection_real_abs (qy2);
      if (coord_scale < INTERSECTION_REAL_0NE)
        coord_scale = INTERSECTION_REAL_0NE;
      abs_tol = rel_tol * coord_scale;

      /* Convert P to canonical s-power form */
      if (p_quadratic)
        {
          p_pt0.x = px0;
          p_pt0.y = py0;
          p_pt3.x = px2;
          p_pt3.y = py2;

          if (p_spower)
            {
              /* px1, py1 is already P1_delta */
              p_q1_delta.x = px1;
              p_q1_delta.y = py1;
            }
          else {
            /* Bernstein control point: P1_delta = 2*P1 - P0 - P2 */
            p_q1_delta.x = INTERSECTION_REAL_TWO * px1 - px0 - px2;
            p_q1_delta.y = INTERSECTION_REAL_TWO * py1 - py0 - py2;
          }
          /* Trivial degree elevation identity: Q1_delta = Q2_delta = P1_delta */
          p_q2_delta = p_q1_delta;
          kp = 1;
        }
      else {
        p_pt0.x = px0;
        p_pt0.y = py0;
        p_pt3.x = px3;
        p_pt3.y = py3;

        if (p_spower)
          {
            p_q1_delta.x = px1;
            p_q1_delta.y = py1;
            p_q2_delta.x = px2;
            p_q2_delta.y = py2;
          }
        else {
          /* Bernstein cubic:
           * Q1_delta = 3*P1 - 2*P0 - P3
           * Q2_delta = 3*P2 - P0 - 2*P3 */
          p_q1_delta.x = INTERSECTION_REAL_THREE * px1 -
            INTERSECTION_REAL_TWO * px0 - px3;
          p_q1_delta.y = INTERSECTION_REAL_THREE * py1 -
            INTERSECTION_REAL_TWO * py0 - py3;
          p_q2_delta.x = INTERSECTION_REAL_THREE * px2 - px0 -
            INTERSECTION_REAL_TWO * px3;
          p_q2_delta.y = INTERSECTION_REAL_THREE * py2 - py0 -
            INTERSECTION_REAL_TWO * py3;
        }
        kp = sp_cubic_subdivisions (p_q1_delta, p_q2_delta, abs_tol);
      }

      /* Convert Q to canonical s-power form */
      if (q_quadratic)
        {
          q_pt0.x = qx0;
          q_pt0.y = qy0;
          q_pt3.x = qx2;
          q_pt3.y = qy2;

          if (q_spower)
            {
              q_q1_delta.x = qx1;
              q_q1_delta.y = qy1;
            }
          else {
            q_q1_delta.x = INTERSECTION_REAL_TWO * qx1 - qx0 - qx2;
            q_q1_delta.y = INTERSECTION_REAL_TWO * qy1 - qy0 - qy2;
          }
          q_q2_delta = q_q1_delta;
          kq = 1;
        }
      else {
        q_pt0.x = qx0;
        q_pt0.y = qy0;
        q_pt3.x = qx3;
        q_pt3.y = qy3;

        if (q_spower)
          {
            q_q1_delta.x = qx1;
            q_q1_delta.y = qy1;
            q_q2_delta.x = qx2;
            q_q2_delta.y = qy2;
          }
        else {
          q_q1_delta.x = (INTERSECTION_REAL_THREE * qx1 -
                          INTERSECTION_REAL_TWO * qx0 - qx3);
          q_q1_delta.y = (INTERSECTION_REAL_THREE * qy1 -
                          INTERSECTION_REAL_TWO * qy0 - qy3);
          q_q2_delta.x = (INTERSECTION_REAL_THREE * qx2 - qx0 -
                          INTERSECTION_REAL_TWO * qx3);
          q_q2_delta.y = (INTERSECTION_REAL_THREE * qy2 - qy0 -
                          INTERSECTION_REAL_TWO * qy3);
        }
        kq = sp_cubic_subdivisions (q_q1_delta, q_q2_delta, abs_tol);
      }

      /* Intersect all pairs of quadratic sub-arcs */
      for (ip = 0; ip < kp && accum.n < 9; ip = ip + 1)
        {
          intersection_real tp0;
          intersection_real tp1;
          intersection_real_quadratic p_sub;

          tp0 = (intersection_real) ip / (intersection_real) kp;
          tp1 = (intersection_real) (ip + 1) / (intersection_real) kp;
          p_sub = sp_cubic_to_quad_sub (p_pt0, p_pt3, p_q1_delta,
                                        p_q2_delta, tp0, tp1);

          for (iq = 0; iq < kq && accum.n < 9; iq = iq + 1)
            {
              intersection_real tq0;
              intersection_real tq1;
              intersection_real_quadratic q_sub;

              tq0 = (intersection_real) iq / (intersection_real) kq;
              tq1 = (intersection_real) (iq + 1) / (intersection_real) kq;
              q_sub = sp_cubic_to_quad_sub (q_pt0, q_pt3, q_q1_delta,
                                            q_q2_delta, tq0, tq1);

              accum = intersection_of_two_quadratics (p_sub, q_sub, rel_tol,
                                                      abs_tol, 0, accum);
            }
        }

      return accum;      
    }
}

#endif /* BEZIER_INTERSECTION_2D_H__INCLUDED__ */

/*
  local variables:
  mode: c
  coding: utf-8
  end:
*/
