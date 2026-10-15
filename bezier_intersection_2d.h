/*
 * Intersection of quadratic and cubic parametric splines in the
 * plane, using the symmetric power (s-power) basis and fat-line
 * Bézier clipping.
 *
 * Written in C99 as pure static inline routines in tail-recursive
 * style.
 *
 * References
 * ----------
 *
 *   Javier Sánchez-Reyes. 1997. The symmetric analogue of the
 *   polynomial power basis. ACM Trans. Graph. 16, 3 (July 1997),
 *   319–357. DOI: https://doi.org/10.1145/256157.256162
 *
 *   Javier Sánchez-Reyes. 2000. Applications of the polynomial
 *   s-power basis in geometry processing. ACM Trans. Graph. 19, 1
 *   (Jan. 2000), 27–55. DOI: https://doi.org/10.1145/343002.343018
 *
 *   Thomas W. Sederberg and Tomoyuki Nishita. 1990. Curve
 *   intersection using Bézier clipping. Comput. Aided Des. 22, 9
 *   (Nov. 1990), 538–549.
 *   https://doi.org/10.1016/0010-4485(90)90039-F
 *
 */

#ifndef BEZIER_INTERSECTION_2D_H__INCLUDED__
#define BEZIER_INTERSECTION_2D_H__INCLUDED__

#include <stdbool.h>
#include <stddef.h>
#include <math.h>

#ifndef INTERSECTION_REAL
#define INTERSECTION_REAL double
#endif

typedef INTERSECTION_REAL intersection_real;

typedef struct
{
  bool degenerate;
  size_t n;
  intersection_real tp[9];  /* p parameters */
  intersection_real tq[9];  /* q parameters */
} intersection_record;

/* A point in the plane. */
typedef struct
{
  intersection_real x;
  intersection_real y;
} intersection_real_point;

/* Quadratic curve in s-power basis:
 *
 *     C(t) = (1 - t)*P0 + t*P2 + (1 - t)*t*P1_delta
 *
 *       for t_min <= t <= t_max
 */
typedef struct
{
  intersection_real_point p0;
  intersection_real_point p2;
  intersection_real_point p1_delta;
  intersection_real t_min;
  intersection_real t_max;
} intersection_real_quad;

static inline intersection_real
sp_nan (void)
{
#if defined(NAN)
  return (intersection_real) NAN;
#else
  return (intersection_real) nan ("");
#endif
}

static inline intersection_real
intersection_real_bound_0_1 (intersection_real x)
{
  return ((x < (intersection_real) 0.0)
          ? ((intersection_real) 0.0)
          : (((intersection_real) 1.0 < x)
             ? ((intersection_real) 1.0)
             : x));
}

static inline intersection_record
// FIXME: GET RID OF THIS ..................................
empty_intersection_record (void)
{
  intersection_record result = {
    .degenerate = false,
    .n = 0,

    .tp[0] = sp_nan (),
    .tp[1] = sp_nan (),
    .tp[2] = sp_nan (),
    .tp[3] = sp_nan (),
    .tp[4] = sp_nan (),
    .tp[5] = sp_nan (),
    .tp[6] = sp_nan (),
    .tp[7] = sp_nan (),
    .tp[8] = sp_nan (),

    .tq[0] = sp_nan (),
    .tq[1] = sp_nan (),
    .tq[2] = sp_nan (),
    .tq[3] = sp_nan (),
    .tq[4] = sp_nan (),
    .tq[5] = sp_nan (),
    .tq[6] = sp_nan (),
    .tq[7] = sp_nan (),
    .tq[8] = sp_nan ()
  };
  return result;
}

static inline bool
intersection_record_contains (const intersection_record *const rec,
                              const intersection_real tp,
                              const intersection_real tq,
                              const intersection_real tol,
                              const size_t i)
{
  if (i == rec->n)
    return false;
  else
    {
      const intersection_real dp =
        (tp < rec->tp[i]) ? (rec->tp[i] - tp) : (tp - rec->tp[i]);
      const intersection_real dq =
        (tq < rec->tq[i]) ? (rec->tq[i] - tq) : (tq - rec->tq[i]);
      if (dp <= tol && dq <= tol)
        return true;
      else
        return intersection_record_contains (rec, tp, tq, tol, i + 1);
    }
}

static inline intersection_record
intersection_record_add (intersection_record rec,
                         intersection_real tp,
                         intersection_real tq,
                         intersection_real tol)
{
  tp = intersection_real_bound_0_1 (tp);
  tq = intersection_real_bound_0_1 (tq);

  if (intersection_record_contains (&rec, tp, tq, tol, 0))
    return rec;
  else if (rec.n < 9)
    {
      rec.tp[rec.n] = tp;
      rec.tq[rec.n] = tq;
      rec.n = rec.n + 1;
      return rec;
    }
  else
    return rec;
}

static inline intersection_real
sp_dist_sq (intersection_real_point a, intersection_real_point b)
{
  intersection_real dx;
  intersection_real dy;
  dx = b.x - a.x;
  dy = b.y - a.y;
  return dx * dx + dy * dy;
}

static inline intersection_real_point
sp_eval_quad (const intersection_real_quad *q, intersection_real u)
{
  intersection_real s;
  intersection_real blend;
  intersection_real_point pt;
  s = (intersection_real) 1.0 - u;
  blend = s * u;
  pt.x = s * q->p0.x + u * q->p2.x + blend * q->p1_delta.x;
  pt.y = s * q->p0.y + u * q->p2.y + blend * q->p1_delta.y;
  return pt;
}

static inline intersection_real_quad
sp_subdivide_quad (const intersection_real_quad *q,
                   intersection_real u0,
                   intersection_real u1)
{
  intersection_real_quad sub;
  intersection_real du;
  intersection_real orig_span;
  intersection_real_point mid;
  intersection_real_point lin_mid;

  du = u1 - u0;
  orig_span = q->t_max - q->t_min;

  sub.t_min = q->t_min + u0 * orig_span;
  sub.t_max = q->t_min + u1 * orig_span;

  sub.p0 = sp_eval_quad (q, u0);
  sub.p2 = sp_eval_quad (q, u1);

  /* Midpoint evaluation to obtain exact s-power deviation vector:
   * P1_delta = 4 * ( C((u0+u1)/2) - 0.5*(P0 + P2) ) */
  mid = sp_eval_quad (q, (intersection_real) 0.5 * (u0 + u1));
  lin_mid.x = (intersection_real) 0.5 * (sub.p0.x + sub.p2.x);
  lin_mid.y = (intersection_real) 0.5 * (sub.p0.y + sub.p2.y);

  sub.p1_delta.x = (intersection_real) 4.0 * (mid.x - lin_mid.x);
  sub.p1_delta.y = (intersection_real) 4.0 * (mid.y - lin_mid.y);

  return sub;
}

/* -------------------------------------------------------------------------
 * Fat-Line Clipping for Two Quadratic Arcs (Tail Recursive)
 * ------------------------------------------------------------------------- */

static inline intersection_record
sp_intersect_quad_pair_tail (intersection_real_quad p,
                             intersection_real_quad q,
                             intersection_real rel_tol,
                             intersection_real abs_tol,
                             int depth,
                             intersection_record acc)
{
  intersection_real p_span;
  intersection_real q_span;
  intersection_real skel_dx;
  intersection_real skel_dy;
  intersection_real skel_len_sq;
  intersection_real skel_len;
  intersection_real d_y;
  intersection_real d_min;
  intersection_real d_max;
  intersection_real q0_perp;
  intersection_real q2_perp;
  intersection_real qd_perp;
  intersection_real a;
  intersection_real b;
  intersection_real c;
  intersection_real v_min;
  intersection_real v_max;
  intersection_real u_min;
  intersection_real u_max;

  if (acc.n >= 9)
    {
      return acc;
    }

  p_span = p.t_max - p.t_min;
  q_span = q.t_max - q.t_min;

  /* If both parameter spans are within relative tolerance, record root */
  if (p_span <= rel_tol && q_span <= rel_tol)
    {
      intersection_real tp_cand;
      intersection_real tq_cand;
      tp_cand = (intersection_real) 0.5 * (p.t_min + p.t_max);
      tq_cand = (intersection_real) 0.5 * (q.t_min + q.t_max);
      return intersection_record_add (acc, tp_cand, tq_cand,
                                      rel_tol * (intersection_real) 2.0);
    }

  /* Limit recursion depth to prevent infinite loops */
  if (depth > 60)
    {
      intersection_real tp_cand;
      intersection_real tq_cand;
      tp_cand = (intersection_real) 0.5 * (p.t_min + p.t_max);
      tq_cand = (intersection_real) 0.5 * (q.t_min + q.t_max);
      return intersection_record_add (acc, tp_cand, tq_cand,
                                      rel_tol * (intersection_real) 2.0);
    }

  /* Skeleton baseline vector of curve P */
  skel_dx = p.p2.x - p.p0.x;
  skel_dy = p.p2.y - p.p0.y;
  skel_len_sq = skel_dx * skel_dx + skel_dy * skel_dy;

  /* If curve P's skeleton has collapsed to a point, swap roles or split */
  if (skel_len_sq <= abs_tol * abs_tol)
    {
      intersection_real q_skel_dx;
      intersection_real q_skel_dy;
      intersection_real q_skel_len_sq;
      q_skel_dx = q.p2.x - q.p0.x;
      q_skel_dy = q.p2.y - q.p0.y;
      q_skel_len_sq = q_skel_dx * q_skel_dx + q_skel_dy * q_skel_dy;

      if (q_skel_len_sq <= abs_tol * abs_tol)
        {
          /* Both curves collapsed to points */
          if (sp_dist_sq (p.p0, q.p0) <= abs_tol * abs_tol)
            {
              return intersection_record_add (acc, p.t_min, q.t_min,
                                              rel_tol *
                                              (intersection_real) 2.0);
            }
          return acc;
        }
      /* Swap P and Q */
      return sp_intersect_quad_pair_tail (q, p, rel_tol, abs_tol,
                                          depth + 1, acc);
    }

  skel_len = sqrt (skel_len_sq);

  /* Transverse deviation of P from its skeleton line:
   * d_y = ( -skel_dy * p1_delta.x + skel_dx * p1_delta.y ) / skel_len */
  d_y = (-skel_dy * p.p1_delta.x + skel_dx * p.p1_delta.y) / skel_len;

  /* Extreme bounds of P along perpendicular axis:
   * The polynomial is (1-t)*0 + t*0 + (1-t)*t*d_y = (1-t)*t*d_y */
  if (d_y >= (intersection_real) 0.0)
    {
      d_min = -(intersection_real) 1e-12;
      d_max = (intersection_real) 0.25 * d_y + (intersection_real) 1e-12;
    }
  else {
    d_min = (intersection_real) 0.25 * d_y - (intersection_real) 1e-12;
    d_max = (intersection_real) 1e-12;
  }

  /* Project Q onto perpendicular distance to P's skeleton line */
  q0_perp = (-skel_dy * (q.p0.x - p.p0.x) + skel_dx *
             (q.p0.y - p.p0.y)) / skel_len;
  q2_perp = (-skel_dy * (q.p2.x - p.p0.x) + skel_dx *
             (q.p2.y - p.p0.y)) / skel_len;
  qd_perp = (-skel_dy * q.p1_delta.x + skel_dx * q.p1_delta.y) /
    skel_len;

  /* Q's perpendicular coordinate polynomial:
   * y(v) = (1 - v)*q0_perp + v*q2_perp + (1 - v)*v*qd_perp
   *      = a*v^2 + b*v + c */
  a = -qd_perp;
  b = q2_perp - q0_perp + qd_perp;
  c = q0_perp;

  /* Compute range of y(v) for v in [0, 1] */
  v_min = (q0_perp < q2_perp) ? q0_perp : q2_perp;
  v_max = (q0_perp > q2_perp) ? q0_perp : q2_perp;

  if (fabs (a) > (intersection_real) 1e-14)
    {
      intersection_real v_vertex;
      v_vertex = -b / ((intersection_real) 2.0 * a);
      if (v_vertex > (intersection_real) 0.0 &&
          v_vertex < (intersection_real) 1.0)
        {
          intersection_real y_vertex;
          y_vertex = a * v_vertex * v_vertex + b * v_vertex + c;
          if (y_vertex < v_min)
            {
              v_min = y_vertex;
            }
          if (y_vertex > v_max)
            {
              v_max = y_vertex;
            }
        }
    }

  /* Disjoint test: if Q's perpendicular band does not overlap P's fat line */
  if (v_min > d_max || v_max < d_min)
    {
      return acc;
    }

  /* Find clipped parameter interval [u_min, u_max] on Q */
  u_min = (intersection_real) 0.0;
  u_max = (intersection_real) 1.0;

  /* Intersect quadratic polynomial with d_min and d_max */
  {
    intersection_real roots[4];
    int nroots;
    int k;
    nroots = 0;

    /* Solve a*v^2 + b*v + (c - d_min) = 0 */
    {
      intersection_real c_sub;
      intersection_real disc;
      c_sub = c - d_min;
      if (fabs (a) < (intersection_real) 1e-14)
        {
          if (fabs (b) > (intersection_real) 1e-14)
            {
              roots[nroots] = -c_sub / b;
              nroots = nroots + 1;
            }
        }
      else {
        disc = b * b - (intersection_real) 4.0 * a * c_sub;
        if (disc >= (intersection_real) 0.0)
          {
            intersection_real s_disc;
            s_disc = sqrt (disc);
            roots[nroots] = (-b - s_disc) /
              ((intersection_real) 2.0 * a);
            nroots = nroots + 1;
            roots[nroots] = (-b + s_disc) /
              ((intersection_real) 2.0 * a);
            nroots = nroots + 1;
          }
      }
    }

    /* Solve a*v^2 + b*v + (c - d_max) = 0 */
    {
      intersection_real c_sub;
      intersection_real disc;
      c_sub = c - d_max;
      if (fabs (a) < (intersection_real) 1e-14)
        {
          if (fabs (b) > (intersection_real) 1e-14)
            {
              roots[nroots] = -c_sub / b;
              nroots = nroots + 1;
            }
        }
      else {
        disc = b * b - (intersection_real) 4.0 * a * c_sub;
        if (disc >= (intersection_real) 0.0)
          {
            intersection_real s_disc;
            s_disc = sqrt (disc);
            roots[nroots] = (-b - s_disc) /
              ((intersection_real) 2.0 * a);
            nroots = nroots + 1;
            roots[nroots] = (-b + s_disc) /
              ((intersection_real) 2.0 * a);
            nroots = nroots + 1;
          }
      }
    }

    /* Candidate parameter bounds from valid roots in [0, 1] */
    {
      intersection_real inside_min;
      intersection_real inside_max;
      bool has_inside;

      inside_min = (intersection_real) 1.0;
      inside_max = (intersection_real) 0.0;
      has_inside = false;

      if (q0_perp >= d_min && q0_perp <= d_max)
        {
          inside_min = (intersection_real) 0.0;
          inside_max = (intersection_real) 0.0;
          has_inside = true;
        }
      if (q2_perp >= d_min && q2_perp <= d_max)
        {
          if (!has_inside || (intersection_real) 1.0 < inside_min)
            {
              inside_min = (intersection_real) 1.0;
            }
          if (!has_inside || (intersection_real) 1.0 > inside_max)
            {
              inside_max = (intersection_real) 1.0;
            }
          has_inside = true;
        }

      for (k = 0; k < nroots; k = k + 1)
        {
          intersection_real r;
          r = roots[k];
          if (r >= (intersection_real) 0.0 &&
              r <= (intersection_real) 1.0)
            {
              if (!has_inside || r < inside_min)
                {
                  inside_min = r;
                }
              if (!has_inside || r > inside_max)
                {
                  inside_max = r;
                }
              has_inside = true;
            }
        }

      if (has_inside)
        {
          u_min = inside_min;
          u_max = inside_max;
        }
    }
  }

  /* If clipping achieved poor progress (reduction < 20%), subdivide the larger curve */
  if ((u_max - u_min) > (intersection_real) 0.8)
    {
      if (p_span >= q_span)
        {
          intersection_real_quad p_left;
          intersection_real_quad p_right;
          p_left = sp_subdivide_quad (&p, (intersection_real) 0.0,
                                      (intersection_real) 0.5);
          p_right = sp_subdivide_quad (&p, (intersection_real) 0.5,
                                       (intersection_real) 1.0);

          acc = sp_intersect_quad_pair_tail (p_left, q, rel_tol,
                                             abs_tol, depth + 1, acc);
          return sp_intersect_quad_pair_tail (p_right, q, rel_tol,
                                              abs_tol, depth + 1, acc);
        }
      else {
        intersection_real_quad q_left;
        intersection_real_quad q_right;
        q_left = sp_subdivide_quad (&q, (intersection_real) 0.0,
                                    (intersection_real) 0.5);
        q_right = sp_subdivide_quad (&q, (intersection_real) 0.5,
                                     (intersection_real) 1.0);

        acc = sp_intersect_quad_pair_tail (p, q_left, rel_tol,
                                           abs_tol, depth + 1, acc);
        return sp_intersect_quad_pair_tail (p, q_right, rel_tol,
                                            abs_tol, depth + 1, acc);
      }
    }

  /* Apply clipped parameter sub-interval to Q and swap roles for rapid convergence */
  {
    intersection_real_quad q_clipped;
    /* Add slight margin */
    intersection_real m0;
    intersection_real m1;
    m0 = u_min - (intersection_real) 0.01;
    m1 = u_max + (intersection_real) 0.01;
    if (m0 < (intersection_real) 0.0)
      {
        m0 = (intersection_real) 0.0;
      }
    if (m1 > (intersection_real) 1.0)
      {
        m1 = (intersection_real) 1.0;
      }
    q_clipped = sp_subdivide_quad (&q, m0, m1);
    /* Swap roles (tail call) */
    return sp_intersect_quad_pair_tail (q_clipped, p, rel_tol, abs_tol,
                                        depth + 1, acc);
  }
}

/* -------------------------------------------------------------------------
 * Cubic Decomposition into Quadratic Arcs in s-power Basis
 * ------------------------------------------------------------------------- */

/* Number of sub-intervals required for cubic degree reduction:
 * E_K = ||Q1_delta - Q2_delta|| / (12 * sqrt(3) * K^2) <= tol */
static inline int
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
  diff_norm = sqrt (diff_x * diff_x + diff_y * diff_y);

  if (diff_norm <= tol || tol <= (intersection_real) 0.0)
    {
      return 1;
    }

  /* 12 * sqrt(3) = 20.784609690826528 */
  k_sq = diff_norm / ((intersection_real) 20.78460969 * tol);
  if (k_sq <= (intersection_real) 1.0)
    {
      return 1;
    }
  k = (int)ceil (sqrt (k_sq));
  if (k > 16)
    {
      k = 16;
    }
  return k;
}

/* Evaluate cubic curve in s-power basis:
 * C(t) = (1-t)*P0 + t*P3 + (1-t)^2*t*Q1_delta + (1-t)*t^2*Q2_delta */
static inline intersection_real_point
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

  s = (intersection_real) 1.0 - t;
  b1 = s * s * t;
  b2 = s * t * t;

  pt.x = s * p0.x + t * p3.x + b1 * q1_delta.x + b2 * q2_delta.x;
  pt.y = s * p0.y + t * p3.y + b1 * q1_delta.y + b2 * q2_delta.y;
  return pt;
}

/* Extract optimal approximating quadratic arc for sub-interval [t0, t1] */
static inline intersection_real_quad
sp_cubic_to_quad_sub (intersection_real_point p0,
                      intersection_real_point p3,
                      intersection_real_point q1_delta,
                      intersection_real_point q2_delta,
                      intersection_real t0,
                      intersection_real t1)
{
  intersection_real_quad q;
  intersection_real_point mid;
  intersection_real_point lin_mid;

  q.t_min = t0;
  q.t_max = t1;
  q.p0 = sp_eval_cubic (p0, p3, q1_delta, q2_delta, t0);
  q.p2 = sp_eval_cubic (p0, p3, q1_delta, q2_delta, t1);

  mid = sp_eval_cubic (p0, p3, q1_delta, q2_delta,
                       (intersection_real) 0.5 * (t0 + t1));
  lin_mid.x = (intersection_real) 0.5 * (q.p0.x + q.p2.x);
  lin_mid.y = (intersection_real) 0.5 * (q.p0.y + q.p2.y);

  q.p1_delta.x = (intersection_real) 4.0 * (mid.x - lin_mid.x);
  q.p1_delta.y = (intersection_real) 4.0 * (mid.y - lin_mid.y);

  return q;
}

/* -------------------------------------------------------------------------
 * Primary Public Entry Point (Static Inline, Pure)
 * ------------------------------------------------------------------------- */

static inline intersection_record
spower_intersect (bool p_quadratic, bool p_spower,
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
  if (rel_tol <= (intersection_real) 0.0)
    return spower_intersect (p_quadratic, p_spower,
                             px0, py0,
                             px1, py1,
                             px2, py2,
                             px3, py3,
                             q_quadratic, q_spower,
                             qx0, qy0,
                             qx1, qy1,
                             qx2, qy2,
                             qx3, qy3,
                             (intersection_real)  1e-7);
  else
    {
      intersection_record acc;
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

      acc = empty_intersection_record ();

      /* Compute coordinate scale for relative tolerance conversion */
      coord_scale = fabs (px0);
      if (fabs (py0) > coord_scale)
        coord_scale = fabs (py0);
      if (fabs (px2) > coord_scale)
        coord_scale = fabs (px2);
      if (fabs (py2) > coord_scale)
        coord_scale = fabs (py2);
      if (fabs (qx0) > coord_scale)
        coord_scale = fabs (qx0);
      if (fabs (qy0) > coord_scale)
        coord_scale = fabs (qy0);
      if (fabs (qx2) > coord_scale)
        coord_scale = fabs (qx2);
      if (fabs (qy2) > coord_scale)
        coord_scale = fabs (qy2);
      if (coord_scale < (intersection_real) 1.0)
        {
          coord_scale = (intersection_real) 1.0;
        }
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
            p_q1_delta.x = (intersection_real) 2.0 * px1 - px0 - px2;
            p_q1_delta.y = (intersection_real) 2.0 * py1 - py0 - py2;
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
          p_q1_delta.x = (intersection_real) 3.0 * px1 -
            (intersection_real) 2.0 * px0 - px3;
          p_q1_delta.y = (intersection_real) 3.0 * py1 -
            (intersection_real) 2.0 * py0 - py3;
          p_q2_delta.x = (intersection_real) 3.0 * px2 - px0 -
            (intersection_real) 2.0 * px3;
          p_q2_delta.y = (intersection_real) 3.0 * py2 - py0 -
            (intersection_real) 2.0 * py3;
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
            q_q1_delta.x = (intersection_real) 2.0 * qx1 - qx0 - qx2;
            q_q1_delta.y = (intersection_real) 2.0 * qy1 - qy0 - qy2;
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
          q_q1_delta.x = ((intersection_real) 3.0 * qx1 -
                          (intersection_real) 2.0 * qx0 - qx3);
          q_q1_delta.y = ((intersection_real) 3.0 * qy1 -
                          (intersection_real) 2.0 * qy0 - qy3);
          q_q2_delta.x = ((intersection_real) 3.0 * qx2 - qx0 -
                          (intersection_real) 2.0 * qx3);
          q_q2_delta.y = ((intersection_real) 3.0 * qy2 - qy0 -
                          (intersection_real) 2.0 * qy3);
        }
        kq = sp_cubic_subdivisions (q_q1_delta, q_q2_delta, abs_tol);
      }

      /* Intersect all pairs of quadratic sub-arcs */
      for (ip = 0; ip < kp && acc.n < 9; ip = ip + 1)
        {
          intersection_real tp0;
          intersection_real tp1;
          intersection_real_quad p_sub;

          tp0 = (intersection_real) ip / (intersection_real) kp;
          tp1 = (intersection_real) (ip + 1) / (intersection_real) kp;
          p_sub = sp_cubic_to_quad_sub (p_pt0, p_pt3, p_q1_delta,
                                        p_q2_delta, tp0, tp1);

          for (iq = 0; iq < kq && acc.n < 9; iq = iq + 1)
            {
              intersection_real tq0;
              intersection_real tq1;
              intersection_real_quad q_sub;

              tq0 = (intersection_real) iq / (intersection_real) kq;
              tq1 = (intersection_real) (iq + 1) / (intersection_real) kq;
              q_sub = sp_cubic_to_quad_sub (q_pt0, q_pt3, q_q1_delta,
                                            q_q2_delta, tq0, tq1);

              acc = sp_intersect_quad_pair_tail (p_sub, q_sub, rel_tol,
                                                 abs_tol, 0, acc);
            }
        }

      return acc;
    }
}

#endif /* BEZIER_INTERSECTION_2D_H__INCLUDED__ */
