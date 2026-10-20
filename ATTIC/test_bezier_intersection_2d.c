/*
 * test_bezier_intersection_2d.c
 * Verification suite for bezier_intersection_2d.h
 */

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "test_bezier_intersection_2d.h"

static intersection_real
evaluate_cubic_spower (const intersection_real c0,
                       const intersection_real c1,
                       const intersection_real c2,
                       const intersection_real c3,
                       const intersection_real t)
{
  const intersection_real s = (intersection_real) 1.0 - t;
  if (t <= (intersection_real) 0.5)
    {
      const intersection_real inner = c1 + (t * (c2 - c1));
      const intersection_real middle = (c3 - c0) + (s * inner);
      return c0 + (t * middle);
    }
  else
    {
      const intersection_real inner = c2 + (s * (c1 - c2));
      const intersection_real middle = (c0 - c3) + (t * inner);
      return c3 + (s * middle);
    }
}

static intersection_real
evaluate_quadratic_spower (const intersection_real c0,
                           const intersection_real c1,
                           const intersection_real c2,
                           const intersection_real not_used,
                           const intersection_real t)
{
  const intersection_real s = (intersection_real) 1.0 - t;
  const intersection_real inner = (c2 - c0) + (c1 * s);
  return c0 + (t * inner);
}

static intersection_real
evaluate_cubic_bernstein (const intersection_real p0,
                          const intersection_real p1,
                          const intersection_real p2,
                          const intersection_real p3,
                          const intersection_real t)
{
  const intersection_real s = (intersection_real) 1.0 - t;
  const intersection_real a0 = (s * p0) + (t * p1);
  const intersection_real a1 = (s * p1) + (t * p2);
  const intersection_real a2 = (s * p2) + (t * p3);
  const intersection_real b0 = (s * a0) + (t * a1);
  const intersection_real b1 = (s * a1) + (t * a2);
  return (s * b0) + (t * b1);
}

static intersection_real
evaluate_quadratic_bernstein (const intersection_real p0,
                              const intersection_real p1,
                              const intersection_real p2,
                              const intersection_real not_used,
                              const intersection_real t)
{
  const intersection_real s = (intersection_real) 1.0 - t;
  const intersection_real a0 = (s * p0) + (t * p1);
  const intersection_real a1 = (s * p1) + (t * p2);
  return (s * a0) + (t * a1);
}

static void
print_record (const char *label, intersection_record rec,
              const intersection_real ax, const intersection_real ay,
              const intersection_real bx, const intersection_real by,
              const intersection_real cx, const intersection_real cy,
              const intersection_real dx, const intersection_real dy,
              intersection_real (*eval_p) (const intersection_real,
                                           const intersection_real,
                                           const intersection_real,
                                           const intersection_real,
                                           const intersection_real),
              const intersection_real ex, const intersection_real ey,
              const intersection_real fx, const intersection_real fy,
              const intersection_real gx, const intersection_real gy,
              const intersection_real hx, const intersection_real hy,
              intersection_real (*eval_q) (const intersection_real,
                                           const intersection_real,
                                           const intersection_real,
                                           const intersection_real,
                                           const intersection_real))
{
  size_t i;
  printf ("=== %s ===\n", label);
  printf ("found %zu intersection(s)\n", rec.n);
  for (i = 0; i < 9; i = i + 1)
    {
      if (isnan (rec.tp[i]))
        printf ("  [%zu] tp = NaN, tq = NaN\n", i);
      else
        printf ("  [%zu] tp = %.8lf (%+lf,%+lf), "
                "tq = %.8lf (%+lf,%+lf)\n",
                i, rec.tp[i],
                eval_p (ax, bx, cx, dx, rec.tp[i]),
                eval_p (ay, by, cy, dy, rec.tp[i]),
                rec.tq[i],
                eval_q (ex, fx, gx, hx, rec.tq[i]),
                eval_q (ey, fy, gy, hy, rec.tq[i]));
    }
}

int
main (void)
{
  /* Test 1: Example 2 from Volume III:
   * P(t) with P0=(-1, 0), P1=(0, 10), P2=(1, 0)
   * Q(t) with Q0=(2, 1), Q1=(-8, 2), Q2=(2, 3)
   * Expected intersection: tp = 0.5, tq = 0.5 at (0, 5) */
  {
    intersection_record r1;
    r1 = bezier_intersection_2d (true, false,
                                 -1.0, 0.0,
                                 0.0, 10.0,
                                 1.0, 0.0,
                                 0.0, 0.0, /* unused */
                                 true, false,
                                 2.0, 1.0,
                                 -8.0, 2.0,
                                 2.0, 3.0,
                                 0.0, 0.0, /* unused */
                                 1e-7);
    print_record ("Test 1: Quadratic vs Quadratic (Bernstein input)",
                  r1,
                  -1.0, 0.0,
                  0.0, 10.0,
                  1.0, 0.0,
                  0.0, 0.0,                /* unused */
                  &evaluate_quadratic_bernstein,
                  2.0, 1.0,
                  -8.0, 2.0,
                  2.0, 3.0,
                  0.0, 0.0,                /* unused */
                  &evaluate_quadratic_bernstein);
  }

  /* Test 2: Same test but passed in s-power form:
   * For P: P0=(-1,0), P2=(1,0), P1_delta = 2*(0,10) - (-1,0) - (1,0) = (0, 20)
   * For Q: Q0=(2,1), Q2=(2,3), Q1_delta = 2*(-8,2) - (2,1) - (2,3) = (-20, 0) */
  {
    intersection_record r2;
    r2 = bezier_intersection_2d (true, true,
                                 -1.0, 0.0,
                                 0.0, 20.0,
                                 1.0, 0.0,
                                 0.0, 0.0,
                                 true, true,
                                 2.0, 1.0,
                                 -20.0, 0.0,
                                 2.0, 3.0,
                                 0.0, 0.0,
                                 1e-7);
    print_record ("Test 2: Quadratic vs Quadratic (s-power input)", r2,
                  -1.0, 0.0,
                  0.0, 20.0,
                  1.0, 0.0,
                  0.0, 0.0,
                  &evaluate_quadratic_spower,
                  2.0, 1.0,
                  -20.0, 0.0,
                  2.0, 3.0,
                  0.0, 0.0,
                  &evaluate_quadratic_spower);
  }

  /* Test 3: Cubic vs Quadratic:
   * Cubic line from (0, -2) to (0, 8) crossed by horizontal quadratic */
  {
    intersection_record r3;
    r3 = bezier_intersection_2d (false, false,
                                 0.0, -2.0,
                                 0.0, 1.0,
                                 0.0, 5.0,
                                 0.0, 8.0,
                                 true, false,
                                 -5.0, 3.0,
                                 0.0, 3.0,
                                 5.0, 3.0,
                                 0.0, 0.0,
                                 1e-7);
    print_record ("Test 3: Cubic vs Quadratic", r3,
                  0.0, -2.0,
                  0.0, 1.0,
                  0.0, 5.0,
                  0.0, 8.0,
                  &evaluate_cubic_bernstein,
                  -5.0, 3.0,
                  0.0, 3.0,
                  5.0, 3.0,
                  0.0, 0.0,
                  &evaluate_quadratic_bernstein);
  }

  return 0;
}
