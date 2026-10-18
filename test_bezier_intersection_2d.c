/*
 * test_bezier_intersection_2d.c
 * Verification suite for bezier_intersection_2d.h
 */

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "test_bezier_intersection_2d.h"

static void
print_record (const char *label, intersection_record rec)
{
  size_t i;
  printf ("=== %s ===\n", label);
  printf ("found %zu intersection(s)\n", rec.n);
  for (i = 0; i < 9; i = i + 1)
    {
      if (isnan (rec.tp[i]))
        {
          printf ("  [%zu] tp = NaN, tq = NaN\n", i);
        }
      else {
          printf ("  [%zu] tp = %.8f, tq = %.8f\n", i, rec.tp[i],
                  rec.tq[i]);
        }
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
                  r1);
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
    print_record ("Test 2: Quadratic vs Quadratic (s-power input)", r2);
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
    print_record ("Test 3: Cubic vs Quadratic", r3);
  }

  return 0;
}
