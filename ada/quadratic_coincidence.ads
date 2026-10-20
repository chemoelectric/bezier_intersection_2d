pragma wide_character_encoding (utf8);
pragma ada_2022;
pragma assertion_policy (check);

package quadratic_coincidence is

   subtype real is float;
   subtype lreal is long_float;

   type point is record
      x, y : real;
   end record;

   type control_point_curve is record
      pt0, pt1, pt2 : point;
   end record;

   type spower_poly is record
      c0, c1, c2 : real;
   end record;

   type spower_curve is record
      x, y : spower_poly;
   end record;

   type opoint is record
      x, y : lreal;
   end record;

   type dpoint is record
      dx, dy : lreal;
   end record;

   type geocurve is record
      pt0    : opoint;
      pt1    : dpoint;
      pt2    : opoint;
      t0, t1 : lreal;
   end record;

   function to_spower_curve
     (ctrl : in control_point_curve) return spower_curve;

   function to_geocurve (scurve : in spower_curve) return geocurve;

   function eval (gcurve : in geocurve; t : in lreal) return opoint;

   procedure check_curve_coincidence
     (c1             : in geocurve;
      c2             : in geocurve;
      tolerance      : in lreal;
      are_coincident : out boolean;
      has_overlap    : out boolean;
      overlap_min_t1 : out lreal;
      overlap_max_t1 : out lreal)
   with
     pre  => tolerance > 0.0,
     post => (if not are_coincident then not has_overlap);

   function are_curves_coincident
     (c1        : in geocurve;
      c2        : in geocurve;
      tolerance : in lreal := 1.0e-5) return boolean;

   function are_curves_coincident
     (p         : in spower_curve;
      q         : in spower_curve;
      tolerance : in lreal := 1.0e-5) return boolean;

end quadratic_coincidence;
