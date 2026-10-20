pragma wide_character_encoding (utf8);
pragma ada_2022;
pragma assertion_policy (check);

with ada.numerics.elementary_functions;
use ada.numerics.elementary_functions;

with ada.numerics.long_elementary_functions;
use ada.numerics.long_elementary_functions;

package body quadratic_coincidence is

   function "-" (pt1, pt2 : in opoint) return dpoint is
   begin
      return (dx => pt1.x - pt2.x, dy => pt1.y - pt2.y);
   end "-";

   function "+" (opt : in opoint; dpt : in dpoint) return opoint is
   begin
      return (x => opt.x + dpt.dx, y => opt.y + dpt.dy);
   end "+";

   function to_spower_curve
     (ctrl : in control_point_curve) return spower_curve
   is
      curve    : spower_curve;
      p0x, p0y : long_float;
      p1x, p1y : long_float;
      p2x, p2y : long_float;
   begin
      curve.x.c0 := ctrl.pt0.x;
      curve.y.c0 := ctrl.pt0.y;

      curve.x.c2 := ctrl.pt2.x;
      curve.y.c2 := ctrl.pt2.y;

      p0x := long_float (ctrl.pt0.x);
      p0y := long_float (ctrl.pt0.y);
      p1x := long_float (ctrl.pt1.x);
      p1y := long_float (ctrl.pt1.y);
      p2x := long_float (ctrl.pt2.x);
      p2y := long_float (ctrl.pt2.y);
      curve.x.c1 := float (p1x + p1x - p0x - p2x);
      curve.y.c1 := float (p1y + p1y - p0y - p2y);

      return curve;
   end to_spower_curve;

   function to_geocurve (scurve : in spower_curve) return geocurve is
      gcurve : geocurve;
   begin
      gcurve.pt0.x := lreal (scurve.x.c0);
      gcurve.pt0.y := lreal (scurve.y.c0);
      gcurve.pt1.dx := lreal (scurve.x.c1);
      gcurve.pt1.dy := lreal (scurve.y.c1);
      gcurve.pt2.x := lreal (scurve.x.c2);
      gcurve.pt2.y := lreal (scurve.y.c2);
      gcurve.t0 := 0.0;
      gcurve.t1 := 1.0;
      return gcurve;
   end to_geocurve;

   function eval (gcurve : in geocurve; t : in lreal) return opoint is
      pt0 : opoint;
      pt1 : dpoint;
      pt2 : opoint;
   begin
      pt0 := gcurve.pt0;
      pt1 := gcurve.pt1;
      pt2 := gcurve.pt2;
      return
        (if t <= 0.5
         then pt0 + t * ((pt2 - pt0) + ((1.0 - t) * pt1))
         else pt2 + (1.0 - t) * ((pt0 - pt2) + (t * pt1)));
   end eval;

   type affine_projection is record
      dx, dy : lreal;
      a, b   : lreal;
      c, d   : lreal;
   end record;

   function project
     (pt : in opoint; trans : in affine_projection) return opoint
   is
      x, y : lreal;
   begin
      x := pt.x + trans.dx;
      y := pt.y + trans.dy;
      return
        (x => trans.a * x + trans.b * y,
         y => trans.c * x + trans.d * y);
   end project;

   function project
     (pt : in dpoint; trans : in affine_projection) return dpoint
   is
      x, y : lreal;
   begin
      x := pt.dx;
      y := pt.dy;
      return
        (dx => trans.a * x + trans.b * y,
         dy => trans.c * x + trans.d * y);
   end project;

   function skeleton_projection
     (clipper : in geocurve) return affine_projection
   is
      a, b  : lreal;
      denom : lreal;
   begin
      a := clipper.pt2.x - clipper.pt0.x;
      b := clipper.pt2.y - clipper.pt0.y;
      denom := (a * a) + (b * b);
      pragma assert (denom /= 0.0);
      a := a / denom;
      b := b / denom;
      return
        (dx => -clipper.pt0.x,
         dy => -clipper.pt0.y,
         a  => a,
         b  => b,
         c  => -b,
         d  => a);
   end skeleton_projection;

   procedure quadratic_formula
     (a, b, c     : in lreal;
      has_smaller : out boolean;
      smaller     : out lreal;
      has_larger  : out boolean;
      larger      : out lreal)
   with post => has_smaller or not has_larger
   is
      bh, e, d        : lreal;
      arbitrary_value : constant lreal := 6.012214076e23;
   begin
      has_smaller := false;
      has_larger := false;
      smaller := arbitrary_value;
      larger := arbitrary_value;

      if a = 0.0 then
         if b /= 0.0 then
            has_smaller := true;
            smaller := -c / b;
         end if;
      elsif c = 0.0 then
         has_smaller := true;
         has_larger := true;
         smaller := 0.0;
         larger := -b / a;
      else
         bh := b / 2.0;
         if abs (bh) >= abs (c) then
            e := 1.0 - (a / bh) * (c / bh);
            d := sqrt (abs (e)) * abs (bh);
         else
            e := (if c < 0.0 then -a else a);
            e := bh * (bh / abs (c)) - e;
            d := sqrt (abs (e)) * sqrt (abs (c));
         end if;
         if e >= 0.0 then
            has_smaller := true;
            has_larger := true;
            if bh >= 0.0 then
               d := -d;
            end if;
            larger := (-bh + d) / a;
            if larger = 0.0 then
               smaller := 0.0;
            else
               smaller := (c / larger) / a;
            end if;
         end if;
      end if;
   end quadratic_formula;

   procedure check_curve_coincidence
     (c1             : in geocurve;
      c2             : in geocurve;
      tolerance      : in lreal;
      are_coincident : out boolean;
      has_overlap    : out boolean;
      overlap_min_t1 : out lreal;
      overlap_max_t1 : out lreal)
   is
      trans          : affine_projection;
      p1_pt1         : dpoint;
      dx1, dy1       : lreal;
      s0, s1, s2     : opoint;
      t0, t2         : lreal;
      t0_found       : boolean;
      t2_found       : boolean;
      t1_expected    : lreal;
      s1_expected    : opoint;
      err_mid_x      : lreal;
      err_mid_y      : lreal;
      tmin, tmax     : lreal;
      has_sm, has_lg : boolean;
      sm, lg         : lreal;

      procedure find_t_on_c1
        (pt : in opoint; found : out boolean; t_param : out lreal)
      is
         err_sm, err_lg : lreal;
         exp_y          : lreal;
         prev_err       : lreal;
      begin
         found := false;
         t_param := 0.0;
         quadratic_formula
           (-dx1, 1.0 + dx1, -pt.x, has_sm, sm, has_lg, lg);
         if has_sm then
            exp_y := sm * (1.0 - sm) * dy1;
            err_sm := abs (pt.y - exp_y);
            if err_sm <= tolerance then
               found := true;
               t_param := sm;
            end if;
         end if;
         if has_lg then
            exp_y := lg * (1.0 - lg) * dy1;
            err_lg := abs (pt.y - exp_y);
            if err_lg <= tolerance then
               if not found then
                  found := true;
                  t_param := lg;
               else
                  prev_err :=
                    abs (pt.y - (t_param * (1.0 - t_param) * dy1));
                  if err_lg < prev_err then
                     t_param := lg;
                  end if;
               end if;
            end if;
         end if;
      end find_t_on_c1;

   begin
      are_coincident := false;
      has_overlap := false;
      overlap_min_t1 := 0.0;
      overlap_max_t1 := 0.0;

      trans := skeleton_projection (c1);
      p1_pt1 := project (c1.pt1, trans);
      dx1 := p1_pt1.dx;
      dy1 := p1_pt1.dy;

      s0 := project (eval (c2, 0.0), trans);
      s1 := project (eval (c2, 0.5), trans);
      s2 := project (eval (c2, 1.0), trans);

      find_t_on_c1 (s0, t0_found, t0);
      find_t_on_c1 (s2, t2_found, t2);

      if t0_found and t2_found then
         t1_expected := 0.5 * (t0 + t2);
         s1_expected.x :=
           t1_expected + t1_expected * (1.0 - t1_expected) * dx1;
         s1_expected.y := t1_expected * (1.0 - t1_expected) * dy1;
         err_mid_x := abs (s1.x - s1_expected.x);
         err_mid_y := abs (s1.y - s1_expected.y);
         if err_mid_x <= tolerance and err_mid_y <= tolerance then
            are_coincident := true;
            tmin := lreal'min (t0, t2);
            tmax := lreal'max (t0, t2);
            overlap_min_t1 := lreal'max (0.0, tmin);
            overlap_max_t1 := lreal'min (1.0, tmax);
            if overlap_min_t1 <= overlap_max_t1 + tolerance then
               has_overlap := true;
            end if;
         end if;
      end if;
   end check_curve_coincidence;

   function are_curves_coincident
     (c1        : in geocurve;
      c2        : in geocurve;
      tolerance : in lreal := 1.0e-5) return boolean
   is
      coincident   : boolean;
      overlap      : boolean;
      t_min, t_max : lreal;
   begin
      check_curve_coincidence
        (c1, c2, tolerance, coincident, overlap, t_min, t_max);
      return coincident;
   end are_curves_coincident;

   function are_curves_coincident
     (p         : in spower_curve;
      q         : in spower_curve;
      tolerance : in lreal := 1.0e-5) return boolean is
   begin
      return
        are_curves_coincident
          (to_geocurve (p), to_geocurve (q), tolerance);
   end are_curves_coincident;

end quadratic_coincidence;
