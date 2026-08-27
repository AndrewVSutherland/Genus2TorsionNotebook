// claude_v7_stratum2.m — derive the [3,1,2,1]-pattern stratum of the order-7
// D_inf locus on the [1,1,1,1,2]-shape chart f = x(x-1)(x-a)(x-b)(x^2+cx+d).
// This is the stratum CONTAINING the RM anchor 152100.eb.2 (all 12 chart
// normalizations show pattern [3,1,2,1]).  Staged derivation:
//   Stage A: CF steps 0..2 over Q(a,b,c,d) (cheap); E0 := lc(Q_2) — the
//            deg-drop condition (deg Q_2 = 1 <=> a_2 has degree 2).
//   Stage B: solve E0 = 0 for one variable, substitute, continue the CF with
//            pattern [.,.,2,1]; closure Q_4 = const gives the last equation.
// Every rational point off the boundary => torsion >= (Z/2)^3 x Z/7 = [2,2,14].
SetColumns(0);
SetMemoryLimit(24*10^9);
Q := Rationals();
F<a,b,c,d> := FunctionField(Q, 4);
P<x> := PolynomialRing(F);

f := x*(x-1)*(x-a)*(x-b)*(x^2+c*x+d);
s0 := x^3;
for k in [1..3] do
  dd := f - s0^2;
  if Degree(dd) le 2 then break; end if;
  s0 := s0 + (Coefficient(dd, 6-k)/(2*Coefficient(s0,3)))*x^(3-k);
end for;

// steps 0,1: generic
Pi := P!0; Qi := P!1;
for i in [0..1] do
  ai := (Pi + s0) div Qi;
  printf "step %o: deg a=%o deg Q=%o\n", i, Degree(ai), Degree(Qi);
  Pn := ai*Qi - Pi;
  Qn := (f - Pn^2) div Qi;
  Pi := Pn; Qi := Qn;
end for;
// now Qi = Q_2 (deg 2 generically); E0 = its leading coefficient
printf "deg Q_2 = %o\n", Degree(Qi);
E0 := Numerator(F!Coefficient(Qi, 2));
R4<aa,bb,cc,dd4> := PolynomialRing(Q, 4);
h4 := hom< Parent(E0) -> R4 | [aa,bb,cc,dd4] >;
G0 := h4(E0);
printf "E0: %o terms, degrees (a,b,c,d) = (%o,%o,%o,%o)\n",
  #Terms(G0), Degree(G0,1), Degree(G0,2), Degree(G0,3), Degree(G0,4);

// anchor sanity: all 12 normalizations satisfy E0 = 0
anchors := [
 [-3/5,16/15,-4/15,-52/75], [-5/3,-16/9,4/9,-52/27], [15/16,-9/16,-1/4,-39/64],
 [8/5,-1/15,-26/15,1/25], [5/8,-1/24,-13/12,1/64], [-15,-24,26,9],
 [8/3,25/9,-22/9,-13/27], [3/8,25/24,-11/12,-13/192], [9/25,24/25,-22/25,-39/625],
 [1/16,25/16,-7/4,9/64], [16,25,-28,36], [16/25,1/25,-28/25,36/625] ];
nok := 0;
for A in anchors do
  if Evaluate(G0, A) eq 0 then nok +:= 1; end if;
end for;
printf "anchors on E0=0: %o / 12\n", nok;

// stage B: solve E0 for d if possible
dd_deg := Degree(G0, 4);
printf "deg_d E0 = %o\n", dd_deg;
if dd_deg eq 1 then
  // d = -coeff0/coeff1 as rational function of (a,b,c)
  c1 := Coefficient(G0, 4, 1); c0 := Coefficient(G0, 4, 0);
  F3<a3,b3,c3> := FunctionField(Q, 3);
  h3n := hom< R4 -> F3 | [a3,b3,c3,0] >;   // only valid on d-free parts
  dsol := - F3!h3n(c0) / F3!h3n(c1);
  printf "d = %o\n", dsol;
  // substituted chart over Q(a,b,c)
  P3<x3> := PolynomialRing(F3);
  f3 := x3*(x3-1)*(x3-a3)*(x3-b3)*(x3^2+c3*x3+dsol);
  s3 := x3^3;
  for k in [1..3] do
    dd := f3 - s3^2;
    if Degree(dd) le 2 then break; end if;
    s3 := s3 + (Coefficient(dd, 6-k)/(2*Coefficient(s3,3)))*x3^(3-k);
  end for;
  Pi3 := P3!0; Qi3 := P3!1;
  pat := [];
  for i in [0..3] do
    ai := (Pi3 + s3) div Qi3;
    Append(~pat, Degree(ai));
    printf "B step %o: deg a=%o deg Q=%o\n", i, Degree(ai), Degree(Qi3);
    Pn := ai*Qi3 - Pi3;
    Qn := (f3 - Pn^2) div Qi3;
    Pi3 := Pn; Qi3 := Qn;
  end for;
  printf "pattern so far %o; deg Q_4 = %o (want 0 for closure)\n", pat, Degree(Qi3);
  if Degree(Qi3) ge 1 then
    Ec := Numerator(F3!Coefficient(Qi3, 1));
    R3<ar,br,cr> := PolynomialRing(Q, 3);
    h3 := hom< Parent(Ec) -> R3 | [ar,br,cr] >;
    Gc := h3(Ec);
    printf "CLOSURE Ec: %o terms, degrees (a,b,c) = (%o,%o,%o)\n",
      #Terms(Gc), Degree(Gc,1), Degree(Gc,2), Degree(Gc,3);
    Write("data/claude_v7_stratum2_eqs.txt",
      Sprintf("dsol := %o;\nEc := %o;\n", dsol, Gc) : Overwrite := true);
    // anchor check on the closure surface: need (a,b,c) with d matching dsol
    for A in anchors do
      dv := Evaluate(Numerator(dsol), [A[1],A[2],A[3]]) / Evaluate(Denominator(dsol), [A[1],A[2],A[3]]);
      if dv eq A[4] then
        va := Evaluate(Gc, [A[1],A[2],A[3]]);
        printf "anchor %o: dsol matches (%o), Ec = %o\n", A, dv eq A[4], va eq 0 select "ZERO" else "nonzero";
      end if;
    end for;
    // also check deg-2 coefficient of Q4 (should vanish identically or be dependent)
    if Degree(Qi3) ge 2 then
      Ec2 := Numerator(F3!Coefficient(Qi3, 2));
      printf "deg-2 coeff of Q4 also present: %o terms\n", #Terms(h3(Ec2));
    end if;
  end if;
else
  print "E0 not linear in d; printing E0 for manual handling";
  Write("data/claude_v7_stratum2_eqs.txt", Sprintf("E0 := %o;\n", G0) : Overwrite := true);
end if;
print "STRATUM2_DONE";
quit;
