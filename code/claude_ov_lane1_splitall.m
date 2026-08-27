// claude_ov_lane1_splitall.m -- Lane 1, overnight 2026-07-25.
//
// RESTRICT THE ORDER-112 ("SPLITALL") CONDITION TO THE u=-1/2 [2,2,14] FAMILY.
//
// On the family (see code/claude_ov_lane1_family.m) the contact-7 quintic
//   Q(v) = v^5 + c4 v^4 + (2c4-1) v^3 + (c4+c0-1/2) v^2 + 2 c0 v + c0
// has the three rational roots  v = -1/2, s, t  with
//   s = (m^2+w)/D, t = (m^2-w)/D, D = (m+1)^2(m-2), w^2 = -m^4+6m^2+4m.
// Write Q(v) = (v+1/2)(v-s)(v-t)(v^2 + B v + C).  Comparing e1 and e5:
//   e1 = -c4  = (-1/2+s+t) - B      =>  B = (-1/2+s+t) + c4
//   e5 = -c0  = (-1/2)*s*t*C        =>  C = -c0/((-1/2)*s*t) = 2c0/(s t)
// The quintic splits completely (=> factor type [1,1,1,1,1] => 2-rank 4 =>
// torsion contains (Z/2)^3 x Z/7, i.e. [2,2,2,14] of ORDER 112, a new record)
// iff  Delta := B^2 - 4C  is a rational square.
//
// This script computes Delta as a function on the rank-1 elliptic curve
//   E : w^2 = -m^4+6m^2+4m   (equivalently y^2 = x^3+6x^2-16, conductor 288)
// and decides the resulting double cover: genus, and a rational-point search.
//
// Usage: magma -b code/claude_ov_lane1_splitall.m

SetColumns(0);
if not assigned MemGB then MemGB := 12; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
if not assigned KMAX then KMAX := 200; elif Type(KMAX) eq MonStgElt then KMAX := StringToInteger(KMAX); end if;
Q := Rationals(); Z := Integers();
Px<x> := PolynomialRing(Q);
Gfun := func<u0 | -(u0^5 - u0^3 - u0^2/2)/(u0+1)^2>;

/* ---- 0. numerical check of the residual-quadratic formulae on a member ---- */
E := EllipticCurve([0,6,0,0,-16]);
Gp := E![4,-12]; Tp := E![-2,0];
Params := function(k, eps)
  Pt := k*Gp + eps*Tp;
  if Pt eq E!0 then return false,_,_,_,_; end if;
  xc := Pt[1]/Pt[3]; yc := Pt[2]/Pt[3];
  if xc eq 0 then return false,_,_,_,_; end if;
  mm := 4/xc; ww := yc*mm^2/4;
  d := (mm+1)^2*(mm-2);
  if d eq 0 then return false,_,_,_,_; end if;
  s0 := (mm^2+ww)/d; t0 := (mm^2-ww)/d;
  return true, mm, ww, s0, t0;
end function;

printf "=== numerical check of B, C on family members ===\n";
for k in [-5..5] do for eps in [0,1] do
  okp, mm, ww, s0, t0 := Params(k,eps);
  if not okp then continue; end if;
  if #{s0,t0,Q!(-1/2)} ne 3 or #{s0^2,t0^2,Q!(1/4)} ne 3 then continue; end if;
  if s0 in {Q|-1,0} or t0 in {Q|-1,0} then continue; end if;
  c4 := (Gfun(s0)-Gfun(t0))/(s0^2-t0^2); c0 := Gfun(s0) - c4*s0^2;
  Qv := x^5 + c4*x^4 + (2*c4-1)*x^3 + (c4+c0-1/2)*x^2 + 2*c0*x + c0;
  ok1 := Evaluate(Qv,-1/2) eq 0 and Evaluate(Qv,s0) eq 0 and Evaluate(Qv,t0) eq 0;
  quo := Qv div ((x+1/2)*(x-s0)*(x-t0));
  B := (-1/2+s0+t0) + c4;  C := 2*c0/(s0*t0);
  ok2 := quo eq x^2 + B*x + C;
  Del := B^2-4*C;
  printf "k=%o eps=%o m=%o : Q has the 3 roots %o ; residual = x^2+Bx+C %o ; Delta=%o square=%o\n",
     k, eps, mm, ok1, ok2, Del, IsSquare(Del);
end for; end for;

/* ---- 1. Delta as a function on E, symbolically ---- */
Fm<m> := FunctionField(Q);
Pw<W> := PolynomialRing(Fm);
gm := -m^4+6*m^2+4*m;
Kw<wq> := ext<Fm | W^2 - gm>;
printf "\nfunction field Kw = Q(m)[w]/(w^2 - (%o)), genus %o\n", gm, Genus(Kw);
Dd := Kw!((m+1)^2*(m-2));
sS := (Kw!(m^2)+wq)/Dd; tS := (Kw!(m^2)-wq)/Dd;
GfunK := func<u0 | -(u0^5 - u0^3 - u0^2/2)/(u0+1)^2>;
c4S := (GfunK(sS)-GfunK(tS))/(sS^2-tS^2);
c0S := GfunK(sS) - c4S*sS^2;
BS := (Kw!(-1/2)+sS+tS) + c4S;
CS := 2*c0S/(sS*tS);
DelS := BS^2 - 4*CS;
printf "Delta (as an element of Kw):\n  %o\n", DelS;

// CLOSED FORM: Delta = 4*g(-m) / ((m-1)^4 (m+2)^2), g(m) = -m^4+6m^2+4m.
// Verified as an IDENTITY in Kw (no evaluation needed).
printf "\n=== closed form of Delta ===\n";
claim := Kw!( 4*(-m^4+6*m^2-4*m) / ((m-1)^4*(m+2)^2) );
printf "Delta == 4*g(-m)/((m-1)^4 (m+2)^2) identically in Kw : %o\n", DelS eq claim;
printf "  => Delta is a SQUARE  <=>  g(-m) = -m^4+6m^2-4m is a square\n";
printf "  => the order-112 locus on this family is the curve\n";
printf "       C : w^2 = g(m),  z^2 = g(-m)   (genus 5) -- see claude_ov_lane1_order112.m\n";

/* ---- 2. the double cover z^2 = Delta ---- */
Pz<Zv> := PolynomialRing(Kw);
try
  Kz := ext<Kw | Zv^2 - DelS>;
  printf "\nSPLITALL cover z^2 = Delta : genus = %o\n", Genus(Kz);
catch e
  printf "\ncould not build/measure the cover directly: %o\n", e`Object;
end try;

/* ---- 3. brute-force search along E(Q) ---- */
printf "\n=== search for SPLITALL along E(Q), |k| <= %o ===\n", KMAX;
found := 0;
for k in [-KMAX..KMAX] do for eps in [0,1] do
  okp, mm, ww, s0, t0 := Params(k,eps);
  if not okp then continue; end if;
  if #{s0,t0,Q!(-1/2)} ne 3 or #{s0^2,t0^2,Q!(1/4)} ne 3 then continue; end if;
  if s0 in {Q|-1,0} or t0 in {Q|-1,0} then continue; end if;
  c4 := (Gfun(s0)-Gfun(t0))/(s0^2-t0^2); c0 := Gfun(s0) - c4*s0^2;
  Bn := (-1/2+s0+t0) + c4; Cn := 2*c0/(s0*t0); Dn := Bn^2-4*Cn;
  if IsSquare(Dn) then
     printf "SPLITALL HIT k=%o eps=%o m=%o s=%o t=%o Delta=%o\n", k,eps,mm,s0,t0,Dn;
     found +:= 1;
  end if;
end for; end for;
printf "SPLITALL hits along E(Q) with |k| <= %o : %o\n", KMAX, found;
printf "LANE1_SPLITALL_DONE\n";
quit;
