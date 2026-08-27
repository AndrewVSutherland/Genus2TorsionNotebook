// claude_ov_lane7_rungs.m -- Lane 7 (overnight 2026-07-25).
//
// MEASURE-BEFORE-BUILD diagnostics for the two "side rungs" on the contact-7
// three-root surface, on all ELEVEN known three-root points.
//
// RUNG (a)  [2,2,28], order 112.  Reaching [2,2,28] from [2,2,14] = (Z/2)^3 x Z/7
// means halving one of the seven nonzero rational 2-torsion classes:
// (Z/2)^3 x Z/7 -> (Z/2)^2 x Z/4 x Z/7 = [2,2,28].  A NECESSARY condition for a
// class T in J(Q)[2] to be 2-divisible in J(Q) is that T mod p lies in 2*J(F_p)
// at EVERY good prime p (reduction is injective on prime-to-p torsion and
// commutes with multiplication).  We measure the per-class local pass rate over
// all good primes in [11,PMAX].  A class that fails at even one good prime can
// never be halved on THAT curve; a class that passes at every prime is a
// genuine candidate and would justify building the squareclass halving cover
// over the surface.  We also record, as a calibration, the same statistic for
// a random-model expectation (the observed rate over all 11 curves x 7 classes).
//
// RUNG (b)  order 168 = [2,2,42].  An INDEPENDENT rational 3-torsion class needs
// 3 | #J(F_p) at every good prime p != 3.  Measured directly.
//
// POSITIVE CONTROLS (rule 3): the same two tests are run on curves where the
// answer is known --
//   * control2 : y^2 = x(x+1)(x+2)... no; we use the [2,2,20] paper witness,
//     whose torsion (Z/2)^2 x Z/20 contains a rational point of order 4 whose
//     double is a rational 2-torsion class, so SOME class must pass at 100%;
//   * control3 : the [6,6] paper witness has full rational 3-torsion, so
//     3 | #J(F_p) must hold at 100% of good primes.
//
// Markers: RUNGA / RUNGB / CONTROL / LANE7_RUNGS_DONE
SetColumns(0);
SetMemoryLimit(3*10^9);
Q := Rationals(); P<x> := PolynomialRing(Q); Z := Integers();

if not assigned PMAX then PMAX := 500; elif Type(PMAX) eq MonStgElt then PMAX := StringToInteger(PMAX); end if;

Gfun := func< v | -(v^5 - v^3 - v^2/2)/(v+1)^2 >;

trips := [
  [-10, -10/7, -1/2], [-5, -15/8, -15/22], [-3, -3/4, -3/5],
  [-15/8, -15/19, -1/2], [-5/18, -10/49, 4/17], [-4/9, -4/25, 4/17],
  [-511/61, -511/625, -1/2], [-165/41, -33/16, -165/289],
  [-164/297, -1/2, 164/361], [-17/50, -34/189, 34/121], [-1/2, -13/49, 13/50]
];

// build the integral model y^2 = F from (s,t)
function BuildF(s,t)
  c4 := (Gfun(s)-Gfun(t))/(s^2-t^2);
  c0 := Gfun(s) - c4*s^2;
  b := c4 - 2; a := 9/2 - c0 - c4;
  h := 1 - (7/2)*x + a*x^2 + b*x^3;
  f := (h^2 + (x-1)^7) div x^2;
  den := LCM([Denominator(co) : co in Coefficients(f)]);
  F := P![ co*den^2 : co in Coefficients(f) ];
  return F, den;
end function;

// ---- generic tester: given integral F (deg 5 or 6) return per-class 2-divisibility
// pass counts over good primes, plus the 3 | #J(F_p) pass count.
procedure Diagnose(F, lbl)
  C := HyperellipticCurve(F);
  J := Jacobian(C);
  lc := Z!LeadingCoefficient(F);
  Dnum := Z!Numerator(Discriminant(F));
  // rational 2-torsion generators: from the rational roots of F
  rts := [ r[1] : r in Roots(F) ];
  gens := [];
  if Degree(F) eq 5 then
    // one rational Weierstrass point at infinity; e_i = [(r_i,0) - inf]
    for r in rts do Append(~gens, [* r *]); end for;
  else
    error "only deg-5 models handled here";
  end if;
  // nonzero elements of the F_2-span: nonempty subsets of the r_i
  subsets := [];
  n := #rts;
  for m in [1..2^n-1] do
    S := [ rts[i] : i in [1..n] | (m div 2^(i-1)) mod 2 eq 1 ];
    Append(~subsets, S);
  end for;
  ncl := #subsets;
  pass := [0 : i in [1..ncl]];
  tot := 0; pass3 := 0;
  badp := [ [] : i in [1..ncl] ];
  for p in PrimesInInterval(11, PMAX) do
    if lc mod p eq 0 or Dnum mod p eq 0 then continue; end if;
    Fp := PolynomialRing(GF(p)) ! F;
    if Degree(Fp) ne 5 or not IsSquarefree(Fp) then continue; end if;
    Cp := HyperellipticCurve(Fp);
    Jp := Jacobian(Cp);
    A, mp := AbelianGroup(Jp);
    twoA := sub< A | [2*g : g in Generators(A)] >;
    tot +:= 1;
    if #A mod 3 eq 0 then pass3 +:= 1; end if;
    PPp<xp> := PolynomialRing(GF(p));
    for i in [1..ncl] do
      up := &*[ xp - GF(p)!r : r in subsets[i] ];
      ok := true;
      try
        Dp := Jp ! [up, PPp!0];
      catch ee
        ok := false;
      end try;
      if not ok then continue; end if;   // should not happen at good p
      if (Dp @@ mp) in twoA then
        pass[i] +:= 1;
      else
        if #badp[i] lt 6 then Append(~badp[i], p); end if;
      end if;
    end for;
  end for;
  printf "RUNGA %o goodprimes=%o\n", lbl, tot;
  for i in [1..ncl] do
    printf "   class{%o} : 2-divisible at %o/%o primes (%o%%)  firstfail=%o\n",
      subsets[i], pass[i], tot, (tot eq 0 select 0 else (100*pass[i]) div tot), badp[i];
  end for;
  printf "RUNGB %o : 3 | #J(F_p) at %o/%o good primes (%o%%)\n",
    lbl, pass3, tot, (tot eq 0 select 0 else (100*pass3) div tot);
end procedure;

print "==== CONTROLS ====";
// control 3: [6,6] paper witness -- full rational 3-torsion, expect 100% on RUNGB
Fc3 := 11389248*x^5 - 18252000*x^4 + 42399396*x^3 - 10288044*x^2 + 29659500*x;
Diagnose(Fc3, "CONTROL66");

// control 2: a curve with a rational 4-torsion point doubling into J(Q)[2].
// y^2 = x(x+1)(x+2)(x+3)(x+4)-type models are split; use the [2,2,20] witness
// factored form of code/claude_endz_certificates.m, deg 6 -> take a deg-5
// alternative: the [2,2,2,8] witness x(x+1)(x+55^2)(x+99^2)(x+125^2) has
// torsion (Z/2)^3 x Z/8, so at least one 2-torsion class IS halvable.
Fc2 := x*(x+1)*(x+55^2)*(x+99^2)*(x+125^2);
Diagnose(Fc2, "CONTROL2228");

print "==== THE ELEVEN THREE-ROOT POINTS ====";
for T in trips do
  F, den := BuildF(T[1], T[2]);
  Diagnose(F, Sprintf("(%o,%o,%o)", T[1], T[2], T[3]));
end for;

print "LANE7_RUNGS_DONE";
quit;
