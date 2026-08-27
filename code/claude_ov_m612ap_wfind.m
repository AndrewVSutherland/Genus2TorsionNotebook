//////////////////////////////////////////////////////////////////////
// claude_ov_m612ap_wfind.m    lane 9 ([6,12]) 2026-07-25
//
// Find an INTRINSIC element of the Prym Pr = (1-iota) J(E8) over Q,
// with no reliance on the (uncertified) bigonal transport.
//
//   E8 : Q8(x,y) = y^8 + A(x) y^4 + B(x) y^2 + C(x) = 0,   genus 4
//   iota : y -> -y,   E8/iota = E4 (genus 2),  Fix(iota) = the two
//   rational boundary places over x = 0.
//
// Q8 is EVEN in y, so Q8(x0,y) = h(x0, y^2) with
//   h(x0,z) = z^4 + A z^2 + B z + C.
// Over Q the irreducible factors of Q8(x0,y) either (a) are even in y --
// then the corresponding place is iota-STABLE and contributes 0 to the
// Prym -- or (b) come in pairs g(y), g(-y) swapped by iota.  Case (b) is
// exactly what we want: for such a factor g the place P_g satisfies
//   W := [P_g - iota P_g] = [P_g - P_{g(-y)}]  in  (1-iota) J(E8)(Q),
// a genuine rational Prym class, given by an explicit rational divisor
// so that its reduction mod any good prime is computable.
//
// This script scans x0 = a/b of bounded height and reports every
// non-even irreducible factor, ranked by degree and coefficient size.
// It ALSO reports rational points (linear factors y - c with c != 0),
// which would be an immediate DISCOVERY (a non-boundary rational point).
//
// Usage: magma -b H:=40 code/claude_ov_m612ap_wfind.m
//////////////////////////////////////////////////////////////////////
SetColumns(0);
SetMemoryLimit(4*10^9);
if not assigned H then H := 40; elif Type(H) eq MonStgElt then H := StringToInteger(H); end if;

QQ := Rationals();
Pxy<x,y> := PolynomialRing(QQ, 2);
A := 216*x^4+72*x^3-24*x^2;
B := -1296*x^6-1728*x^5-432*x^4+64*x^3;
C := -3888*x^8-2592*x^7+432*x^6+288*x^5-48*x^4;
Q8 := y^8 + A*y^4 + B*y^2 + C;

Py<t> := PolynomialRing(QQ);

printf "=== claude_ov_m612ap_wfind : intrinsic Prym class search, H = %o ===\n", H;

// helper: is a univariate polynomial even in t?
IsEvenPoly := function(g)
  return &and[ Coefficient(g,i) eq 0 : i in [1..Degree(g) by 2] ];
end function;

xs := [];
for b in [1..H] do
  for a in [-H*b .. H*b] do
    if GCD(Abs(a), b) ne 1 then continue; end if;
    if a eq 0 then continue; end if;
    if Abs(a) gt H*b then continue; end if;
    Append(~xs, QQ!a/b);
  end for;
end for;
xs := Sort(Setseq(Seqset(xs)), func<u,v | Max(Abs(Numerator(u)),Denominator(u)) - Max(Abs(Numerator(v)),Denominator(v))>);
printf "scanning %o rational x-values\n", #xs;

nrat := 0; nodd := 0; cnt := 0;
best := [];
for x0 in xs do
  cnt +:= 1;
  Av := Evaluate(A, [x0,0]); Bv := Evaluate(B, [x0,0]); Cv := Evaluate(C, [x0,0]);
  g := t^8 + Av*t^4 + Bv*t^2 + Cv;
  if Cv eq 0 then continue; end if;          // y | Q8, degenerate fibre
  if Discriminant(g) eq 0 then continue; end if;  // ramified fibre, skip
  fac := Factorization(g);
  for fe in fac do
    gg := fe[1];
    if IsEvenPoly(gg) then continue; end if;
    d := Degree(gg);
    // size measure
    hh := Max([ Max(Abs(Numerator(c)), Denominator(c)) : c in Coefficients(gg) ]);
    if d eq 1 then
      nrat +:= 1;
      printf "*** RATIONAL POINT on E8: x = %o, y = %o\n", x0, -Coefficient(gg,0);
    end if;
    nodd +:= 1;
    Append(~best, <d, hh, x0, gg>);
  end for;
  if cnt mod 500 eq 0 then printf "PROGRESS %o/%o  found %o\n", cnt, #xs, #best; end if;
end for;

printf "non-even factors found: %o (rational points: %o)\n", #best, nrat;
best := Sort(best, func<u,v | u[1] ne v[1] select u[1]-v[1] else (u[2] lt v[2] select -1 else (u[2] gt v[2] select 1 else 0))>);
for i in [1..Min(40, #best)] do
  b := best[i];
  printf "CAND deg %o  x0 = %o  g = %o\n", b[1], b[3], b[4];
end for;
print "WFIND_DONE";
quit;
