//////////////////////////////////////////////////////////////////////
// claude_ov_m612prym_mwsieve.m   (lane 9, overnight 2026-07-25)
//
// THE ABEL-PRYM MORDELL-WEIL SIEVE for E8 (the [6,12] gatekeeper).
// Attack 2 of notes/contact6_m612_relative3_s3_quotients_2026_07_11.md,
// for which no code existed anywhere in the repository.
//
// SETUP.  E8 : Q8(x,y) = 0 is the genus-4 curve
//   y^8 + (216x^4+72x^3-24x^2)y^4 + (-1296x^6-1728x^5-432x^4+64x^3)y^2
//       + (-3888x^8-2592x^7+432x^6+288x^5-48x^4),
// (x = e, y = e*w, w^2 = Y, mp8(w,e) = mp4(w^2,e)).  It is EVEN in y, and
// the Prym involution is iota : y -> -y, with quotient E4 (genus 2).
//
//   rank J(E4)(Q) = 3, rank Prym(E8/E4)(Q) = 1, so rank J(E8)(Q) = 4 = g.
//
// THE MAP.  alpha : E8 -> Pr := (1-iota)J(E8),  P |-> [P - iota P].
// Both known rational points of E8 (the two boundary places over e = 0)
// are FIXED by iota, hence alpha(b1) = alpha(b2) = 0.  Conversely
// alpha(P) = 0 forces P = iota P (a degree-1 function would exist on a
// genus-4 curve otherwise), i.e. P is one of the finitely many
// ramification points of E8 -> E4.  So:
//
//   proving  alpha(E8(Q)) = {0}  proves  E8(Q) = Fix(iota)(Q),
//
// an explicit finite set.  Since Pr(Q) has rank 1, alpha(P) = n*G for an
// integer n (mod torsion), and the reduction of that ONE integer at
// several good primes is exactly what an MW sieve constrains.
//
// WHAT THIS SCRIPT COMPUTES, for one prime q:
//   1. the function field of E8 / F_q; asserts genus 4 (good reduction);
//   2. iota as a field automorphism, and its action on degree-1 places,
//      via a fingerprint of valuations of uniformisers (rigorous: the
//      fingerprint injectivity is ASSERTED);
//   3. Cl^0 = J(E8)(F_q) via ClassGroup;
//   4. Zbar := the reduction of the transported rank-1 Prym generator
//      Z = [D1 - Dinf] (the bigonal transport of the J(C2'_min)(Q)
//      generator (x^2+2x+4, 5x+5)), with an anti-invariance check
//      iota(Zbar) = -Zbar;
//   5. m := order of Zbar, and for every degree-1 place P the residue
//      n(P) in Z/m with alpha(P) = n(P)*Zbar  (or NOTINGP);
//   6. ell-divisibility of Zbar in J(E8)(F_q) for small ell -- a
//      NON-divisibility at any single q proves ell does not divide the
//      index [Pr(Q) : Z*Z + tors], which is what makes the sieve valid.
//
// Requires q = 1 mod 3 (so sqrt(-3) in F_q) and the bigonal discriminant
// to be a square mod q; otherwise the generator is not F_q-rational and
// the prime is skipped (reported as SKIP).
//
// Run:  magma -b q:=61 out:=/path/sieve_q61.txt claude_ov_m612prym_mwsieve.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned q then q := 61; else q := StringToInteger(q); end if;
if not assigned out then out := ""; end if;
if not assigned maxenum then maxenum := 20000000; else maxenum := StringToInteger(maxenum); end if;

lines := [];
rec := procedure(s)
  printf "%o\n", s;
end procedure;

printf "=== claude_ov_m612prym_mwsieve  q = %o ===\n", q;
if q in {2,3} then printf "SKIP %o small\n", q; quit; end if;
if q mod 3 ne 1 then printf "SKIP %o (q mod 3 = %o, need 1 so that sqrt(-3) in F_q)\n", q, q mod 3; quit; end if;

k := GF(q);
Fx<X> := RationalFunctionField(k);
Py<Y> := PolynomialRing(Fx);
Q8 := Y^8 + (216*X^4+72*X^3-24*X^2)*Y^4
    + (-1296*X^6-1728*X^5-432*X^4+64*X^3)*Y^2
    + (-3888*X^8-2592*X^7+432*X^6+288*X^5-48*X^4);
t0 := Cputime();
FF<yy> := FunctionField(Q8);
gg := Genus(FF);
printf "genus = %o (%o s)\n", gg, Cputime(t0);
if gg ne 4 then printf "SKIP %o BAD REDUCTION (genus %o)\n", q, gg; quit; end if;

iota := hom< FF -> FF | -yy >;
xf := FF!X;
assert iota(iota(yy)) eq yy;
assert iota(xf) eq xf;

pl1 := Places(FF, 1);
np := #pl1;
printf "#E8(F_%o) (degree-1 places) = %o\n", q, np;

// ---- iota on degree-1 places, via valuation fingerprints
unis := [UniformizingElement(P) : P in pl1];
rows := [ [ Valuation(f, P) : f in unis ] : P in pl1 ];
assert #Set(rows) eq np;             // fingerprint is injective: RIGOROUS
irows := [ [ Valuation(iota(f), P) : f in unis ] : P in pl1 ];
pos := AssociativeArray();
for i in [1..np] do pos[rows[i]] := i; end for;
iotaidx := [];
for i in [1..np] do
  assert IsDefined(pos, irows[i]);
  iotaidx[i] := pos[irows[i]];
end for;
assert [iotaidx[iotaidx[i]] : i in [1..np]] eq [1..np];   // involution
fixed := [i : i in [1..np] | iotaidx[i] eq i];
printf "iota-fixed degree-1 places: %o of %o (indices %o)\n", #fixed, np, fixed;

// ---- boundary places over x = 0
z0 := Zeros(Fx.1)[1];
dec := Decomposition(FF, z0);
bd := [];
for P in dec do
  if Degree(P) eq 1 then Append(~bd, P); end if;
end for;
printf "degree-1 places over x=0: %o with v(x) = %o\n",
    #bd, [Valuation(xf, P) : P in bd];

// ---- class group
t0 := Cputime();
Cl, mCl := ClassGroup(FF);
printf "ClassGroup invariants %o  (%o s)\n", Invariants(Cl), Cputime(t0);
h0 := &*[ i : i in Invariants(Cl) | i ne 0 ];
printf "#J(E8)(F_%o) = %o\n", q, h0;

// ---- the transported Prym generator, reduced mod q
// (exactly the char-0 transport used in the Coleman driver)
QQ := Rationals();
Kt<tt> := FunctionField(QQ);
xtf := -72/(tt^2+3); vtf := -72*tt/(tt^2+3);
utf := 6*(xtf^2-6*xtf-36)*vtf/(xtf^2*(xtf+24));
Atf := (xtf^4/4 - 15*xtf^3/2 + 324*xtf + 972)/(xtf^3*(xtf+24)/24);
Stf := 2*Atf + 2*utf;
numS := Numerator(Stf); denS := Denominator(Stf);
f2t := numS*denS;
f2sft := &*[ fe[1]^(fe[2] mod 2) : fe in Factorization(f2t) ] * LeadingCoefficient(f2t);
sq2 := f2t / f2sft;
assert Denominator(sq2) eq 1;
sqpoly := Sqrt(Numerator(sq2));

bad := false;
lamsq := k!(-3);
if not IsSquare(lamsq) then printf "SKIP %o : -3 not a square\n", q; quit; end if;
_, lam := IsSquare(lamsq);

D1 := [];      // <xval, yval> of the four transported generator points
for i in [1,2] do
  sgnl := (i eq 1) select 1 else -1;
  t0v := k!(-1) + sgnl*2*lam;
  s0m := (k!5)*(sgnl*lam)/3;
  if t0v^2 + 3 eq 0 then bad := true; break; end if;
  sqv  := &+[ k | k!Coefficient(sqpoly,j)*t0v^j : j in [0..Degree(sqpoly)] ];
  denv := &+[ k | k!Coefficient(denS,j)*t0v^j : j in [0..Degree(denS)] ];
  if denv eq 0 then bad := true; break; end if;
  s0 := s0m * sqv/denv;
  x0 := -72/(t0v^2+3);
  v0 := -72*t0v/(t0v^2+3);
  if x0 eq 0 or x0+24 eq 0 then bad := true; break; end if;
  u0 := 6*(x0^2-6*x0-36)*v0/(x0^2*(x0+24));
  if s0^2 - 2*((x0^4/4 - 15*x0^3/2 + 324*x0 + 972)/(x0^3*(x0+24)/24)) - 2*u0 ne 0 then
    printf "SKIP %o : transport consistency s^2 = 2A+2u FAILS mod q\n", q; quit;
  end if;
  dsc := s0^2 - 4*u0;
  if not IsSquare(dsc) then
    printf "SKIP %o : bigonal discriminant not a square mod q\n", q; quit;
  end if;
  _, sq := IsSquare(dsc);
  for sgn in [1,-1] do
    w0 := (s0 + sgn*sq)/2;
    if w0^2 - x0 eq 0 then bad := true; break; end if;
    e0 := (x0/3 + 2)/(w0^2 - x0);
    Append(~D1, <e0, e0*w0>);
  end for;
  if bad then break; end if;
end for;
if bad then printf "SKIP %o : degenerate transport mod q\n", q; quit; end if;

Dinf := [ <k!(-2)/3, 2*lam/3>, <k!(-2)/3, -2*lam/3> ];

// locate each (x,y) as a degree-1 place
findplace := function(xv, yv)
  zs := Zeros(Fx.1 - xv);
  if #zs eq 0 then return 0; end if;
  dc := Decomposition(FF, zs[1]);
  cand := [];
  for P in dc do
    if Degree(P) ne 1 then continue; end if;
    ok := true;
    try
      val := Evaluate(yy, P);
    catch err
      ok := false;
    end try;
    if ok and val eq yv then Append(~cand, P); end if;
  end for;
  if #cand ne 1 then return #cand eq 0 select 0 else -1; end if;
  return cand[1];
end function;

Dgen := DivisorGroup(FF)!0;
okall := true;
for pt in D1 do
  P := findplace(pt[1], pt[2]);
  if Type(P) eq RngIntElt then okall := false; break; end if;
  Dgen +:= Divisor(P);
end for;
for pt in Dinf do
  P := findplace(pt[1], pt[2]);
  if Type(P) eq RngIntElt then okall := false; break; end if;
  Dgen -:= Divisor(P);
end for;
if not okall then printf "SKIP %o : could not locate a transported place uniquely\n", q; quit; end if;
for P in bd do Dgen -:= Divisor(P); end for;
printf "generator divisor degree = %o (must be 0)\n", Degree(Dgen);
assert Degree(Dgen) eq 0;

Zcl := Dgen @@ mCl;
m := Order(Zcl);
printf "Zbar order m = %o  (h0/m = %o)\n", m, h0 div m;

// anti-invariance check: iota(Zbar) = -Zbar
DG := DivisorGroup(FF);
iotadiv := function(D)
  E := DG!0;
  for t in Support(D) do
    P := t;
    i := pos[[ Valuation(f, P) : f in unis ]];
    E +:= Valuation(D, P)*Divisor(pl1[iotaidx[i]]);
  end for;
  return E;
end function;
sup := Support(Dgen);
alldeg1 := &and[ Degree(P) eq 1 : P in sup ];
if alldeg1 then
  Ziot := iotadiv(Dgen) @@ mCl;
  printf "ANTIINVARIANCE iota(Zbar) + Zbar = %o (must be 0)\n", Ziot + Zcl;
else
  printf "ANTIINVARIANCE: support not all degree 1, skipped\n";
end if;

// ---- alpha(P) for every degree-1 place, and its residue n mod m
tab := AssociativeArray();
enum := (m le maxenum);
if enum then
  t0 := Cputime();
  cur := Cl!0;
  for n in [0..m-1] do
    tab[cur] := n;
    cur +:= Zcl;
  end for;
  printf "enumerated <Zbar> (%o elements, %o s)\n", m, Cputime(t0);
end if;

resid := [];
for i in [1..np] do
  P := pl1[i];
  D := Divisor(P) - Divisor(pl1[iotaidx[i]]);
  c := D @@ mCl;
  if enum then
    if IsDefined(tab, c) then
      Append(~resid, <i, tab[c]>);
    else
      Append(~resid, <i, -1>);
    end if;
  else
    Append(~resid, <i, -2>);
  end if;
end for;
inset := [r : r in resid | r[2] ge 0];
printf "ALPHA-RESIDUES q=%o m=%o : %o of %o places have alpha(P) in <Zbar>\n",
    q, m, #inset, np;
printf "RESIDUESET %o\n", Sort([r[2] : r in inset]);
for r in resid do
  P := pl1[r[1]];
  xv := "inf";
  try
    xv := Sprint(Evaluate(xf, P));
  catch err
    xv := "inf";
  end try;
  printf "PLACE %o x=%o iota=%o n=%o\n", r[1], xv, iotaidx[r[1]], r[2];
end for;

// ---- ell-divisibility of Zbar in J(E8)(F_q)  (saturation evidence)
gens := [Cl.i : i in [1..Ngens(Cl)]];
for ell in [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47] do
  if h0 mod ell ne 0 then continue; end if;
  H := sub< Cl | [ell*g : g in gens] >;
  if Zcl in H then
    printf "DIVTEST q=%o ell=%o : Zbar IS in ell*J(F_q) (no information)\n", q, ell;
  else
    printf "DIVTEST q=%o ell=%o : Zbar NOT in ell*J(F_q) => ell does NOT divide the index\n", q, ell;
  end if;
end for;

printf "MWSIEVE_DONE q=%o\n", q;
quit;
