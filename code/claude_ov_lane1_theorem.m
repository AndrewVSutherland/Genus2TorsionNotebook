// claude_ov_lane1_theorem.m -- Lane 1, 2026-07-25 (resumed session).
//
// THE FAMILY THEOREM for the u = -1/2 contact-7 [2,2,14] family.
//
// Setting.  E : y^2 = x^3+6x^2-16 (conductor 288, disc 2^12*3^3),
// E(Q) = <G=(4,-12)> x <T=(-2,0)> = Z x Z/2.  For P in E(Q) put
//    m = 4/x(P),  w = y(P)*m^2/4     (so w^2 = -m^4+6m^2+4m),
//    s = (m^2+w)/d, t = (m^2-w)/d, u = -1/2,  d = (m+1)^2(m-2),
//    c4 = (G7(s)-G7(t))/(s^2-t^2), c0 = G7(s)-c4*s^2, b = c4-2, a = 9/2-c0-c4,
//    h = 1-(7/2)x+a x^2+b x^3,   C_P : y^2 = f(x) := (h(x)^2+(x-1)^7)/x^2,
// with G7(v) = -(v^5-v^3-v^2/2)/(v+1)^2.
//
// Unconditional containment (no computation needed, valid for EVERY member):
//   * f(1-s^2) = f(1-t^2) = f(1-u^2) = 0, three distinct rational roots, so
//     f has >= 4 irreducible factors and the rational 2-rank is >= 3;
//   * on C_P, (x*y-h)(x*y+h) = x^2 f - h^2 = (x-1)^7 and x*y-h vanishes only
//     at Pt=(1,h(1)) (order 7) with a pole of order 7 at infinity, so
//     7*[Pt-infty] = 0; the class is nonzero (a genus-2 curve carries no
//     degree-1 rational function), hence has order exactly 7.
//   => J(Q)_tors contains (Z/2)^3 x Z/7 = [2,2,14].
//   * 2-rank is EXACTLY 3: the residual quadratic splits iff g(-m) is a
//     square (results/claude_ov_lane1_splitall.log), and the genus-2
//     Chabauty computation in code/claude_ov_lane1_order112.m shows that
//     happens only at the degenerate m in {0,+-1,+-2}.
//
// Exactness + genericity on an infinite subfamily.  Everything above is a
// rational function of P, so the reduction of C_P mod an odd prime p of good
// reduction for E depends ONLY on the image of P in E(F_p).  Hence for a
// finite set S of primes and N := lcm_{p in S} ord(Gbar in E(F_p)), all the
// members in the coset P_0 + N*<G> share their #J(F_p) and Frobenius charpoly
// chi_p for every p in S.  Choosing S so that
//     gcd_{p in S} #J(F_p) = 56                (torsion injects: tors | 56)
//     and S carries a strict disjoint certificate pair (p0,q0)
// proves that the WHOLE coset -- an infinite set -- consists of curves with
// J(Q)_tors exactly [2,2,14] and End(J_Qbar) = Z.
//
// This script searches base members and prime sets for the smallest N,
// verifies that the construction really descends to E(F_p), and verifies the
// congruence prediction at coset members by direct computation.
//
// Usage: magma -b KSCAN:=12 PMAX:=500 code/claude_ov_lane1_theorem.m

SetColumns(0);
if not assigned MemGB then MemGB := 8;   elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
if not assigned KSCAN then KSCAN := 12;  elif Type(KSCAN) eq MonStgElt then KSCAN := StringToInteger(KSCAN); end if;
if not assigned PMAX  then PMAX  := 500; elif Type(PMAX)  eq MonStgElt then PMAX  := StringToInteger(PMAX);  end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals(); Z := Integers();
P<x> := PolynomialRing(Q);
E := EllipticCurve([0,6,0,0,-16]);
Gp := E![4,-12]; Tp := E![-2,0];
dE := Z!Discriminant(E);
printf "E = %o ; disc = %o = %o ; conductor %o ; E(Q) = Z x Z/2\n", E, dE, Factorization(dE), Conductor(E);

function BuildF(RG, mm, ww)      // the construction over any field, char <> 2
  FF := BaseRing(RG); xx := RG.1;
  d := (mm+1)^2*(mm-2);
  if d eq 0 then return false,_,_,_,_,_; end if;
  s0 := (mm^2+ww)/d; t0 := (mm^2-ww)/d; u0 := FF!(-1/2);
  vs := [s0,t0,u0];
  if #{v : v in vs} ne 3 or #{v^2 : v in vs} ne 3 then return false,_,_,_,_,_; end if;
  if (FF!(-1) in vs) or (FF!0 in vs) then return false,_,_,_,_,_; end if;
  G7 := func<u | -(u^5 - u^3 - u^2/2)/(u+1)^2>;
  c4 := (G7(s0)-G7(t0))/(s0^2-t0^2); c0 := G7(s0) - c4*s0^2;
  b := c4 - 2; a := 9/2 - c0 - c4;
  h := 1 - (FF!7/2)*xx + a*xx^2 + b*xx^3;
  num := h^2 + (xx-1)^7;
  if num mod xx^2 ne 0 then return false,_,_,_,_,_; end if;
  f := num div xx^2;
  if Degree(f) ne 5 then return false,_,_,_,_,_; end if;
  return true, f, h, s0, t0, u0;
end function;

function MemberQ(k, eps)
  Pt := k*Gp + eps*Tp;
  if Pt eq E!0 then return false,_,_,_,_,_,_; end if;
  xc := Pt[1]/Pt[3];
  if xc eq 0 then return false,_,_,_,_,_,_; end if;
  mm := 4/xc; ww := (Pt[2]/Pt[3])*mm^2/4;
  ok, f, h, s0, t0, u0 := BuildF(P, mm, ww);
  if not ok then return false,_,_,_,_,_,_; end if;
  return true, f, h, mm, s0, t0, u0;
end function;

function PrimeData(f, p)
  Fp := GF(p); RP<xp> := PolynomialRing(Fp);
  cs := Coefficients(f);
  for c in cs do if Valuation(Denominator(c), p) gt 0 then return false,_,_,_; end if; end for;
  fp := RP![Fp!c : c in cs];
  if Degree(fp) ne 5 or not IsSquarefree(fp) then return false,_,_,_; end if;
  Lp := EulerFactor(Jacobian(HyperellipticCurve(fp)));
  Np := Z!Evaluate(Lp, 1);
  chi := P!Reverse(Coefficients(Lp));
  certok := false;
  if IsIrreducible(chi) then
    K := NumberField(chi); pi := K.1; certok := true;
    for nn in [2..12] do
      if Degree(MinimalPolynomial(pi^nn)) lt 4 then certok := false; break; end if;
    end for;
  end if;
  return true, Np, chi, certok;
end function;

ordG := AssociativeArray();
plist := [];
for p in PrimesInInterval(5, PMAX) do
  if dE mod p eq 0 then continue; end if;
  ordG[p] := Order(ChangeRing(E,GF(p))!Gp);
  Append(~plist, p);
end for;

// ---- scan base members for the smallest congruence modulus N ----
gbest := <10^9, 0, 0, [], 0, 0, 0, 0>;
// KSCAN := 0 short-circuits the (slow) base-member scan and installs the
// minimum it already found: base k0=3, eps=0, S = {43,263,307}, N = 11.
if KSCAN eq 0 then gbest := <11, 3, 0, [43,263,307], 43, 263, 43, 307>; end if;
for k0 in [1..KSCAN] do
 for eps0 in [0,1] do
  ok, fb := MemberQ(k0, eps0);
  if not ok then continue; end if;
  pool := [];
  for p in plist do
    okp, Np, chi, certok := PrimeData(fb, p);
    if okp then Append(~pool, <p, ordG[p], Np, chi, certok>); end if;
  end for;
  sfd := AssociativeArray();
  for e in pool do if e[5] then sfd[e[1]] := SplittingField(e[4]); end if; end for;
  best := <10^9,0,0,0,0>;
  for i in [1..#pool] do
   for j in [i..#pool] do
    if GCD(pool[i][3], pool[j][3]) ne 56 then continue; end if;
    for a in [1..#pool] do
     if not pool[a][5] then continue; end if;
     for bb in [a+1..#pool] do
      if not pool[bb][5] then continue; end if;
      SS := {pool[i][1], pool[j][1], pool[a][1], pool[bb][1]};
      NN := 1; for p in SS do NN := LCM(NN, ordG[p]); end for;
      if NN ge best[1] then continue; end if;
      d0 := Degree(sfd[pool[a][1]]); d1 := Degree(sfd[pool[bb][1]]);
      if Degree(SplittingField(pool[a][4]*pool[bb][4])) ne d0*d1 then continue; end if;
      best := <NN, pool[i][1], pool[j][1], pool[a][1], pool[bb][1]>;
     end for;
    end for;
   end for;
  end for;
  if best[1] lt 10^9 then
    SS := Sort(Setseq({best[2],best[3],best[4],best[5]}));
    printf "base k0=%-3o eps=%o : best N = %-6o  torsion witnesses (%o,%o)  cert (%o,%o)  S=%o\n",
       k0, eps0, best[1], best[2], best[3], best[4], best[5], SS;
    if best[1] lt gbest[1] then
      gbest := <best[1], k0, eps0, SS, best[2], best[3], best[4], best[5]>;
    end if;
  end if;
 end for;
end for;

NN := gbest[1]; K0 := gbest[2]; EPS0 := gbest[3]; S := gbest[4];
printf "\n=== BEST CONGRUENCE FAMILY ===\n";
printf "  base member k0 = %o, eps = %o ; prime set S = %o ; N = %o\n", K0, EPS0, S, NN;
printf "  torsion witnesses (p1,p2) = (%o,%o) ; End=Z certificate (p0,q0) = (%o,%o)\n",
   gbest[5], gbest[6], gbest[7], gbest[8];

okb, fb, hb, mb, sb, tb, ub := MemberQ(K0, EPS0);
printf "  m0 = %o ; (s,t,u) = (%o,%o,%o)\n", mb, sb, tb, ub;
printf "  IDENTITY x^2*f - h^2 = (x-1)^7 : %o ; h(1) = %o\n", x^2*fb - hb^2 eq (x-1)^7, Evaluate(hb,1);
printf "  the three rational roots 1-v^2 are roots of f : %o\n",
   &and[Evaluate(fb, 1-v^2) eq 0 : v in [sb,tb,ub]];
for p in S do
  okp, Np, chi, certok := PrimeData(fb, p);
  printf "  p=%-4o ordGbar=%-4o #J(F_p)=%-8o cert=%-5o chi=%o\n", p, ordG[p], Np, certok, chi;
end for;
gS := 0;
for p in S do okp, Np := PrimeData(fb, p); gS := GCD(gS, Np); end for;
printf "  gcd over S = %o\n", gS;

printf "\n=== the construction descends to E(F_p) (rebuild from Pbar alone) ===\n";
Pt0 := K0*Gp + EPS0*Tp;
for p in S do
  Fp := GF(p); RP<xp> := PolynomialRing(Fp);
  Rb := ChangeRing(E, Fp)!Pt0;
  mbar := 4/(Rb[1]/Rb[3]); wbar := (Rb[2]/Rb[3])*mbar^2/4;
  okp, fpc := BuildF(RP, mbar, wbar);
  fred := RP![Fp!c : c in Coefficients(fb)];
  printf "  p=%-4o Pbar=%o : built-from-Pbar equals reduction of f : %o\n", p, Rb, okp and (fpc eq fred);
end for;

printf "\n=== direct verification of the congruence prediction ===\n";
baseD := [];
for p in S do okp, Np, chi := PrimeData(fb, p); Append(~baseD, <p, Np, chi>); end for;
for n in [-2,-1,1,2] do
  kk := K0 + n*NN;
  t0 := Cputime();
  ok, ff, hh, mm, s1, t1, u1 := MemberQ(kk, EPS0);
  if not ok then printf "  n=%-3o k=%-6o : DEGENERATE (unexpected)\n", n, kk; continue; end if;
  agree := true; gg := 0;
  for tri in baseD do
    okp, Np, chi := PrimeData(ff, tri[1]);
    if (not okp) or Np ne tri[2] or chi ne tri[3] then agree := false; end if;
    if okp then gg := GCD(gg, Np); end if;
  end for;
  roots3 := &and[Evaluate(ff, 1-v^2) eq 0 : v in [s1,t1,u1]] and #{s1,t1,u1} eq 3;
  gneg := -mm^4 + 6*mm^2 - 4*mm;                       // 2-rank 4 iff this is a square
  hgt := Maximum([Abs(Numerator(cc)) + Abs(Denominator(cc)) : cc in Coefficients(ff)]);
  printf "  n=%-3o k=%-6o log10ht=%-7o 3 rational roots: %o ; g(-m) square (2-rank 4): %o ; same (N_p,chi_p) on S: %o ; gcd=%o => tors EXACTLY [2,2,14] and End=Z: %o  [%os]\n",
     n, kk, Ilog(10,Maximum(2,hgt)), roots3, IsSquare(gneg), agree, gg,
     (roots3 and (not IsSquare(gneg)) and agree and gg eq 56), Cputime(t0);
end for;

printf "\n=== absolute Igusa invariant j1 (Qbar-isomorphism invariant), small members ===\n";
igs := [];
for kk in [2..8] do
  ok, ff := MemberQ(kk, 1);
  if not ok then continue; end if;
  // NOTE: this Magma has no AbsoluteIgusaInvariants; G2Invariants returns the
  // three absolute (Qbar-isomorphism) invariants of a genus-2 curve.
  iv := G2Invariants(HyperellipticCurve(ff));
  Append(~igs, iv[1]);
  printf "  k=%o eps=1  j1 = %o\n", kk, iv[1];
end for;
printf "  pairwise distinct: %o  => m |-> C_P is finite-to-one, so the coset gives infinitely many Qbar-isomorphism classes\n", #Seqset(igs) eq #igs;

printf "\nLANE1_THEOREM_DONE\n";
quit;
