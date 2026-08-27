// claude_ov_lane1_cover.m -- Lane 1, 2026-07-25 (third session).
//
// UPGRADE of Theorem C: from "ONE coset k = 3 mod 11" to "EVERY member except
// an explicit handful of congruence classes".
//
// Setting (identical to code/claude_ov_lane1_theorem.m).
//   E : y^2 = x^3+6x^2-16 (conductor 288), E(Q) = <G=(4,-12)> x <T=(-2,0)>.
//   For P = k*G + eps*T put m = 4/x(P), w = y(P)*m^2/4, d = (m+1)^2(m-2),
//   s = (m^2+w)/d, t = (m^2-w)/d, u = -1/2, G7(v) = -(v^5-v^3-v^2/2)/(v+1)^2,
//   c4 = (G7(s)-G7(t))/(s^2-t^2), c0 = G7(s)-c4 s^2, b = c4-2, a = 9/2-c0-c4,
//   h = 1-(7/2)x+a x^2+b x^3, C_P : y^2 = f_P(x) = (h^2+(x-1)^7)/x^2.
//
// KEY LEMMA.  Every ingredient is a rational function of P with Z[1/2]
// coefficients.  Let p >= 5 be of good reduction for E, Pbar the image of P in
// E(F_p).  If the SAME construction run over F_p starting from Pbar succeeds
// (Pbar != O, x(Pbar) != 0, d != 0, s,t,u pairwise distinct with pairwise
// distinct squares, none of them 0 or -1, f_bar squarefree of degree 5), then
// every denominator is a p-unit at P, f_P is p-integral, f_P = f_bar mod p, and
// C_P has GOOD reduction at p.  Hence (p odd)
//      J_P(Q)_tors  ->  J_{f_bar}(F_p)   is INJECTIVE,
// and the Frobenius charpoly of J_P mod p is chi_{f_bar}.  Both depend only on
// Pbar, i.e. only on (k mod ord(Gbar_p), eps).   [verified numerically below]
//
// COVERING ARGUMENT.  Fix L and take S = all good p <= PMAX with n_p :=
// ord(Gbar_p) dividing L.  The class (k mod L, eps) then determines Pbar at
// EVERY p in S simultaneously.  Write T = J_P(Q)_tors and g = gcd over the
// usable p of #J(F_p).  By Theorem A, T contains (Z/2)^3 x Z/7 and has 2-rank
// EXACTLY 3, so T_2 = Z/2^a x Z/2^b x Z/2^c with 1 <= a <= b <= c.  A class is
// CERTIFIED (T = [2,2,14] exactly, End(J_Qbar) = Z) when
//   (i)  odd(g) = 7                                     -> T_odd = Z/7 ;
//   (ii) v_2(g) = 3  OR  some usable p has J(F_p)_2 ELEMENTARY ABELIAN, i.e.
//        v_2(#J(F_p)) = 2-rank J(F_p) = (#irreducible factors of f_bar) - 1.
//        In the second case T_2 has exponent 2, hence T_2 = (Z/2)^3.  This is
//        strictly stronger than demanding g = 56 and it is free: the factor
//        count is a by-product of building f_bar.       -> T_2 = (Z/2)^3 ;
//   (iii) two usable p carry strict Frobenius quartics with linearly disjoint
//        splitting fields                               -> End(J_Qbar) = Z.
// The classes that no prime can ever see are exactly the congruence classes of
// the eight GLOBALLY degenerate points {O, +-G, +-G+T, +-2G, T}.
//
// Note: for a FIXED n, only finitely many p have ord(Gbar_p) = n (they divide
// the n-th elliptic divisibility term), which is why we must combine many n | L
// rather than take many primes of one order.  Cost is controlled by tabulating
// (#J(F_p), chi_p) once per (p, r mod n_p, eps) -- 2*sum(n_p) Euler factors --
// and then doing the 2L classes by lookup.
//
// Usage:
//   code/claude_magma_slot.sh -b PMAX:=2500 LTARGET:=3960 MemGB:=6 \
//        code/claude_ov_lane1_cover.m > results/claude_ov_lane1_cover.log 2>&1 &

SetColumns(0);
if not assigned MemGB   then MemGB   := 6;    elif Type(MemGB)   eq MonStgElt then MemGB   := StringToInteger(MemGB);   end if;
if not assigned PMAX    then PMAX    := 2500; elif Type(PMAX)    eq MonStgElt then PMAX    := StringToInteger(PMAX);    end if;
if not assigned LTARGET then LTARGET := 3960; elif Type(LTARGET) eq MonStgElt then LTARGET := StringToInteger(LTARGET); end if;
if not assigned NMAX    then NMAX    := 120;  elif Type(NMAX)    eq MonStgElt then NMAX    := StringToInteger(NMAX);    end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals(); Z := Integers();
P<x> := PolynomialRing(Q);
E := EllipticCurve([0,6,0,0,-16]);
Gp := E![4,-12]; Tp := E![-2,0];
dE := Z!Discriminant(E);
printf "E = %o ; disc = %o ; conductor %o ; E(Q) = <G=(4,-12)> x <T=(-2,0)> = Z x Z/2\n",
   E, Factorization(dE), Conductor(E);
printf "PMAX = %o ; LTARGET = %o\n", PMAX, LTARGET;

// ---------------------------------------------------------------- construction
function BuildF(RG, mm, ww)      // over Q and over F_p, char <> 2
  FF := BaseRing(RG); xx := RG.1;
  d := (mm+1)^2*(mm-2);
  if d eq 0 then return false,_; end if;
  s0 := (mm^2+ww)/d; t0 := (mm^2-ww)/d; u0 := FF!(-1/2);
  vs := [s0,t0,u0];
  if #{v : v in vs} ne 3 or #{v^2 : v in vs} ne 3 then return false,_; end if;
  if (FF!(-1) in vs) or (FF!0 in vs) then return false,_; end if;
  G7 := func<v | -(v^5 - v^3 - v^2/2)/(v+1)^2>;
  c4 := (G7(s0)-G7(t0))/(s0^2-t0^2); c0 := G7(s0) - c4*s0^2;
  b := c4 - 2; a := 9/2 - c0 - c4;
  h := 1 - (FF!7/2)*xx + a*xx^2 + b*xx^3;
  num := h^2 + (xx-1)^7;
  if num mod xx^2 ne 0 then return false,_; end if;
  f := num div xx^2;
  if Degree(f) ne 5 or not IsSquarefree(f) then return false,_; end if;
  return true, f;
end function;

function MemberQ(k, eps)
  Pt := k*Gp + eps*Tp;
  if Pt eq E!0 then return false,_; end if;
  xc := Pt[1]/Pt[3];
  if xc eq 0 then return false,_; end if;
  mm := 4/xc; ww := (Pt[2]/Pt[3])*mm^2/4;
  return BuildF(P, mm, ww);
end function;

// ------------------------------------------------- the globally degenerate P's
printf "\n=== globally degenerate members (no genus-2 curve) ===\n";
degcl := [];
for k in [-6..6] do
  for eps in [0,1] do
    Pt := k*Gp + eps*Tp;
    if Pt eq E!0 then Append(~degcl, <k,eps,"O">); continue; end if;
    ok := MemberQ(k, eps);
    if not ok then Append(~degcl, <k,eps,Sprint(Pt[1]/Pt[3])>); end if;
  end for;
end for;
printf "  (k,eps,x(P)) with |k|<=6 : %o\n", degcl;
printf "  => the degenerate set is D = {O, +-G, +-G+T, +-2G, T}, i.e. k in {0,+-1,+-2}\n";

// ---------------------------------------------------------------- prime pool
printf "\n=== prime pool: good p <= %o with n_p = ord(Gbar_p) dividing L = %o ===\n", PMAX, LTARGET;
S := []; np := AssociativeArray();
for p in PrimesInInterval(5, PMAX) do
  if dE mod p eq 0 then continue; end if;
  n := Order(ChangeRing(E,GF(p))!Gp);
  if LTARGET mod n eq 0 and n le NMAX then Append(~S, p); np[p] := n; end if;
end for;
L := 1; for p in S do L := LCM(L, np[p]); end for;
printf "  %o primes ; lcm of their n_p is L = %o\n", #S, L;
printf "  (p, n_p) = %o\n", [<p, np[p]> : p in S];
printf "  total Euler factors to tabulate = 2*sum n_p = %o\n", 2*&+[np[p] : p in S];

// -------------------------------------------------- tabulate (N_p, chi_p, strict)
// tab[p] is a sequence indexed by (2*r + eps + 1), r = 0..n_p-1
t0 := Cputime();
tab := AssociativeArray();
sfdeg := AssociativeArray();     // cache: SplittingField degree, keyed by chi
for p in S do
  Fp := GF(p); Ep := ChangeRing(E, Fp); Gb := Ep!Gp; Tb := Ep!Tp;
  RP := PolynomialRing(Fp);
  row := [];
  for r in [0..np[p]-1] do
    for eps in [0,1] do
      Pb := r*Gb + eps*Tb;
      ok := false; Np := 0; chi := P!0; strict := false; r2 := 0;
      if Pb ne Ep!0 and Pb[3] ne 0 then
        xc := Pb[1]/Pb[3];
        if xc ne 0 then
          mm := 4/xc; ww := (Pb[2]/Pb[3])*mm^2/4;
          okb, fp := BuildF(RP, mm, ww);
          if okb then
            ok := true;
            // 2-rank of J(F_p): deg f = 5, so it is (#irreducible factors) - 1
            r2 := #Factorization(fp) - 1;
            Lp := EulerFactor(Jacobian(HyperellipticCurve(fp)));
            Np := Z!Evaluate(Lp, 1);
            chi := P!Reverse(Coefficients(Lp));
            if IsIrreducible(chi) then
              K := NumberField(chi); pi := K.1; strict := true;
              for nn in [2..12] do
                if Degree(MinimalPolynomial(pi^nn)) lt 4 then strict := false; break; end if;
              end for;
            end if;
            if strict and not IsDefined(sfdeg, chi) then sfdeg[chi] := Degree(SplittingField(chi)); end if;
          end if;
        end if;
      end if;
      Append(~row, <ok, Np, chi, strict, r2>);
    end for;
  end for;
  tab[p] := row;
end for;
printf "  tabulated in %os ; %o distinct strict chi cached\n", Cputime(t0), #Keys(sfdeg);

// ------------------------------------------------------------ class sweep
// NOTE (Magma scoping): a `function` may NOT assign to a variable of the
// enclosing scope ("Imported environment value ... cannot be used as a local"),
// so the memo table for the linear-disjointness test is maintained inline in
// the sweep loop rather than inside a helper function.
disjcache := AssociativeArray();

printf "\n=== sweep of all %o classes (k mod %o, eps) ===\n", 2*L, L;
t0 := Cputime();
certified := {};  noprime := []; badgcd := []; nocert := []; nsaved := 0;
gcdhist := AssociativeArray();
for k in [0..L-1] do
 for eps in [0,1] do
  g := 0; strictps := []; elemab := false;
  for p in S do
    e := tab[p][2*(k mod np[p]) + eps + 1];
    if not e[1] then continue; end if;
    g := GCD(g, e[2]);
    // J(F_p)_2 elementary abelian  <=>  v_2(#J(F_p)) = 2-rank of J(F_p).
    // If that happens at ONE usable prime, J(Q)_tors_2 has exponent 2, so with
    // 2-rank exactly 3 (Theorem A) it IS (Z/2)^3 -- no Z/4, no order 112.
    if Valuation(e[2], 2) eq e[5] then elemab := true; end if;
    if e[4] then Append(~strictps, e[3]); end if;
  end for;
  if g eq 0 then Append(~noprime, <k,eps>); continue; end if;
  if not IsDefined(gcdhist, g) then gcdhist[g] := 0; end if;
  gcdhist[g] +:= 1;
  // odd part must be exactly 7 (T_odd contains Z/7 and divides odd(g));
  // 2-part must be (Z/2)^3, either because v_2(g) = 3 or by the elementary
  // abelian witness above.
  oddok := (g div 2^Valuation(g,2)) eq 7;
  twook := (Valuation(g,2) eq 3) or elemab;
  if not (oddok and twook) then Append(~badgcd, <k,eps,g>); continue; end if;
  if g ne 56 then nsaved +:= 1; end if;
  okc := false;
  for i in [1..Minimum(#strictps,8)] do
   for j in [i+1..Minimum(#strictps,8)] do
     c1 := strictps[i]; c2 := strictps[j];
     key := (c1 lt c2) select <c1,c2> else <c2,c1>;
     if IsDefined(disjcache, key) then
       dj := disjcache[key];
     else
       dj := Degree(SplittingField(c1*c2)) eq sfdeg[c1]*sfdeg[c2];
       disjcache[key] := dj;
     end if;
     if dj then okc := true; break; end if;
   end for;
   if okc then break; end if;
  end for;
  if okc then Include(~certified, <k,eps>); else Append(~nocert, <k,eps>); end if;
 end for;
end for;
printf "  swept in %os\n", Cputime(t0);
printf "  CERTIFIED (exact [2,2,14] and End=Z) : %o of %o classes\n", #certified, 2*L;
printf "  no usable prime at all               : %o classes  %o\n", #noprime, (#noprime le 600 select noprime else "(list suppressed)");
printf "  odd(g) <> 7 or 2-part not settled    : %o classes  %o\n", #badgcd, (#badgcd le 600 select badgcd else "(list suppressed)");
printf "  certified DESPITE gcd > 56, by the elementary-abelian 2-part test : %o classes\n", nsaved;
printf "  criteria (i),(ii) ok but no disjoint strict pair  : %o classes  %o\n", #nocert, (#nocert le 600 select nocert else "(list suppressed)");
printf "  gcd histogram: %o\n", [<g, gcdhist[g]> : g in Sort([gg : gg in Keys(gcdhist)])];

// the undecided classes -- are they exactly the degenerate ones?
undec := noprime cat [<e[1],e[2]> : e in badgcd] cat nocert;
degclasses := [];
for k in [-2..2] do
 for eps in [0,1] do
  ok := (k*Gp + eps*Tp eq E!0) select false else MemberQ(k, eps);
  if not ok then Append(~degclasses, <(k mod L), eps>); end if;
 end for;
end for;
printf "\ncongruence classes of the 8 degenerate points, mod L=%o : %o\n", L, Sort(degclasses);
printf "undecided classes = degenerate classes ? %o\n", Seqset(undec) eq Seqset(degclasses);
printf "undecided classes (k mod %o, eps) : %o\n", L, (#undec le 600 select Sort(undec) else "(suppressed)");

printf "\nDENSITY: %o of %o classes undecided = %o\n", #undec, 2*L, RealField(6)!(#undec/(2*L));

// smallest |k| that is NOT certified but IS a genuine member
printf "\n=== smallest genuine members not covered ===\n";
cnt := 0;
for k in [-40..40] do
 for eps in [0,1] do
  if <(k mod L), eps> in certified then continue; end if;
  if k*Gp + eps*Tp eq E!0 then continue; end if;
  ok := MemberQ(k, eps);
  if ok and cnt lt 20 then printf "   k=%-6o eps=%o : genuine member, NOT certified by this run\n", k, eps; end if;
  if ok then cnt +:= 1; end if;
 end for;
end for;
if cnt eq 0 then printf "   NONE with |k| <= 40: EVERY genuine member in that range is certified\n";
else printf "   (%o such members with |k| <= 40)\n", cnt; end if;

// representatives of the residual classes (NB: only the eight degenerate
// classes are guaranteed to contain no member of small height; the others do)
printf "\n=== representatives of the residual (undecided) classes ===\n";
for e in Sort(undec) do
  kk := e[1]; eps := e[2];
  reps := [kk - L, kk, kk + L];
  printf "   class (k = %o mod %o, eps = %o) : representatives %o\n", kk, L, eps, reps;
end for;

// ------------------------------------------- validation against the rational side
printf "\n=== validation: f_P mod p really equals the curve built from Pbar ===\n";
for tri in [<3,0>,<5,1>,<8,0>,<14,1>,<25,0>] do
  k := tri[1]; eps := tri[2];
  ok, fQ := MemberQ(k, eps);
  if not ok then continue; end if;
  agree := true; tested := 0; gg := 0;
  for p in S do
    Fp := GF(p); RP := PolynomialRing(Fp); Ep := ChangeRing(E, Fp);
    e := tab[p][2*(k mod np[p]) + eps + 1];
    if not e[1] then continue; end if;
    bad := false;
    for c in Coefficients(fQ) do if Valuation(Denominator(c), p) gt 0 then bad := true; end if; end for;
    if bad then agree := false; continue; end if;
    fred := RP![Fp!c : c in Coefficients(fQ)];
    if not IsSquarefree(fred) or Degree(fred) ne 5 then agree := false; continue; end if;
    Lp := EulerFactor(Jacobian(HyperellipticCurve(fred)));
    if Z!Evaluate(Lp,1) ne e[2] or P!Reverse(Coefficients(Lp)) ne e[3] then agree := false; end if;
    gg := GCD(gg, e[2]); tested +:= 1;
  end for;
  printf "  k=%-4o eps=%o : %o usable primes, (#J(F_p),chi_p) from the CLASS TABLE match the actual curve : %o ; gcd = %o\n",
     k, eps, tested, agree, gg;
end for;

printf "\nLANE1_COVER_DONE\n";
quit;
