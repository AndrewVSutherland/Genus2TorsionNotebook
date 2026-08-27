//////////////////////////////////////////////////////////////////////
// claude_ov_m612ap_sieve.m      lane 9 (overnight 2026-07-25)
//
// THE ABEL-PRYM MORDELL-WEIL SIEVE for the [6,12] gatekeeper curve
//
//     E8 :  Q8(x,y) = 0,   x = e,   y = e*w,   genus 4,
//     iota : y -> -y   (the Prym involution; E8/iota = E4, genus 2)
//
// Attack 2 of notes/contact6_m612_route.md.  No code for it has ever
// existed in this repo.  It is INDEPENDENT of the Coleman computation.
//
// THE ARGUMENT.
//   * Pr := (1-iota) J(E8) is the Prym, dim 2, rank Pr(Q) = 1 (proved,
//     code/contact6_m612_prym_rank_verifier.m).
//   * The Abel-Prym map  alpha : E8(Q) -> Pr(Q),  P |-> [P - iota P],
//     is defined over Q and kills exactly the iota-fixed points: if P is
//     not iota-fixed and [P - iota P] = 0 then E8 carries a degree-1 map
//     to P^1, absurd in genus 4.
//   * Pr(Q) = Z.G + T, T finite, so every rational point P has an integer
//     n(P) with alpha(P) = n(P).G (mod T).
//   * For a good prime q, reduction gives c_q(P) = n(P).Gbar_q, so
//     n(P) mod ord(Gbar_q) lies in
//        N_q := { n : n.Gbar_q in {[Pbar - iota Pbar] : Pbar in E8(F_q)} }.
//   * Intersecting the N_q over many q through a smooth common modulus
//     bounds n(P).  Combined with the Prym-Chabauty statement "at most one
//     rational point per residue disk at p = 37" this closes E8(Q).
//
// THE GENERATOR.  The class Z of the transported generator of
// J(C2'_min)(Q) = Z (Mumford (x^2+2x+4, 5x+5)), transported through the
// bigonal correspondence exactly as in Coleman stage 2:
//     t0 = 2*xm+1, xm = -1 +- sqrt(-3), s0m = (5*xm+5)/3,
//     s0 = s0m*sq(t0)/den(t0),  x0 = -72/(t0^2+3),  v0 = -72*t0/(t0^2+3),
//     u0 = 6(x0^2-6x0-36)v0/(x0^2(x0+24)),  w = (s0 +- sqrt(s0^2-4u0))/2,
//     X = (x0/3+2)/(w^2-x0),  Y = X*w.
// D1 = those four points; Dinf = the two transported infinities of C2'
// (X = -2/3, Y = +-2*lam/3) plus the two boundary places over x = 0.
// Z := D1 - Dinf, degree 0.  The script CHECKS iota(Z) = -Z in the class
// group (the Prym-membership certificate) rather than assuming it, and
// falls back to Z - iota(Z) if the check fails.
//
// Usage:  magma -b qmin:=7 qmax:=200 wdir:=... tag:=... claude_ov_m612ap_sieve.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned qmax then qmax := 200; else qmax := StringToInteger(qmax); end if;
if not assigned qmin then qmin := 7;   else qmin := StringToInteger(qmin); end if;
if not assigned wdir then wdir := ".";  end if;
if not assigned tag  then tag := "ap";  end if;
if not assigned nfork then nfork := 48; else nfork := StringToInteger(nfork); end if;
if not assigned elltop then elltop := 100; else elltop := StringToInteger(elltop); end if;

printf "=== claude_ov_m612ap_sieve : Abel-Prym MW sieve on E8 ===\n";
printf "primes in [%o, %o], fork width %o, wdir %o, tag %o\n",
    qmin, qmax, nfork, wdir, tag;

QQ := Rationals();

// ---- the transport auxiliaries, computed once over Q ----------------
Kt<tt> := FunctionField(QQ);
xtf := -72/(tt^2+3); vtf := -72*tt/(tt^2+3);
utf := 6*(xtf^2-6*xtf-36)*vtf/(xtf^2*(xtf+24));
Atf := (xtf^4/4 - 15*xtf^3/2 + 324*xtf + 972)/(xtf^3*(xtf+24)/24);
Stf := 2*Atf + 2*utf;
numS := Numerator(Stf); denS := Denominator(Stf);
f2t := numS*denS;
f2sft := &*[ fe[1]^(fe[2] mod 2) : fe in Factorization(f2t) ] * LeadingCoefficient(f2t);
sq2 := f2t / f2sft;
error if Denominator(sq2) ne 1, "sq2 is not a polynomial";
sqpoly := Sqrt(Numerator(sq2));
denpoly := denS;
printf "transport: deg den = %o, deg sq = %o\n", Degree(denpoly), Degree(sqpoly);
sqcoef := Coefficients(sqpoly);
dncoef := Coefficients(denpoly);

writelines := procedure(fn, L)
  fh := Open(fn, "w");
  for s in L do Puts(fh, s); end for;
  delete fh;
end procedure;

readlines := function(fn)
  ok := true; s := "";
  try s := Read(fn); catch err ok := false; end try;
  if not ok then return false, []; end if;
  return true, Split(s, "\n");
end function;

seqstr := function(S)
  if #S eq 0 then return ""; end if;
  return &cat[ IntegerToString(u) cat " " : u in S ];
end function;

// =====================================================================
// per-prime worker
// =====================================================================
worker := function(q, elltop, sqcoef, dncoef)
  L := [];
  Fq := GF(q);
  A<X> := RationalFunctionField(Fq);
  PY<Y> := PolynomialRing(A);
  Pq := PolynomialRing(Fq);
  Q8 := Y^8 + (216*X^4+72*X^3-24*X^2)*Y^4
        + (-1296*X^6-1728*X^5-432*X^4+64*X^3)*Y^2
        + (-3888*X^8-2592*X^7+432*X^6+288*X^5-48*X^4);
  if not IsIrreducible(Q8) then return ["SKIP reducible"]; end if;
  F<yy> := FunctionField(Q8);
  g := Genus(F);
  if g ne 4 then return [Sprintf("SKIP genus %o", g)]; end if;
  xx := F!(A.1);

  ok, lam := IsSquare(Fq!(-3));
  if not ok then return ["SKIP no sqrt(-3)"]; end if;
  sqp := Pq![Fq!c : c in sqcoef];
  dnp := Pq![Fq!c : c in dncoef];

  D1 := [];
  for i in [1,2] do
    xm := (i eq 1) select -1+lam else -1-lam;
    t0 := 2*xm+1; s0m := (5*xm+5)/3;
    dv := Evaluate(dnp, t0); sv := Evaluate(sqp, t0);
    if dv eq 0 or t0^2+3 eq 0 then return ["SKIP degenerate t0"]; end if;
    s0 := s0m*sv/dv;
    x0 := -72/(t0^2+3);
    v0 := -72*t0/(t0^2+3);
    if x0 eq 0 or x0 eq -24 then return ["SKIP degenerate x0"]; end if;
    u0 := 6*(x0^2-6*x0-36)*v0/(x0^2*(x0+24));
    dsc := s0^2 - 4*u0;
    bq, sq := IsSquare(dsc);
    if not bq then return ["SKIP dsc nonsquare"]; end if;
    for sgn in [1,-1] do
      w0 := (s0 + sgn*sq)/2;
      if w0^2 - x0 eq 0 then return ["SKIP degenerate w"]; end if;
      e0 := (x0/3 + 2)/(w0^2 - x0);
      Append(~D1, <e0, e0*w0>);
    end for;
  end for;
  Dinffin := [ <Fq!(-2)/3, 2*lam/3>, <Fq!(-2)/3, -2*lam/3> ];

  ev := function(X0, Y0)
    s := Fq!0;
    for i in [0..8] do
      ci := Coefficient(Q8, i);
      nv := Evaluate(Pq![Fq!c : c in Coefficients(Numerator(ci))], X0);
      dvv := Evaluate(Pq![Fq!c : c in Coefficients(Denominator(ci))], X0);
      if dvv eq 0 then return Fq!1; end if;
      s +:= (nv/dvv)*Y0^i;
    end for;
    return s;
  end function;
  for P in D1 cat Dinffin do
    if ev(P[1], P[2]) ne 0 then return ["SKIP transported point off the model"]; end if;
  end for;

  // ---- degree-1 places, keys, iota permutation ----
  pl1 := Places(F, 1);
  n1 := #pl1;
  keys := [];
  for P in pl1 do
    vx := Valuation(xx, P);
    if vx lt 0 then kx := <1, Fq!0, vx>;
    else cx := Fq!Evaluate(xx, P); kx := <0, cx, Valuation(xx - F!cx, P)>; end if;
    vy := Valuation(yy, P);
    if vy lt 0 then ky := <1, Fq!0, vy>;
    else cy := Fq!Evaluate(yy, P); ky := <0, cy, Valuation(yy - F!cy, P)>; end if;
    Append(~keys, <kx, ky>);
  end for;
  if #Set(keys) ne n1 then
    return [Sprintf("SKIP key collision (%o places, %o distinct keys)", n1, #Set(keys))];
  end if;
  idx := AssociativeArray();
  for i in [1..n1] do idx[keys[i]] := i; end for;
  iotaperm := [];
  for i in [1..n1] do
    k := keys[i]; ky := k[2];
    kk := <k[1], <ky[1], -ky[2], ky[3]>>;
    okk, j := IsDefined(idx, kk);
    if not okk then return [Sprintf("SKIP iota image of place %o missing", i)]; end if;
    Append(~iotaperm, j);
  end for;
  for i in [1..n1] do
    if iotaperm[iotaperm[i]] ne i then return ["SKIP iota not an involution"]; end if;
  end for;

  findpt := function(X0, Y0)
    hits := [i : i in [1..n1] | keys[i][1][1] eq 0 and keys[i][1][2] eq X0
                                and keys[i][2][1] eq 0 and keys[i][2][2] eq Y0];
    if #hits ne 1 then return false, 0; end if;
    return true, hits[1];
  end function;
  D1idx := []; Dinfidx := [];
  for P in D1 do
    okf, i := findpt(P[1], P[2]);
    if not okf then return ["SKIP D1 point has no unique degree-1 place"]; end if;
    Append(~D1idx, i);
  end for;
  for P in Dinffin do
    okf, i := findpt(P[1], P[2]);
    if not okf then return ["SKIP Dinf point has no unique degree-1 place"]; end if;
    Append(~Dinfidx, i);
  end for;
  bidx := [i : i in [1..n1] | keys[i][1][1] eq 0 and keys[i][1][2] eq 0
                              and keys[i][2][1] eq 0 and keys[i][2][2] eq 0];
  Append(~L, Sprintf("BOUNDARY %o vx %o", #bidx, seqstr([keys[i][1][3] : i in bidx])));
  if #bidx ne 2 then return L cat ["SKIP boundary place count wrong"]; end if;
  // both boundary places must be iota-fixed
  if iotaperm[bidx[1]] ne bidx[1] or iotaperm[bidx[2]] ne bidx[2] then
    return L cat ["SKIP boundary places not iota-fixed"];
  end if;

  // ---- class group ----
  tcl := Cputime();
  Cl, mp, cmp := ClassGroup(F);
  Tor := TorsionSubgroup(Cl);
  Append(~L, Sprintf("CLASSNUMBER %o", #Tor));
  Append(~L, Sprintf("INVARIANTS %o", seqstr(Invariants(Tor))));
  Append(~L, Sprintf("NPLACES1 %o", n1));
  Append(~L, Sprintf("CLTIME %o", Cputime(tcl)));

  Zdiv := &+[Divisor(pl1[i]) : i in D1idx]
          - &+[Divisor(pl1[i]) : i in Dinfidx]
          - &+[Divisor(pl1[i]) : i in bidx];
  if Degree(Zdiv) ne 0 then return L cat ["SKIP Zdiv degree nonzero"]; end if;
  Zc := Tor!cmp(Zdiv);
  iZdiv := &+[Divisor(pl1[iotaperm[i]]) : i in D1idx]
           - &+[Divisor(pl1[iotaperm[i]]) : i in Dinfidx]
           - &+[Divisor(pl1[iotaperm[i]]) : i in bidx];
  iZc := Tor!cmp(iZdiv);
  antiok := (iZc eq -Zc);
  Append(~L, Sprintf("ANTIINV %o", antiok select 1 else 0));
  Append(~L, Sprintf("ZORDER %o", Order(Zc)));

  useZ := antiok select Zc else Zc - iZc;
  m := Order(useZ);
  Append(~L, Sprintf("USEDOUBLE %o", antiok select 0 else 1));
  Append(~L, Sprintf("MODULUS %o", m));
  if m eq 1 then return L cat ["SKIP generator reduces to 0"]; end if;

  Zm := AbelianGroup([m]);
  hm := hom< Zm -> Tor | useZ >;

  ns := {}; offs := 0; ptl := [];
  for i in [1..n1] do
    c := Tor!cmp(Divisor(pl1[i]) - Divisor(pl1[iotaperm[i]]));
    okp, w := HasPreimage(c, hm);
    nn := -1;
    if okp then nn := (Integers()!Eltseq(w)[1]) mod m; Include(~ns, nn);
    else offs +:= 1; end if;
    k := keys[i];
    Append(~ptl, Sprintf("PT %o %o %o %o %o %o", i,
        k[1][1] eq 1 select "inf" else Sprint(k[1][2]),
        k[1][3],
        k[2][1] eq 1 select "inf" else Sprint(k[2][2]),
        k[2][3], nn));
  end for;
  Append(~L, Sprintf("OFFLINE %o", offs));
  Append(~L, "RESIDUES " cat seqstr(Sort(SetToSequence(ns))));

  // ---- ell-divisibility of the generator (saturation / index control) ----
  dvs := [];
  for ell in PrimesInInterval(2, elltop) do
    if #Tor mod ell ne 0 then continue; end if;
    hl := hom< Tor -> Tor | [ ell*Tor.i : i in [1..Ngens(Tor)] ] >;
    if HasPreimage(useZ, hl) then Append(~dvs, ell); end if;
  end for;
  Append(~L, "DIVISIBLE " cat seqstr(dvs));
  dvt := [];
  for ell in PrimesInInterval(2, elltop) do
    if #Tor mod ell eq 0 then Append(~dvt, ell); end if;
  end for;
  Append(~L, "TESTED " cat seqstr(dvt));

  return L cat ptl;
end function;

// =====================================================================
// driver: fork one child per prime
// =====================================================================
primes := [q : q in PrimesInInterval(qmin, qmax) | q mod 3 eq 1];
printf "candidate primes (q = 1 mod 3): %o\n", primes;

running := 0;
for q in primes do
  if Fork() eq 0 then
    tq := Cputime();
    L := [];
    try
      L := worker(q, elltop, sqcoef, dncoef);
    catch err
      L := ["SKIP worker error"];
    end try;
    Append(~L, Sprintf("TIME %o", Cputime(tq)));
    writelines(wdir cat "/" cat tag cat Sprintf("_q%o.txt", q), L);
    quit;
  end if;
  running +:= 1;
  if running mod nfork eq 0 then WaitForAllChildren(); end if;
end for;
WaitForAllChildren();
print "workers finished";

// ---------------------------------------------------- collect and combine
info := [];
for q in primes do
  okr, lines := readlines(wdir cat "/" cat tag cat Sprintf("_q%o.txt", q));
  if not okr then printf "q = %o : NO OUTPUT\n", q; continue; end if;
  m := 0; res := []; skip := ""; nplc := 0; anti := -1; ordz := 0; off := 0;
  inv := ""; tm := "?"; dvs := ""; tst := ""; pts := []; hnum := 0; usedbl := 0;
  for s in lines do
    if s eq "" then continue; end if;
    t := Split(s, " ");
    h := t[1];
    if h eq "SKIP" then skip := s;
    elif h eq "MODULUS" then m := StringToInteger(t[2]);
    elif h eq "RESIDUES" then
      nt := #t; res := [StringToInteger(t[i]) : i in [2..nt]];
    elif h eq "NPLACES1" then nplc := StringToInteger(t[2]);
    elif h eq "ANTIINV" then anti := StringToInteger(t[2]);
    elif h eq "ZORDER" then ordz := StringToInteger(t[2]);
    elif h eq "OFFLINE" then off := StringToInteger(t[2]);
    elif h eq "INVARIANTS" then inv := s;
    elif h eq "CLASSNUMBER" then hnum := StringToInteger(t[2]);
    elif h eq "USEDOUBLE" then usedbl := StringToInteger(t[2]);
    elif h eq "DIVISIBLE" then dvs := s;
    elif h eq "TESTED" then tst := s;
    elif h eq "TIME" then tm := t[2];
    elif h eq "PT" then Append(~pts, s);
    end if;
  end for;
  if skip ne "" then printf "q = %o : %o [%o s]\n", q, skip, tm; continue; end if;
  printf "q = %3o : #E8(F_q) = %3o , h = %o , anti-inv %o , double %o , m = %o , |N_q| = %o , off-line %o  [%o s]\n",
      q, nplc, hnum, anti, usedbl, m, #res, off, tm;
  printf "         %o\n         %o\n         %o\n", inv, dvs, tst;
  Append(~info, <q, m, res, pts, nplc, hnum>);
end for;

printf "usable primes: %o\n", #info;
if #info eq 0 then print "M612AP_SIEVE_DONE"; quit; end if;

// gcd of the class numbers bounds J(E8)(Q)_tors
gh := 0;
for t in info do gh := GCD(gh, t[6]); end for;
printf "GCD of the class numbers (bounds #J(E8)(Q)_tors): %o\n", gh;

print "=== combination through a smooth modulus ===";
cands := {};
for t in info do
  for fe in Factorization(t[2]) do
    for k in [1..fe[2]] do Include(~cands, fe[1]^k); end for;
  end for;
end for;
cands := Sort(SetToSequence(cands));
parts := [];
for pk in cands do
  users := [t : t in info | t[2] mod pk eq 0];
  if #users eq 0 then continue; end if;
  allowed := {0..pk-1};
  for t in users do
    allowed := allowed meet { n mod pk : n in t[3] };
  end for;
  printf "modulus %6o : %2o primes contribute , |allowed| = %o %o\n",
      pk, #users, #allowed,
      #allowed le 12 select Sprint(Sort(SetToSequence(allowed))) else "";
  Append(~parts, <pk, allowed, #users>);
end for;

best := AssociativeArray();
for t in parts do
  ell := Factorization(t[1])[1][1];
  okb, cur := IsDefined(best, ell);
  if not okb then best[ell] := t;
  elif #t[2] lt #cur[2] or (#t[2] eq #cur[2] and t[1] gt cur[1]) then best[ell] := t;
  end if;
end for;
Mtot := 1;
for ell in Sort(SetToSequence(Keys(best))) do
  t := best[ell];
  printf "  ell = %o : best modulus %o , allowed %o\n", ell, t[1],
      #t[2] le 12 select Sprint(Sort(SetToSequence(t[2]))) else Sprint(#t[2]);
  if t[2] eq {0} then Mtot *:= t[1]; end if;
end for;
printf "SIEVE CONCLUSION: n(P) = 0 mod %o for every rational point P of E8\n", Mtot;
print "M612AP_SIEVE_DONE";
quit;
