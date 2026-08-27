//////////////////////////////////////////////////////////////////////
// claude_ov_m612ap_sieve3.m      lane 9 ([6,12])   2026-07-25
//
// THE ABEL-PRYM MORDELL-WEIL SIEVE on the [6,12] gatekeeper curve E8,
// third revision.  Same mathematics as claude_ov_m612ap_sieve2.m; three
// changes, each of which fixes a defect of the sieve2 run:
//
//  (a) PER-GENERATOR SKIP.  sieve2 discarded a whole prime as soon as ONE
//      of the three fibres x = x0 became ramified mod q.  The bad primes
//      are  gen1 (x0=-2/3): 2,3,5 ; gen2 (x0=-2/5): 2,3,5,19 ;
//      gen3 (x0=2/15): 2,3,5,7,11  (discriminants factored in PARI).
//      Now each generator is handled separately, so q = 7, 11, 19 become
//      usable primes instead of being thrown away.
//
//  (b) TORSION-ROBUST RESIDUES.  sieve2 recorded only
//          RES_j = { d : [R - iota R] = d * Wbar_j },
//      which is the correct constraint ONLY if Pr(Q) is torsion free.
//      The gcd of #J(E8)(F_q) over the sieve2 primes is 9, so Pr(Q) may
//      have 3-torsion, and then alpha(P) = d*W_j + t with 3t = 0 is NOT
//      detected by RES_j at all -- the sieve2 conclusion was unsound.
//      Since 9*t = 0 for any t of order dividing 9, we also record
//          RES9_j = { e : 9*[R - iota R] = e * Wbar_j },
//      and the sieve on   9d = e (mod ord Wbar_j)   is valid whatever the
//      (order-dividing-9) torsion is.
//
//  (c) the collection stage of sieve2 died on a Magma parse error
//      ("return X select Y else Z" needs parentheses), so the sieve2 run
//      produced per-prime files and NO conclusion.  Fixed here; the
//      combination is done by code/claude_ov_m612ap_combine2.py anyway.
//
// Also recorded: TBORD, the order of [b1 - b2] in J(E8)(F_q) (b1,b2 are
// the two rational iota-fixed boundary places), which lower-bounds the
// rational torsion actually present.
//
//   E8 : Q8(x,y) = y^8 + A(x) y^4 + B(x) y^2 + C(x) = 0,   genus 4,
//        A = 216x^4+72x^3-24x^2,
//        B = -1296x^6-1728x^5-432x^4+64x^3,
//        C = -3888x^8-2592x^7+432x^6+288x^5-48x^4.
//   iota : y -> -y ;  E8/iota = E4 (genus 2) ;  Riemann-Hurwitz => iota
//   has exactly 2 fixed points, the two rational places over x = 0.
//   alpha : P |-> [P - iota P] in Pr := (1-iota) J(E8) ;  alpha(P) = 0
//   iff P in {b1,b2} (else a degree-1 function on a genus-4 curve).
//   rank Pr(Q) = 1  (code/contact6_m612_prym_rank_verifier.m).
//
// Usage:
//   code/claude_magma_slot.sh -b qmin:=5 qmax:=103 nfork:=8 memgb:=12 \
//        wdir:=results/ap3 tag:=ap3 code/claude_ov_m612ap_sieve3.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned qmax   then qmax := 200;   elif Type(qmax) eq MonStgElt then qmax := StringToInteger(qmax); end if;
if not assigned qmin   then qmin := 5;     elif Type(qmin) eq MonStgElt then qmin := StringToInteger(qmin); end if;
if not assigned nfork  then nfork := 8;    elif Type(nfork) eq MonStgElt then nfork := StringToInteger(nfork); end if;
if not assigned memgb  then memgb := 12;   elif Type(memgb) eq MonStgElt then memgb := StringToInteger(memgb); end if;
if not assigned wdir   then wdir := "results/ap3"; end if;
if not assigned tag    then tag := "ap3";  end if;
if not assigned elltop then elltop := 200; elif Type(elltop) eq MonStgElt then elltop := StringToInteger(elltop); end if;
if not assigned mtop   then mtop := 24;    elif Type(mtop) eq MonStgElt then mtop := StringToInteger(mtop); end if;
if not assigned single then single := 0;   elif Type(single) eq MonStgElt then single := StringToInteger(single); end if;

SetMemoryLimit(memgb*10^9);

QQ := Rationals();
Pt<t> := PolynomialRing(QQ);

// the three intrinsic Prym classes: <x0, g>
GENS := [
  < QQ!(-2)/3,  t^3 - 2*t^2 + (QQ!4/3)*t + QQ!8/3 >,
  < QQ!(-2)/5,  t^3 - (QQ!6/5)*t^2 + (QQ!12/25)*t - QQ!152/125 >,
  < QQ!(2)/15,  t^3 - (QQ!2/5)*t^2 + (QQ!4/75)*t - QQ!88/375 >
];

Acf := func< X | 216*X^4+72*X^3-24*X^2 >;
Bcf := func< X | -1296*X^6-1728*X^5-432*X^4+64*X^3 >;
Ccf := func< X | -3888*X^8-2592*X^7+432*X^6+288*X^5-48*X^4 >;
for gg in GENS do
  x0 := gg[1]; g := gg[2];
  Q8t := t^8 + Acf(x0)*t^4 + Bcf(x0)*t^2 + Ccf(x0);
  error if Q8t mod g ne 0, "generator polynomial does not divide the fibre";
  error if Discriminant(Q8t) eq 0, "ramified fibre";
  error if IsIrreducible(g) eq false, "generator polynomial not irreducible";
  gm := Evaluate(g, -t);
  error if gm eq g or gm eq -g, "generator polynomial is even/odd in y";
  error if Q8t mod gm ne 0, "mirror polynomial does not divide the fibre";
end for;
printf "=== claude_ov_m612ap_sieve3 : intrinsic Abel-Prym MW sieve on E8 ===\n";
printf "generators verified over Q: %o cubic place pairs\n", #GENS;
printf "memory limit per process: %o GB, nfork = %o\n", memgb, nfork;

writelines := procedure(fn, L)
  fh := Open(fn, "w");
  for s in L do Puts(fh, s); end for;
  Flush(fh);
  delete fh;
end procedure;

seqstr := function(S)
  if #S eq 0 then return "-"; end if;
  return &cat[ IntegerToString(u) cat "," : u in S ];
end function;

// =====================================================================
// per-prime worker.  Returns a list of tagged lines.
// =====================================================================
worker := function(q)
  L := [];
  kq := GF(q);
  Fx<X> := RationalFunctionField(kq);
  Py<Y> := PolynomialRing(Fx);
  Q8 := Y^8 + (216*X^4+72*X^3-24*X^2)*Y^4
      + (-1296*X^6-1728*X^5-432*X^4+64*X^3)*Y^2
      + (-3888*X^8-2592*X^7+432*X^6+288*X^5-48*X^4);
  if not IsIrreducible(Q8) then return ["SKIP reducible model"]; end if;
  F<yy> := FunctionField(Q8);
  if Genus(F) ne 4 then return [Sprintf("SKIP genus %o", Genus(F))]; end if;
  xx := F!X;
  Pk<TT> := PolynomialRing(kq);

  // ---------------- boundary places over x = 0 ----------------------
  z0 := Zeros(Fx.1)[1];
  dec0 := Decomposition(F, z0);
  dec0 := [ Type(D) eq Tup select D[1] else D : D in dec0 ];
  if #dec0 ne 2 then return [Sprintf("SKIP %o places over x=0", #dec0)]; end if;
  if {Degree(D) : D in dec0} ne {1} then return ["SKIP boundary places not rational"]; end if;
  rams := Sort([RamificationIndex(D) : D in dec0]);
  Append(~L, Sprintf("BOUNDARY rams %o", seqstr(rams)));
  if rams ne [2,6] then return L cat ["SKIP boundary ramification not 2,6"]; end if;

  // ---------------- degree-1 places and the iota permutation ---------
  pl1 := Places(F, 1);
  n1 := #pl1;
  Append(~L, Sprintf("NPLACES1 %o", n1));
  yox := yy/xx;
  depth := 6;
  mkkey := function(R, neg)
    vx := Valuation(xx, R);
    if vx lt 0 then
      inff := 1; u := 1/xx; phi := yox; cxs := "inf";
    else
      inff := 0; cx := Evaluate(xx, R); u := xx - F!cx; phi := yy; cxs := Sprint(cx);
    end if;
    ex := Valuation(u, R);
    if ex eq 1 then
      cs := []; cur := phi;
      for j in [0..depth] do
        c := Evaluate(cur, R);
        Append(~cs, neg select -c else c);
        cur := (cur - F!c)/u;
      end for;
      return Sprintf("%o|%o|%o|A|%o", inff, cxs, ex, cs);
    else
      c0 := Evaluate(phi, R);
      v0 := Valuation(phi - F!c0, R);
      return Sprintf("%o|%o|%o|B|%o|%o", inff, cxs, ex, neg select -c0 else c0, v0);
    end if;
  end function;
  sig := [ mkkey(R, false) : R in pl1 ];
  if #Set(sig) ne n1 then return L cat ["SKIP signature collision"]; end if;
  idx := AssociativeArray();
  for i in [1..n1] do idx[sig[i]] := i; end for;
  iot := [];
  for i in [1..n1] do
    okk, j := IsDefined(idx, mkkey(pl1[i], true));
    if not okk then return L cat [Sprintf("SKIP iota image of place %o missing", i)]; end if;
    Append(~iot, j);
  end for;
  for i in [1..n1] do
    if iot[iot[i]] ne i then return L cat ["SKIP iota not an involution"]; end if;
  end for;
  nfix := #[i : i in [1..n1] | iot[i] eq i];
  Append(~L, Sprintf("IOTAFIX %o", nfix));
  if nfix ne 2 then return L cat [Sprintf("SKIP iota fixes %o degree-1 places (want 2)", nfix)]; end if;
  bidx := [];
  for D in dec0 do
    ok0 := false;
    for i in [1..n1] do if pl1[i] eq D then ok0 := (iot[i] eq i); Append(~bidx, i); end if; end for;
    if not ok0 then return L cat ["SKIP a boundary place is not iota-fixed"]; end if;
  end for;

  // ---------------- class group -------------------------------------
  tcl := Cputime();
  Cl, mCl := ClassGroup(F);
  J := TorsionSubgroup(Cl);
  Append(~L, Sprintf("CLTIME %o", Cputime(tcl)));
  Append(~L, Sprintf("JORDER %o", #J));
  Append(~L, Sprintf("JINV %o", seqstr(Invariants(J))));

  clsof := function(D)
    return J ! (D @@ mCl);
  end function;

  // order of [b1 - b2]: a rational torsion class of J(E8) coming from the
  // two rational boundary points.  (It is iota-INVARIANT, so it is not in
  // the Prym; recorded to see the rational torsion that is really there.)
  tb := clsof(Divisor(pl1[bidx[1]]) - Divisor(pl1[bidx[2]]));
  Append(~L, Sprintf("TBORD %o", Order(tb)));

  // ---------------- reduce the global Prym classes -------------------
  // per-generator: unusable generators get order 0 and no residues.
  ws := []; gok := [];
  for gi in [1..#GENS] do
    x0 := GENS[gi][1]; g := GENS[gi][2];
    bad := false;
    dn := LCM([Denominator(c) : c in Coefficients(g)] cat [Denominator(x0)]);
    if q le 5 or dn mod q eq 0 then bad := true; end if;
    if not bad then
      Q8t := TT^8 + (kq!Acf(x0))*TT^4 + (kq!Bcf(x0))*TT^2 + (kq!Ccf(x0));
      if Discriminant(Q8t) eq 0 then bad := true; end if;
    end if;
    if bad then Append(~ws, J!0); Append(~gok, 0); continue; end if;
    x0q := kq!Numerator(x0) / kq!Denominator(x0);
    gq  := Pk![kq!c : c in Coefficients(g)];
    gmq := Evaluate(gq, -TT);
    zx := Zeros(Fx.1 - x0q)[1];
    decx := Decomposition(F, zx);
    decx := [ Type(D) eq Tup select D[1] else D : D in decx ];
    Dp := DivisorGroup(F) ! 0; Dm := DivisorGroup(F) ! 0; dgp := 0; dgm := 0;
    for P in decx do
      RF := ResidueClassField(P);
      rr := Evaluate(yy, P);
      vp := &+[ (RF!Coefficient(gq,i))  * rr^i : i in [0..Degree(gq)] ];
      vm := &+[ (RF!Coefficient(gmq,i)) * rr^i : i in [0..Degree(gmq)] ];
      if vp eq 0 then Dp +:= Divisor(P); dgp +:= Degree(P); end if;
      if vm eq 0 then Dm +:= Divisor(P); dgm +:= Degree(P); end if;
    end for;
    if dgp ne 3 or dgm ne 3 then
      Append(~ws, J!0); Append(~gok, 0); continue;
    end if;
    Append(~ws, clsof(Dp - Dm)); Append(~gok, 1);
  end for;
  Append(~L, Sprintf("GENOK %o", seqstr(gok)));
  if &+gok eq 0 then return L cat ["SKIP no usable generator"]; end if;
  ords := [ gok[i] eq 1 select Order(ws[i]) else 0 : i in [1..#ws] ];
  Append(~L, Sprintf("WORDERS %o", seqstr(ords)));

  S := sub< J | [ws[i] : i in [1..#ws] | gok[i] eq 1] >;
  Append(~L, Sprintf("SORDER %o", #S));
  Append(~L, Sprintf("SINV %o", seqstr(Invariants(S))));

  dlog := function(w, c)
    m := Order(w);
    if m eq 1 then return (c eq J!0) select 0 else -1; end if;
    okp, e := HasPreimage(c, hom< AbelianGroup([m]) -> J | w >);
    if not okp then return -1; end if;
    return (Integers()!Eltseq(e)[1]) mod m;
  end function;

  for i in [1..#ws] do
    if gok[i] eq 0 then Append(~L, Sprintf("DLOGW%o -", i)); continue; end if;
    Append(~L, Sprintf("DLOGW%o %o", i,
      seqstr([ gok[j] eq 1 select dlog(ws[i], ws[j]) else -2 : j in [1..#ws] ])));
  end for;

  // ---------------- Abel-Prym classes of the F_q points --------------
  inS := []; minS := [];
  res  := [ [] : i in [1..#ws] ];
  res3 := [ [] : i in [1..#ws] ];
  res9 := [ [] : i in [1..#ws] ];
  for i in [1..n1] do
    if iot[i] eq i then continue; end if;
    c := clsof(Divisor(pl1[i]) - Divisor(pl1[iot[i]]));
    c3 := 3*c; c9 := 9*c;
    if c in S then Append(~inS, i); end if;
    for j in [1..#ws] do
      if gok[j] eq 0 then continue; end if;
      d := dlog(ws[j], c);
      if d ge 0 then Append(~res[j], d); end if;
      d3 := dlog(ws[j], c3);
      if d3 ge 0 then Append(~res3[j], d3); end if;
      d9 := dlog(ws[j], c9);
      if d9 ge 0 then Append(~res9[j], d9); end if;
    end for;
    mm := 0;
    for m in [1..mtop] do
      if (m*c) in S then mm := m; break; end if;
    end for;
    Append(~minS, mm);
  end for;
  Append(~L, Sprintf("NONFIXED %o", n1 - nfix));
  Append(~L, Sprintf("INS %o", seqstr(inS)));
  // the two iota-fixed places b1,b2 are rational points of E8 with
  // alpha = 0, so residue 0 is always attained.
  for j in [1..#ws] do
    if gok[j] eq 0 then
      Append(~L, Sprintf("RES%o -", j)); Append(~L, Sprintf("RES3%o -", j));
      Append(~L, Sprintf("RES9%o -", j)); continue;
    end if;
    Append(~L, Sprintf("RES%o %o", j, seqstr(Sort(Setseq(Seqset(res[j]) join {0})))));
    Append(~L, Sprintf("RES3%o %o", j, seqstr(Sort(Setseq(Seqset(res3[j]) join {0})))));
    Append(~L, Sprintf("RES9%o %o", j, seqstr(Sort(Setseq(Seqset(res9[j]) join {0})))));
  end for;
  Append(~L, Sprintf("MINMULT %o", seqstr(Sort(minS))));

  // ---------------- Prym order diagnostic ---------------------------
  try
    Q4 := Y^4 + (216*X^2+72*X-24)*Y^2 + (-1296*X^3-1728*X^2-432*X+64)*Y
        + (-3888*X^4-2592*X^3+432*X^2+288*X-48);
    F4 := FunctionField(Q4);
    if Genus(F4) eq 2 then
      Cl4 := ClassGroup(F4);
      h4 := #TorsionSubgroup(Cl4);
      Append(~L, Sprintf("JE4ORD %o", h4));
      if #J mod h4 eq 0 then Append(~L, Sprintf("PRYMORD %o", #J div h4)); end if;
    else
      Append(~L, Sprintf("JE4ORD bad-genus-%o", Genus(F4)));
    end if;
  catch err
    Append(~L, "JE4ORD error");
  end try;

  // ---------------- ell-divisibility (saturation) of the W_i ---------
  // ell can divide [Pr(Q) : Lambda] only if EVERY usable generator is
  // ell-divisible in J(E8)(F_q); a single failure at a single q kills ell.
  dvs := []; tst := [];
  for ell in PrimesInInterval(2, elltop) do
    if #J mod ell ne 0 then continue; end if;
    Append(~tst, ell);
    hl := hom< J -> J | [ ell*J.i : i in [1..Ngens(J)] ] >;
    alld := true;
    for i in [1..#ws] do
      if gok[i] eq 0 then continue; end if;
      if not HasPreimage(ws[i], hl) then alld := false; break; end if;
    end for;
    if alld then Append(~dvs, ell); end if;
  end for;
  Append(~L, "DIVISIBLE " cat seqstr(dvs));
  Append(~L, "TESTED " cat seqstr(tst));
  return L;
end function;

// =====================================================================
// single-prime verbose mode
// =====================================================================
if single ne 0 then
  tq := Cputime();
  L := worker(single);
  for s in L do printf "q=%o %o\n", single, s; end for;
  printf "q=%o TIME %o\n", single, Cputime(tq);
  print "AP3_SINGLE_DONE";
  quit;
end if;

// =====================================================================
// driver: fork one child per prime
// =====================================================================
primes := [q : q in PrimesInInterval(qmin, qmax)];
printf "primes: %o\n", primes;

running := 0;
for q in primes do
  if Fork() eq 0 then
    tq := Cputime();
    L := [];
    try
      L := worker(q);
    catch err
      L := ["SKIP worker error " cat Sprint(err`Object)];
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
print "AP3_SIEVE_DONE";
quit;
