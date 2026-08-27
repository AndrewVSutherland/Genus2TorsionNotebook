//////////////////////////////////////////////////////////////////////
// claude_ov_m612ap_sieve2.m      lane 9 ([6,12])   2026-07-25
//
// THE ABEL-PRYM MORDELL-WEIL SIEVE for the [6,12] gatekeeper curve E8.
// Rewritten from claude_ov_m612ap_sieve.m so that the rank-1 Prym class
// is INTRINSIC (explicit rational divisors on E8) instead of transported
// through the uncertified bigonal correspondence.
//
//   E8 : Q8(x,y) = y^8 + A(x) y^4 + B(x) y^2 + C(x) = 0,   genus 4,
//        A = 216x^4+72x^3-24x^2,
//        B = -1296x^6-1728x^5-432x^4+64x^3,
//        C = -3888x^8-2592x^7+432x^6+288x^5-48x^4,
//        x = e, y = e*w  (mp8(w,e) = mp4(w^2,e), data/contact6_m612_E8_mp8Q.txt)
//   iota : y -> -y  is the Prym involution, E8/iota = E4 (genus 2).
//   Riemann-Hurwitz: 2*4-2 = 2*(2*2-2) + r  =>  r = 2, so iota has
//   EXACTLY two fixed points; they are the two rational places over e=0
//   (ramification 2 and 6 in E8 -> P^1_e, hence individually iota-stable).
//
//   alpha : E8 -> Pr := (1-iota) J(E8),  P |-> [P - iota P].
//   alpha(P) = 0 and P != iota P would give a degree-1 function on a
//   genus-4 curve; so  alpha(P) = 0  <=>  P in {b1,b2}.  Therefore
//        alpha(E8(Q)) = {0}   <=>   E8(Q) = {b1, b2}.
//   rank Pr(Q) = 1 (code/contact6_m612_prym_rank_verifier.m), so
//   alpha(P) = n(P)*G + torsion for one integer n(P).
//
// THE INTRINSIC PRYM CLASSES.  Q8 is even in y: Q8(x0,y) = h(x0,y^2).
// An irreducible factor g of Q8(x0,y) over Q that is NOT even in y is
// swapped with g(-y) by iota, so, with P_g the corresponding place,
//        W(x0,g) := [P_g - P_{g(-y)}] = (1-iota)[P_g]  in  Pr(Q).
// code/claude_ov_m612ap_wfind.m scanned every x0 = a/b of height <= 30
// and found exactly three such x0, each with a pair of cubic factors:
//        x0 = -2/3  : t^3 - 2t^2 + (4/3)t + 8/3
//        x0 = -2/5  : t^3 - (6/5)t^2 + (12/25)t - 152/125
//        x0 = 2/15  : t^3 - (2/5)t^2 + (4/75)t - 88/375
// Lambda := <W1,W2,W3>  is an explicitly given subgroup of Pr(Q); it is
// reduced mod q by reducing the DIVISORS (fibres over x0 are etale, so
// each place reduces to the set of places whose y-residue kills g mod q).
//
// PER PRIME q THIS COMPUTES
//   1. genus 4 check (good reduction),
//   2. the two boundary places over x=0 and the iota-permutation of the
//      degree-1 places (signature-based, with an injectivity assert),
//   3. Cl^0 = J(E8)(F_q) = TorsionSubgroup(ClassGroup(F)),
//   4. W1,W2,W3 in J(E8)(F_q) and S_q := <W1,W2,W3>,
//   5. for every degree-1 place R : c_R := [R - iota R], whether
//      m*c_R in S_q for m = 1..mtop, and the discrete log in <W1> if any,
//   6. ell-divisibility of each W_i in J(E8)(F_q) for ell <= elltop
//      (a single NON-divisibility proves ell does not divide the index
//      [Pr(Q)/tors : Lambda/tors]).
//
// Usage:
//   magma -b qmin:=5 qmax:=200 nfork:=12 wdir:=results tag:=ap2 \
//         code/claude_ov_m612ap_sieve2.m
//   magma -b single:=31 code/claude_ov_m612ap_sieve2.m     (one prime, verbose)
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(4*10^9);

if not assigned qmax   then qmax := 200;   elif Type(qmax) eq MonStgElt then qmax := StringToInteger(qmax); end if;
if not assigned qmin   then qmin := 5;     elif Type(qmin) eq MonStgElt then qmin := StringToInteger(qmin); end if;
if not assigned nfork  then nfork := 12;   elif Type(nfork) eq MonStgElt then nfork := StringToInteger(nfork); end if;
if not assigned wdir   then wdir := "results"; end if;
if not assigned tag    then tag := "ap2";  end if;
if not assigned elltop then elltop := 200; elif Type(elltop) eq MonStgElt then elltop := StringToInteger(elltop); end if;
if not assigned mtop   then mtop := 24;    elif Type(mtop) eq MonStgElt then mtop := StringToInteger(mtop); end if;
if not assigned single then single := 0;   elif Type(single) eq MonStgElt then single := StringToInteger(single); end if;

QQ := Rationals();
Pt<t> := PolynomialRing(QQ);

// the three intrinsic Prym classes: <x0, g>
GENS := [
  < QQ!(-2)/3,  t^3 - 2*t^2 + (QQ!4/3)*t + QQ!8/3 >,
  < QQ!(-2)/5,  t^3 - (QQ!6/5)*t^2 + (QQ!12/25)*t - QQ!152/125 >,
  < QQ!(2)/15,  t^3 - (QQ!2/5)*t^2 + (QQ!4/75)*t - QQ!88/375 >
];

// sanity over Q: each g must divide Q8(x0,y) and be non-even
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
printf "=== claude_ov_m612ap_sieve2 : intrinsic Abel-Prym MW sieve on E8 ===\n";
printf "generators verified over Q: %o cubic place pairs\n", #GENS;

writelines := procedure(fn, L)
  fh := Open(fn, "w");
  for s in L do Puts(fh, s); end for;
  Flush(fh);
  delete fh;
end procedure;

readlines := function(fn)
  ok := true; s := "";
  try s := Read(fn); catch err ok := false; end try;
  if not ok then return false, []; end if;
  return true, Split(s, "\n");
end function;

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
  // An iota-EQUIVARIANT key.  The "position" data (x-value / infinity,
  // ramification index over the x-line) is iota-invariant; the function
  // phi (= y at a finite place, = y/x at an infinite one) is
  // anti-invariant, and the uniformiser u (= x-cx, resp. 1/x) is
  // invariant.  Hence expanding phi = sum c_j u^j and negating every c_j
  // yields exactly the key of iota(R).  (mkkey(R,true) = key of iota(R).)
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
  // Riemann-Hurwitz forces exactly two geometric fixed points, and both
  // are the rational boundary places; any other count means the key is
  // not faithful (or the reduction is bad).
  if nfix ne 2 then return L cat [Sprintf("SKIP iota fixes %o degree-1 places (want 2)", nfix)]; end if;
  for D in dec0 do
    ok0 := false;
    for i in [1..n1] do if pl1[i] eq D then ok0 := (iot[i] eq i); end if; end for;
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

  // ---------------- reduce the three global Prym classes -------------
  ws := [];
  for gg in GENS do
    x0 := gg[1]; g := gg[2];
    if q eq 2 or q eq 3 or q eq 5 then return L cat ["SKIP tiny prime"]; end if;
    dn := LCM([Denominator(c) : c in Coefficients(g)] cat [Denominator(x0)]);
    if dn mod q eq 0 then return L cat ["SKIP q divides a generator denominator"]; end if;
    x0q := kq!Numerator(x0) / kq!Denominator(x0);
    Q8t := TT^8 + (kq!Acf(x0))*TT^4 + (kq!Bcf(x0))*TT^2 + (kq!Ccf(x0));
    if Discriminant(Q8t) eq 0 then return L cat ["SKIP generator fibre ramified mod q"]; end if;
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
      return L cat [Sprintf("SKIP reduction degree %o/%o (want 3/3)", dgp, dgm)];
    end if;
    Append(~ws, clsof(Dp - Dm));
  end for;
  ords := [Order(w) : w in ws];
  Append(~L, Sprintf("WORDERS %o", seqstr(ords)));

  S := sub< J | ws >;
  Append(~L, Sprintf("SORDER %o", #S));
  Append(~L, Sprintf("SINV %o", seqstr(Invariants(S))));

  // discrete-log helper in the cyclic group <w>
  dlog := function(w, c)
    m := Order(w);
    if m eq 1 then return (c eq J!0) select 0 else -1; end if;
    okp, e := HasPreimage(c, hom< AbelianGroup([m]) -> J | w >);
    if not okp then return -1; end if;
    return (Integers()!Eltseq(e)[1]) mod m;
  end function;

  // relations among the three global generators, read mod q
  for i in [1..#ws] do
    Append(~L, Sprintf("DLOGW%o %o", i,
      seqstr([ dlog(ws[i], ws[j]) : j in [1..#ws] ])));
  end for;

  // ---------------- Abel-Prym classes of the F_q points --------------
  inS := []; minS := []; res := [ [] : i in [1..#ws] ];
  for i in [1..n1] do
    if iot[i] eq i then continue; end if;
    c := clsof(Divisor(pl1[i]) - Divisor(pl1[iot[i]]));
    if c in S then Append(~inS, i); end if;
    for j in [1..#ws] do
      d := dlog(ws[j], c);
      if d ge 0 then Append(~res[j], d); end if;
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
    Append(~L, Sprintf("RES%o %o", j, seqstr(Sort(Setseq(Seqset(res[j]) join {0})))));
  end for;
  Append(~L, Sprintf("M1 %o", ords[1]));
  Append(~L, Sprintf("MINMULT %o", seqstr(Sort(minS))));

  // ---------------- Prym order diagnostic ---------------------------
  // E4 : Z^4 + (216x^2+72x-24) Z^2 + (-1296x^3-1728x^2-432x+64) Z
  //           + (-3888x^4-2592x^3+432x^2+288x-48) = 0,  Z = y^2/x.
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
  dvs := []; tst := [];
  for ell in PrimesInInterval(2, elltop) do
    if #J mod ell ne 0 then continue; end if;
    Append(~tst, ell);
    hl := hom< J -> J | [ ell*J.i : i in [1..Ngens(J)] ] >;
    // ell divides the index only if EVERY generator of S is ell-divisible
    alld := true;
    for w in ws do
      if not HasPreimage(w, hl) then alld := false; break; end if;
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
  print "AP2_SINGLE_DONE";
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

// =====================================================================
// collect
// =====================================================================
info := [];
for q in primes do
  okr, lines := readlines(wdir cat "/" cat tag cat Sprintf("_q%o.txt", q));
  if not okr then printf "q = %-4o NO OUTPUT\n", q; continue; end if;
  D := AssociativeArray();
  skip := "";
  for s in lines do
    if s eq "" then continue; end if;
    tk := Split(s, " ");
    if tk[1] eq "SKIP" then skip := s; continue; end if;
    D[tk[1]] := #tk ge 2 select tk[2] else "";
  end for;
  if skip ne "" then
    printf "q = %-4o %o  [%o s]\n", q, skip, IsDefined(D,"TIME") select D["TIME"] else "?";
    continue;
  end if;
  gt := function(k) return IsDefined(D,k) select D[k] else "-"; end function;
  printf "q = %-4o #E8(F_q)=%-4o #J=%-16o #Prym=%-10o S=%-12o nonfix=%-4o inS=%-10o div=%-6o [%o s]\n",
      q, gt("NPLACES1"), gt("JORDER"), gt("PRYMORD"), gt("SORDER"), gt("NONFIXED"),
      gt("INS"), gt("DIVISIBLE"), gt("TIME");
  printf "        JINV=%o  WORDERS=%o  SINV=%o  DLOGW1=%o DLOGW2=%o DLOGW3=%o\n",
      gt("JINV"), gt("WORDERS"), gt("SINV"), gt("DLOGW1"), gt("DLOGW2"), gt("DLOGW3");
  printf "        RES1=%o  RES2=%o  RES3=%o  MINMULT=%o  TESTED=%o\n",
      gt("RES1"), gt("RES2"), gt("RES3"), gt("MINMULT"), gt("TESTED");
  Append(~info, <q, D>);
end for;

printf "usable primes: %o\n", #info;
if #info eq 0 then print "AP2_SIEVE_DONE"; quit; end if;

// gcd of #J(F_q) bounds #J(E8)(Q)_tors, hence the Prym torsion
gh := 0; mlcm := 1; allempty := true; divever := {};
for tt in info do
  gh := GCD(gh, StringToInteger(tt[2]["JORDER"]));
  mlcm := LCM(mlcm, StringToInteger(tt[2]["M1"]));
  if tt[2]["INS"] ne "-" then allempty := false; end if;
  if tt[2]["DIVISIBLE"] ne "-" then
    for s in Split(tt[2]["DIVISIBLE"], ",") do Include(~divever, StringToInteger(s)); end for;
  end if;
end for;
printf "GCD of #J(E8)(F_q) over usable primes (bounds #J(E8)(Q)_tors): %o\n", gh;
printf "LCM of ord(W1 mod q): %o\n", mlcm;
printf "primes ell for which SOME q could not rule out ell | index: %o\n", Sort(SetToSequence(divever));
if allempty then
  print "SIEVE CONCLUSION: NO non-iota-fixed point of E8(F_q) has its Abel-Prym class in Lambda mod q, for any usable q.";
  print "  => for every rational point P of E8 outside {b1,b2}, alpha(P) is NOT in Lambda.";
  print "  => if Lambda is saturated in Pr(Q) (index 1) then E8(Q) = {b1,b2}.";
else
  print "SIEVE CONCLUSION: some residue classes survive; see the INS/RESIDUES columns.";
end if;
print "AP2_SIEVE_DONE";
quit;
