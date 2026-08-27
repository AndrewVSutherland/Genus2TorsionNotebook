// claude_ov_88lift_xT.m  -- Lane 5 (overnight 2026-07-25): the [8,8] lift layer on
// Nicholls' Lambda_334 family, done with the even-degree x-T machinery.
//
// STRUCTURAL POINT BEING TESTED (new this session):
//   C1 : y^2 = d2*f1,  f1 = l1*l2*l3.  On this family l2 = -x^2+bc SPLITS (bc is a
//   square identically), so the factor type of the sextic is [1,1,2,2].  The rational
//   2-torsion is still exactly (Z/2)^2 = {0,[L1],[L2],[L3]} (Galois-stable even
//   subsets of orbit sizes 1,1,2,2, modulo complement = 4 classes), i.e. 2-rank 2.
//   Therefore the 2-primary part of J1(Q)_tors has AT MOST TWO cyclic factors, so
//     J1(Q)_tors  contains  Z/8 x Z/8   =>   J1(Q)[4] = Sigma = (Z/4)^2
//                                       =>   ALL 12 order-4 elements of Sigma
//                                            lie in 2*J1(Q).
//   Hence "every one of the 12 order-4 elements passes the x-T test" is NECESSARY
//   for [8,8].  The 2026-07-18 diagnostic saw only 4 of 12 pass even the norm
//   presieve; if that is family-wide, [8,8] is impossible on this chart.
//
// x-T machinery.  L = Q[x]/f1 = prod_i Q[x]/g_i over the IRREDUCIBLE factors g_i.
//   D in 2 J1(Q)  =>  exists lambda in Q* with u_D(theta_i) = lambda*(square in K_i)
//                     for every i.
//   * deg g_i = 1, g_i = x-r :  K_i = Q, lambda == u_D(r) exactly mod squares.
//   * deg g_i = 2 irreducible:  N_i := Norm(u_D(theta_i)) must be a SQUARE; then
//     Hilbert 90 gives the rational representative
//       lambda_i = n_i / Norm(1 + u_D(theta_i)/n_i),  n_i^2 = N_i
//     (flip the sign of n_i if the argument vanishes), well defined modulo squares
//     and modulo d_i = squarefree disc(K_i) (since sqrt(d_i) is a square in K_i).
//   Lift exists  <=>  exists e_i in {1,d_i} with all lambda_i*e_i equal mod squares.
//
// Params:  MODE = "controls" | "list" | "box" | "split" ; see bottom of file.
SetColumns(0);
SetMemoryLimit(12*10^9);
Q := Rationals();  P<x> := PolynomialRing(Q);
Z := Integers();

if not assigned MODE then MODE := "controls"; end if;

// ---------------------------------------------------------------- family
function MemberSTV(sv, tv, vv)
  if sv eq 0 or tv eq 0 or tv^2 eq 1 then return false, _, _, _; end if;
  Av := sv^2 - tv^4 + tv^2;
  if Av eq 0 then return false, _, _, _; end if;
  den := -sv^2*tv*Av*vv^2 + tv;
  if den eq 0 then return false, _, _, _; end if;
  uv := (-sv^2*Av*vv^2 - 2*Av*vv - 1)/den;
  if uv^2*sv^2 + 1 - tv^2 eq 0 then return false, _, _, _; end if;
  av := Av/(1 - tv^2);
  bv := Av/(uv^2*sv^2 + 1 - tv^2);
  cv := tv^2;
  d2 := Av*(sv^2*uv^2 + tv^4 - 2*tv^2 + 1)
          *(sv^4*uv^2 - sv^2*tv^2*uv^2 + sv^2*uv^2 - tv^6 + 3*tv^4 - 3*tv^2 + 1);
  if d2 eq 0 then return false, _, _, _; end if;
  l1 := (-av + bv + cv - 1)*x^2 + (2*av - 2*bv*cv)*x + (av*bv*cv - av*bv - av*cv + bv*cv);
  l2 := -x^2 + bv*cv;
  l3 := x^2 - av;
  f1 := l1*l2*l3;
  if Degree(f1) ne 6 then return false, _, _, _; end if;
  f := d2*f1/LeadingCoefficient(f1);
  if Discriminant(f) eq 0 then return false, _, _, _; end if;
  dd := LCM([Denominator(cc) : cc in Coefficients(f)]);
  fint := P![ cc*dd^2 : cc in Coefficients(f) ];
  cont := GCD([Z!cc : cc in Coefficients(fint)]);
  sqf := 1;
  for pf in Factorization(cont) do sqf *:= pf[1]^(2*(pf[2] div 2)); end for;
  fint := P![ (Z!cc) div sqf : cc in Coefficients(fint) ];
  Ls := [ l1/LeadingCoefficient(l1), l2/LeadingCoefficient(l2), l3/LeadingCoefficient(l3) ];
  return true, fint, Ls, [av,bv,cv,Av,uv,d2];
end function;

function MemberMNV(mv, nv, vv)
  al := (mv^2-1)/(2*mv);  tv := (mv^2+1)/(2*mv);
  if nv eq 0 then return false, _, _, _; end if;
  sv := (nv^2 + al^4)/(2*nv);
  return MemberSTV(sv, tv, vv);
end function;

// ---------------------------------------------------------------- square classes
function SqFree(r)
  assert r ne 0;
  m := Numerator(r)*Denominator(r);
  s := Sign(m);  m := AbsoluteValue(m);
  res := s;
  for pf in Factorization(m) do
    if IsOdd(pf[2]) then res *:= pf[1]; end if;
  end for;
  return res;
end function;

// ---------------------------------------------------------------- the x-T test
// facs: list of irreducible monic factors of f1 (degrees 1 or 2).
// Returns  status, normpat (string over the three ORIGINAL quadratics),
//          lams, ds (per irreducible component), compat, witness
function xTData(uD, facs, Ls)
  // (a) norm presieve on the three original quadratics: Res(Lj, uD) must be a square
  normsq := [];
  for j in [1..3] do
    R := Resultant(Ls[j], uD);
    if R eq 0 then return "degenerate", "", [], [], false, 0; end if;
    Append(~normsq, IsSquare(R));
  end for;
  pat := &cat[ (b select "1" else "0") : b in normsq ];
  if not &and normsq then return "ok", pat, [], [], false, 0; end if;
  // (b) per-component rational representative
  lams := []; ds := [];
  for g in facs do
    if Degree(g) eq 1 then
      r := -Coefficient(g,0);
      gam := Evaluate(uD, r);
      if gam eq 0 then return "degenerate", pat, [], [], false, 0; end if;
      Append(~lams, SqFree(gam));  Append(~ds, 1);
    else
      K := NumberField(g);
      th := K.1;
      gam := Evaluate(uD, th);
      if gam eq 0 then return "degenerate", pat, [], [], false, 0; end if;
      Ng := Norm(gam);
      isq, ng := IsSquare(Ng);
      if not isq then return "ok", pat, [], [], false, 0; end if;  // shouldn't happen
      del := 1 + gam/ng;
      if del eq 0 then ng := -ng; del := 1 + gam/ng; end if;
      lam := ng/Norm(del);
      assert gam eq lam*del^2;
      Append(~lams, SqFree(lam));
      Append(~ds, SqFree(Discriminant(g)));
    end if;
  end for;
  // (c) compatibility over the choices e_i in {1,d_i}
  k := #lams;
  compat := false; wit := 0;
  for mask in [0 .. 2^k - 1] do
    cs := [ SqFree(lams[i]*(((mask div 2^(i-1)) mod 2) eq 0 select 1 else ds[i])) : i in [1..k] ];
    if #Seqset(cs) eq 1 then compat := true; wit := cs[1]; end if;
  end for;
  return "ok", pat, lams, ds, compat, wit;
end function;

// full report for one member
procedure Report(tag, fint, Ls)
  C := HyperellipticCurve(fint);
  J := Jacobian(C);
  fac := Factorization(fint);
  facs := [ g[1]/LeadingCoefficient(g[1]) : g in fac ];
  shape := Sort([ Degree(g) : g in facs ]);
  T, mT := TorsionSubgroup(J);
  inv := Invariants(T);
  printf "MEMBER %o factorshape=%o J1tors=%o\n", tag, shape, inv;
  n4 := 0; npass := 0; nnorm := 0; pats := [];
  for g in T do
    if Order(g) ne 4 then continue; end if;
    n4 +:= 1;
    D := mT(g);  D2 := mT(2*g);
    uD := D[1];
    if Degree(uD) ne D[3] then
      printf "  ORD4 %o DEG-MISMATCH deg(u)=%o d=%o\n", n4, Degree(uD), D[3];
      Append(~pats, "XXX");
      continue;
    end if;
    st, pat, lams, ds, compat, wit := xTData(uD, facs, Ls);
    Append(~pats, st eq "ok" select pat else "DEG");
    if st ne "ok" then printf "  ORD4 %o DEGENERATE pat=%o\n", n4, pat; continue; end if;
    if pat eq "111" then nnorm +:= 1; end if;
    if compat then npass +:= 1; end if;
    printf "  ORD4 %o 2D_u=%o normpat=%o lam=%o d=%o compat=%o wit=%o\n",
           n4, D2[1], pat, lams, ds, compat, wit;
  end for;
  ps := {* p : p in pats *};
  printf "SUMMARY %o order4=%o normpass111=%o xTpass=%o patterns=%o J1tors=%o\n",
         tag, n4, nnorm, npass, ps, inv;
end procedure;

// ---------------------------------------------------------------- positive control
if MODE eq "controls" then
  printf "=== POSITIVE CONTROL A: x-T test on D = 2E, random y^2=l1*l2*l3 (2-rank 2) ===\n";
  ncontrol := 0; nok := 0; nEfalse := 0;
  for trial in [1..60] do
    repeat
      l1 := x^2 + Random(-5,5)*x + Random(-5,5);
      l2 := x^2 + Random(-5,5)*x + Random(-5,5);
      l3 := x^2 + Random(-5,5)*x + Random(-5,5);
      f := l1*l2*l3;
    until Discriminant(f) ne 0 and #Roots(f) eq 0 and IsIrreducible(l1)
          and IsIrreducible(l2) and IsIrreducible(l3);
    C := HyperellipticCurve(f);  J := Jacobian(C);
    facs := [l1,l2,l3];  Ls := [l1,l2,l3];
    pts := Points(J : Bound := 40);
    for E in pts do
      if E eq J!0 then continue; end if;
      D := 2*E;
      if D eq J!0 or Degree(D[1]) ne D[3] then continue; end if;
      st, pat, lams, ds, compat, wit := xTData(D[1], facs, Ls);
      if st ne "ok" then continue; end if;
      ncontrol +:= 1;
      if compat then nok +:= 1; else
        printf "CONTROL_FAIL f=%o E=%o pat=%o lam=%o d=%o\n", f, E, pat, lams, ds;
      end if;
      if Degree(E[1]) eq E[3] then
        st2, p2, l2s, d2s, c2, w2 := xTData(E[1], facs, Ls);
        if st2 eq "ok" and not c2 then nEfalse +:= 1; end if;
      end if;
      if ncontrol ge 60 then break trial; end if;
    end for;
  end for;
  printf "CONTROL_A: %o classes D=2E, x-T said divisible for %o; %o of the E themselves NOT divisible\n",
         ncontrol, nok, nEfalse;

  printf "=== POSITIVE CONTROL B: factor type [1,1,2,2] (the shape of the Lambda_334 sextic) ===\n";
  nB := 0; nBok := 0; nBfalse := 0;
  for trial in [1..200] do
    repeat
      r := Random(-6,6);
      l2 := x^2 - r^2;                            // splits, like -x^2+bc with bc square
      l1 := x^2 + Random(-5,5)*x + Random(-5,5);
      l3 := x^2 + Random(-5,5)*x + Random(-5,5);
      f := l1*l2*l3;
    until r ne 0 and Discriminant(f) ne 0 and IsIrreducible(l1) and IsIrreducible(l3);
    C := HyperellipticCurve(f);  J := Jacobian(C);
    facs := [ g[1]/LeadingCoefficient(g[1]) : g in Factorization(f) ];
    Ls := [l1,l2,l3];
    pts := Points(J : Bound := 30);
    for E in pts do
      if E eq J!0 then continue; end if;
      D := 2*E;
      if D eq J!0 or Degree(D[1]) ne D[3] then continue; end if;
      st, pat, lams, ds, compat, wit := xTData(D[1], facs, Ls);
      if st ne "ok" then continue; end if;
      nB +:= 1;
      if compat then nBok +:= 1; else
        printf "CONTROL_B_FAIL f=%o E=%o pat=%o lam=%o d=%o\n", f, E, pat, lams, ds;
      end if;
      if Degree(E[1]) eq E[3] then
        st2, p2, l2s, d2s, c2, w2 := xTData(E[1], facs, Ls);
        if st2 eq "ok" and not c2 then nBfalse +:= 1; end if;
      end if;
      if nB ge 60 then break trial; end if;
    end for;
  end for;
  printf "CONTROL_B: %o classes D=2E, x-T said divisible for %o; %o of the E themselves NOT divisible\n",
         nB, nBok, nBfalse;
  printf "CONTROLS_DONE\n";
  quit;
end if;

// ---------------------------------------------------------------- member list mode
if MODE eq "list" then
  if not assigned LIST then printf "need LIST\n"; quit; end if;
  for chunk in Split(LIST, ";") do
    prts := Split(chunk, ",");
    mv := Q!eval(prts[1]); nv := Q!eval(prts[2]); vv := Q!eval(prts[3]);
    ok, fint, Ls := MemberMNV(mv, nv, vv);
    if not ok then printf "MEMBER (%o,%o,%o) INVALID\n", mv,nv,vv; continue; end if;
    tag := Sprintf("(m,n,v)=(%o,%o,%o)", mv, nv, vv);
    try
      Report(tag, fint, Ls);
    catch e
      printf "MEMBER %o ERROR %o\n", tag, e`Object;
    end try;
  end for;
  printf "LIST_DONE\n";
  quit;
end if;

// ---------------------------------------------------------------- (s,t,v) box mode
if MODE eq "box" then
  if not assigned HB then HB := 3; elif Type(HB) eq MonStgElt then HB := StringToInteger(HB); end if;
  if not assigned SZLIM then SZLIM := 3000; elif Type(SZLIM) eq MonStgElt then SZLIM := StringToInteger(SZLIM); end if;
  cnt := 0;
  vals := [Q| ];
  for a in [1..HB] do for b in [1..HB] do
    if GCD(a,b) eq 1 then Include(~vals, Q!a/b); Include(~vals, -Q!a/b); end if;
  end for; end for;
  pos := [ w : w in vals | w gt 0 ];
  for sv in pos do for tv in pos do for vv in vals do
    ok, fint, Ls := MemberSTV(sv, tv, vv);
    if not ok then continue; end if;
    if #Sprint(Coefficients(fint)) gt SZLIM then continue; end if;
    cnt +:= 1;
    tag := Sprintf("(s,t,v)=(%o,%o,%o)", sv, tv, vv);
    try
      Report(tag, fint, Ls);
    catch e
      printf "MEMBER %o ERROR %o\n", tag, e`Object;
    end try;
    printf "PROGRESS box count=%o\n", cnt;
  end for; end for; end for;
  printf "BOX_DONE count=%o\n", cnt;
  quit;
end if;

// -------------------------------------------------- is l2 = -x^2+bc always split?
if MODE eq "split" then
  if not assigned HB then HB := 5; elif Type(HB) eq MonStgElt then HB := StringToInteger(HB); end if;
  vals := [Q| ];
  for a in [1..HB] do for b in [1..HB] do
    if GCD(a,b) eq 1 then Include(~vals, Q!a/b); Include(~vals, -Q!a/b); end if;
  end for; end for;
  tot := 0; nsplit := 0; nAsq := 0; nbc := 0; shapes := {* *};
  for sv in vals do for tv in vals do for vv in vals do
    if sv le 0 or tv le 0 then continue; end if;
    ok, fint, Ls, aux := MemberSTV(sv, tv, vv);
    if not ok then continue; end if;
    av := aux[1]; bv := aux[2]; cv := aux[3]; Av := aux[4];
    tot +:= 1;
    if IsSquare(bv*cv) then nsplit +:= 1; end if;
    if IsSquare(bv) then nbc +:= 1; end if;
    if IsSquare(av) then nAsq +:= 1; end if;
    Include(~shapes, Sort([Degree(g[1]) : g in Factorization(fint)]));
  end for; end for; end for;
  printf "SPLITCHECK total=%o bc_square=%o b_square=%o a_square=%o shapes=%o\n",
         tot, nsplit, nbc, nAsq, shapes;
  printf "SPLIT_DONE\n";
  quit;
end if;

printf "unknown MODE\n";
quit;
