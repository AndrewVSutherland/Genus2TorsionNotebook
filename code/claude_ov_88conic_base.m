// claude_ov_88conic_base.m -- Lane 5 (2026-07-25).  ONE (m,n) base of the Nicholls
// Lambda_334 double-stage-1 locus: find the v's, then extract the FULL x-T lambda
// profile of all 12 order-4 elements of Sigma on J1 -- i.e. the data of the
// "common-lambda compatibility" conic layer that the [8,8] route needs.
//
// Geometry of the component algebra (established this session):
//   f1 = l1*l2*l3 with l2 = -x^2+bc and bc a SQUARE identically on the family,
//   so L = Q[x]/f1 = Q x Q x K1 x K3 with K1 = Q(sqrt(D1)), D1 = disc(l1),
//   K3 = Q(sqrt(a)).  Because the two Q-components admit NO square-class ambiguity,
//   the common lambda is FORCED: lambda = lambda_2 := u_D(rho) mod squares
//   (and the norm presieve at l2 is exactly lambda_+ = lambda_-).  The lift then
//   exists iff
//        def1 := lambda_2*lambda_1 in {1, D1}   and   def3 := lambda_2*lambda_3 in {1, a}
//   modulo squares -- two conic conditions.  This script measures def1, def3.
//
// Params: m, n (strings), OUT (file), PSB (point-search bound, default 0 = use MW),
//         NV (max v per base, default 16)
SetColumns(0);
SetMemoryLimit(3*10^9);
Q := Rationals();  P<x> := PolynomialRing(Q);
Z := Integers();
if not assigned m or not assigned n then printf "need m,n\n"; quit; end if;
if not assigned NV then NV := 16; elif Type(NV) eq MonStgElt then NV := StringToInteger(NV); end if;

mv := Q!eval(m);  nv := Q!eval(n);

function SqFree(r)
  if r eq 0 then return 0; end if;
  q := Numerator(r)*Denominator(r);
  s := Sign(q);  q := AbsoluteValue(q);
  res := s;
  for pf in Factorization(q) do if IsOdd(pf[2]) then res *:= pf[1]; end if; end for;
  return res;
end function;

// returns ok, fint, [l1,l2,l3] monic, aux
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

// lambda for a single component: g irreducible of degree 1 or 2
// returns ok, lambda (squarefree rational), d (squarefree)
function CompLambda(uD, g)
  if Degree(g) eq 1 then
    r := -Coefficient(g,0)/LeadingCoefficient(g);
    gam := Evaluate(uD, r);
    if gam eq 0 then return false, 0, 0; end if;
    return true, SqFree(gam), 1;
  end if;
  K := NumberField(g);  th := K.1;
  gam := Evaluate(uD, th);
  if gam eq 0 then return false, 0, 0; end if;
  isq, ng := IsSquare(Norm(gam));
  if not isq then return false, 0, SqFree(Discriminant(g)); end if;
  del := 1 + gam/ng;
  if del eq 0 then ng := -ng; del := 1 + gam/ng; end if;
  lam := ng/Norm(del);
  return true, SqFree(lam), SqFree(Discriminant(g));
end function;

al := (mv^2-1)/(2*mv); tv := (mv^2+1)/(2*mv);
if nv eq 0 or mv in [Q|0,1,-1] then printf "BASE (%o,%o) degenerate\n", m, n; quit; end if;
sv := (nv^2 + al^4)/(2*nv);
Av := sv^2 - tv^4 + tv^2;
if Av eq 0 or sv eq 0 then printf "BASE (%o,%o) degenerate\n", m, n; quit; end if;

E0 := EllipticCurve([0, sv^2+Av, 0, sv^2*Av, 0]);
Em, mp := MinimalModel(E0);
r1, r2 := RankBounds(Em);
printf "BASE (%o,%o) s=%o t=%o rankbounds [%o,%o]\n", m, n, sv, tv, r1, r2;
if r2 eq 0 then printf "BASEDONE (%o,%o) rank0\n", m, n; quit; end if;
gens := [];
try
  G, mpG := MordellWeilGroup(Em);
  gens := [ mpG(G.i) : i in [1..Ngens(G)] | Order(G.i) eq 0 ];
catch e printf "BASE (%o,%o) MW failed\n", m, n; end try;
if #gens eq 0 then printf "BASEDONE (%o,%o) nogens\n", m, n; quit; end if;
tor2 := [ Em!0 ] cat [ Em![rr[1],0] : rr in Roots(DivisionPolynomial(Em,2)) ];
zs := {Q| };
for g in gens do
  for k in [1,-1,2,-2,3,-3] do
    for T in tor2 do
      Pt := k*g + T;
      if Pt eq Em!0 then continue; end if;
      P0 := Inverse(mp)(Pt);
      if P0[3] eq 0 then continue; end if;
      zz := P0[1];
      if zz ne 0 and zz ne -sv^2 and zz ne -Av then Include(~zs, zz); end if;
      if #zs ge NV then break g; end if;
    end for;
  end for;
end for;
printf "BASE (%o,%o) zvals=%o\n", m, n, #zs;

for zz in zs do
  vv := 1/zz;
  ok, fint, Ls, aux := MemberSTV(sv, tv, vv);
  if not ok then continue; end if;
  if #Sprint(Coefficients(fint)) gt 6000 then printf "MEMBER (%o,%o) v=%o TOOBIG\n", m,n,vv; continue; end if;
  try
    av := aux[1];
    // canonical component order: the two linear factors of l2, then l1, then l3
    facs := []; owner := [];
    okbc, rho := IsSquare(aux[2]*aux[3]);
    if not okbc then printf "MEMBER (%o,%o) v=%o BC_NOT_SQUARE\n", m,n,vv; continue; end if;
    Append(~facs, x - rho); Append(~owner, 2);
    Append(~facs, x + rho); Append(~owner, 2);
    for j in [1,3] do
      L := Ls[j];
      rts := Roots(L);
      if #rts eq 0 then Append(~facs, L); Append(~owner, j);
      else for rr in rts do Append(~facs, x - rr[1]); Append(~owner, j); end for; end if;
    end for;
    C := HyperellipticCurve(fint);  J := Jacobian(C);
    T, mT := TorsionSubgroup(J);
    inv := Invariants(T);
    D1 := SqFree(Discriminant(Ls[1]));
    A3 := SqFree(av);
    n4 := 0; nnorm := 0; npass := 0;
    lines := [];
    for gg in T do
      if Order(gg) ne 4 then continue; end if;
      n4 +:= 1;
      D := mT(gg); D2 := mT(2*gg);
      uD := D[1];
      if Degree(uD) ne D[3] or Degree(uD) ne 2 then Append(~lines, "  D4 DEGENERATE"); continue; end if;
      u2 := D2[1]/LeadingCoefficient(D2[1]);
      wj := 0; for j in [1..3] do if u2 eq Ls[j] then wj := j; end if; end for;
      nsq := []; resc := [];
      bad := false;
      for j in [1..3] do
        R := Resultant(Ls[j], uD);
        if R eq 0 then bad := true; break; end if;
        Append(~nsq, IsSquare(R)); Append(~resc, SqFree(R));
      end for;
      if bad then Append(~lines, "  D4 DEGENERATE"); continue; end if;
      pat := &cat[ (b select "1" else "0") : b in nsq ];
      if pat ne "111" then
        Append(~lines, Sprintf("  D4 above=T%o pat=%o resclass=%o", wj, pat, resc));
        continue;
      end if;
      nnorm +:= 1;
      lams := []; dds := []; okall := true;
      for i in [1..#facs] do
        okc, lam, dd := CompLambda(uD, facs[i]);
        if not okc then okall := false; break; end if;
        Append(~lams, lam); Append(~dds, dd);
      end for;
      if not okall then Append(~lines, "  D4 LAMBDA_FAIL"); continue; end if;
      // the two Q-components force lambda exactly
      lam2 := lams[1];
      consistent2 := (lams[1] eq lams[2]);
      def := []; okdef := true;
      for i in [3..#facs] do
        dfi := SqFree(lam2*lams[i]);
        Append(~def, dfi);
        if not (dfi eq 1 or dfi eq dds[i]) then okdef := false; end if;
      end for;
      cp := consistent2 and okdef;
      if cp then npass +:= 1; end if;
      Append(~lines, Sprintf("  D4 above=T%o pat=111 lam=%o d=%o lam2ok=%o def=%o compat=%o",
             wj, lams, dds, consistent2, def, cp));
    end for;
    printf "MEMBER (%o,%o) v=%o tors=%o D1=%o a=%o order4=%o norm111=%o xTpass=%o\n",
           m, n, vv, inv, D1, A3, n4, nnorm, npass;
    for l in lines do printf "%o\n", l; end for;
    if npass gt 0 then
      printf "LIFTCANDIDATE (%o,%o) v=%o tors=%o f=%o\n", m, n, vv, inv, fint;
    end if;
    n8 := #[ i : i in inv | i mod 8 eq 0 ];
    if n8 ge 2 then printf "JACKPOT_88 (%o,%o) v=%o tors=%o f=%o\n", m,n,vv,inv,fint;
    elif n8 eq 1 then printf "UPGRADE_8 (%o,%o) v=%o tors=%o f=%o\n", m,n,vv,inv,fint; end if;
  catch e
    printf "MEMBER (%o,%o) v=%o ERROR %o\n", m, n, vv, e`Object;
  end try;
end for;
printf "BASEDONE (%o,%o)\n", m, n;
quit;
