// claude_ov_88conic_base2.m -- Lane 5 (2026-07-25, resumed session).
//
// Finishes the conic-base sweep of the Nicholls Lambda_334 DOUBLE-STAGE-1 locus that
// claude_ov_88conic_base.m stalled on.  Two changes:
//
//   (1) SetClassGroupBounds("GRH") -- the fix Magma itself printed in
//       results/claude_ov_88conic_bases_partial.log.  Anything reached only through
//       the MordellWeilGroup route is therefore GRH-CONDITIONAL (route=MW below).
//   (2) An UNCONDITIONAL primary route: a naive point search on
//         E_{s,t} : y^2 = z(z+s^2)(z+A),   A = s^2 - t^4 + t^2
//       (route=PS below).  Every member found this way is an honest rational point,
//       so every MEMBER/D4 line with route=PS is unconditional; only the CLAIM
//       "these are all the members reachable from this base" would need GRH.
//
// Per base (m,n):  t=(m^2+1)/(2m), alpha=(m^2-1)/(2m), s=(n^2+alpha^4)/(2n).
// Per member v = 1/z:  build C1 : y^2 = d2*l1*l2*l3, factor type [1,1,2,2],
// exact TorsionSubgroup, then for each of the 12 order-4 elements D of Sigma=(Z/4)^2
//   * the norm presieve  Res(l_j^monic, u_D) square, j=1,2,3   (pattern "pat")
//   * the per-component rational representatives lambda_i (Hilbert 90)
//   * def_i := lambda_2 * lambda_i for the two quadratic components; the lift exists
//     iff def_1 in {1, disc l1} and def_3 in {1, a} modulo squares.
// D in 2J(Q) with ord(D)=4  <=>  Z/8 above D; [8,8]  <=>  all 12 pass.
//
// Params: BASEFILE (one "m n" per line), NCH (children, default 14),
//         PSB (point-search bound, default 3000), NV (max z per base, default 24),
//         USEMW (default 1), OUTDIR.
SetColumns(0);
SetMemoryLimit(4*10^9);
SetClassGroupBounds("GRH");
Q := Rationals();  P<x> := PolynomialRing(Q);  Z := Integers();

if not assigned BASEFILE then printf "need BASEFILE\n"; quit; end if;
if not assigned NCH then NCH := 14; elif Type(NCH) eq MonStgElt then NCH := StringToInteger(NCH); end if;
if not assigned PSB then PSB := 3000; elif Type(PSB) eq MonStgElt then PSB := StringToInteger(PSB); end if;
if not assigned NV then NV := 24; elif Type(NV) eq MonStgElt then NV := StringToInteger(NV); end if;
if not assigned USEMW then USEMW := 1; elif Type(USEMW) eq MonStgElt then USEMW := StringToInteger(USEMW); end if;
if not assigned OUTDIR then
  OUTDIR := "/tmp/claude-1000/-home-claude-torsion-jac/1c22f408-33bb-486e-870f-c5fca84bb522/scratchpad/conic2";
end if;
System("mkdir -p " cat OUTDIR);

function SqFree(r)
  if r eq 0 then return 0; end if;
  q := Numerator(r)*Denominator(r);
  s := Sign(q);  q := AbsoluteValue(q);
  res := s;
  for pf in Factorization(q) do if IsOdd(pf[2]) then res *:= pf[1]; end if; end for;
  return res;
end function;

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

// -------- one base: returns a list of output strings
function DoBase(mstr, nstr)
  out := [];
  mv := Q!eval(mstr);  nv := Q!eval(nstr);
  if nv eq 0 or mv in [Q|0,1,-1] then return [Sprintf("BASE (%o,%o) degenerate", mstr, nstr)]; end if;
  al := (mv^2-1)/(2*mv);  tv := (mv^2+1)/(2*mv);
  sv := (nv^2 + al^4)/(2*nv);
  Av := sv^2 - tv^4 + tv^2;
  if Av eq 0 or sv eq 0 then return [Sprintf("BASE (%o,%o) degenerate", mstr, nstr)]; end if;
  E0 := EllipticCurve([0, sv^2+Av, 0, sv^2*Av, 0]);
  Em, mp := MinimalModel(E0);
  route := "PS";
  gens := [];
  try
    pts := Points(Em : Bound := PSB);
    gens := [ p : p in pts | Order(p) eq 0 ];
  catch e ;
  end try;
  if #gens eq 0 and USEMW eq 1 then
    route := "MW";
    try
      G, mpG := MordellWeilGroup(Em);
      gens := [ mpG(G.i) : i in [1..Ngens(G)] | Order(G.i) eq 0 ];
    catch e
      Append(~out, Sprintf("BASE (%o,%o) MW_FAILED", mstr, nstr));
    end try;
  end if;
  Append(~out, Sprintf("BASE (%o,%o) s=%o t=%o route=%o gens=%o", mstr, nstr, sv, tv, route, #gens));
  if #gens eq 0 then Append(~out, Sprintf("BASEDONE (%o,%o) nogens", mstr, nstr)); return out; end if;
  // NB: Em![x,0] is WRONG whenever the minimal model has a1 or a3 nonzero -- it
  // raised "Point does not lie in scheme" on 14 of 8000 bases.  Take the 2-torsion
  // from the torsion subgroup instead.
  T2g, m2g := TorsionSubgroup(Em);
  tor2 := [ m2g(g) : g in T2g | 2*g eq T2g!0 ];
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
  Append(~out, Sprintf("BASE (%o,%o) zvals=%o", mstr, nstr, #zs));
  for zz in zs do
    vv := 1/zz;
    ok, fint, Ls, aux := MemberSTV(sv, tv, vv);
    if not ok then continue; end if;
    if #Sprint(Coefficients(fint)) gt 6000 then
      Append(~out, Sprintf("MEMBER (%o,%o) v=%o TOOBIG", mstr,nstr,vv)); continue; end if;
    try
      av := aux[1];
      okbc, rho := IsSquare(aux[2]*aux[3]);
      if not okbc then
        Append(~out, Sprintf("MEMBER (%o,%o) v=%o BC_NOT_SQUARE", mstr,nstr,vv)); continue; end if;
      facs := [x - rho, x + rho];  owner := [2,2];
      for j in [1,3] do
        L := Ls[j];  rts := Roots(L);
        if #rts eq 0 then Append(~facs, L); Append(~owner, j);
        else for rr in rts do Append(~facs, x - rr[1]); Append(~owner, j); end for; end if;
      end for;
      C := HyperellipticCurve(fint);  J := Jacobian(C);
      T, mT := TorsionSubgroup(J);
      inv := Invariants(T);
      D1 := SqFree(Discriminant(Ls[1]));  A3 := SqFree(av);
      n4 := 0; nnorm := 0; npass := 0;  lines := [];
      for gg in T do
        if Order(gg) ne 4 then continue; end if;
        n4 +:= 1;
        D := mT(gg); D2 := mT(2*gg);  uD := D[1];
        if Degree(uD) ne D[3] or Degree(uD) ne 2 then
          Append(~lines, "  D4 DEGENERATE"); continue; end if;
        u2 := D2[1]/LeadingCoefficient(D2[1]);
        wj := 0; for j in [1..3] do if u2 eq Ls[j] then wj := j; end if; end for;
        nsq := []; resc := []; bad := false;
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
        Append(~lines, Sprintf("  D4 above=T%o pat=111 d=%o lam2ok=%o def=%o compat=%o",
               wj, dds, consistent2, def, cp));
      end for;
      Append(~out, Sprintf("MEMBER (%o,%o) v=%o route=%o tors=%o D1=%o a=%o order4=%o norm111=%o xTpass=%o",
             mstr, nstr, vv, route, inv, D1, A3, n4, nnorm, npass));
      for l in lines do Append(~out, l); end for;
      if npass gt 0 then
        Append(~out, Sprintf("LIFTCANDIDATE (%o,%o) v=%o tors=%o f=%o", mstr, nstr, vv, inv, fint));
      end if;
      n8 := #[ i : i in inv | i mod 8 eq 0 ];
      if n8 ge 2 then Append(~out, Sprintf("JACKPOT_88 (%o,%o) v=%o tors=%o f=%o", mstr,nstr,vv,inv,fint));
      elif n8 eq 1 then Append(~out, Sprintf("UPGRADE_8 (%o,%o) v=%o tors=%o f=%o", mstr,nstr,vv,inv,fint)); end if;
    catch e
      Append(~out, Sprintf("MEMBER (%o,%o) v=%o ERROR %o", mstr, nstr, vv, e`Object));
    end try;
  end for;
  Append(~out, Sprintf("BASEDONE (%o,%o)", mstr, nstr));
  return out;
end function;

bases := [];
for ln in Split(Read(BASEFILE), "\n") do
  if ln eq "" or ln[1] eq "#" then continue; end if;
  tk := Split(ln, " ");
  if #tk ge 2 then Append(~bases, <tk[1], tk[2]>); end if;
end for;
printf "BASES_LOADED %o children=%o PSB=%o USEMW=%o\n", #bases, NCH, PSB, USEMW;

for c in [0..NCH-1] do
  pid := Fork();
  if pid eq 0 then
    F := Open(OUTDIR cat "/part" cat IntegerToString(c) cat ".txt", "w");
    i := c+1;
    while i le #bases do
      b := bases[i];
      try
        for l in DoBase(b[1], b[2]) do Puts(F, l); end for;
      catch e
        Puts(F, Sprintf("BASE (%o,%o) FATAL %o", b[1], b[2], e`Object));
      end try;
      Flush(F);
      i +:= NCH;
    end while;
    Puts(F, Sprintf("CHILD_DONE %o", c));
    Flush(F); delete F;
    quit;
  end if;
end for;
WaitForAllChildren();
nl := 0;
for c in [0..NCH-1] do
  for l in Split(Read(OUTDIR cat "/part" cat IntegerToString(c) cat ".txt"), "\n") do
    if #l gt 0 then printf "%o\n", l; nl +:= 1; end if;
  end for;
end for;
printf "CONIC2_DONE bases=%o lines=%o\n", #bases, nl;
quit;
