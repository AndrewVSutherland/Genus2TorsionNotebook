// claude_ov_88extend.m -- Lane 5 (2026-07-25, third session).
//
// UNCONDITIONAL extension of the DS1 sample for [8,8] on Nicholls' Lambda_334.
//
// The 108-base harvest cost 4 h 10 m and was killed at 25% because its generator
// supply came from MordellWeilGroup(E_{s,t}) -- heavy-tailed, and GRH-conditional.
// It also turned out to be badly redundant: 108 (m,n) bases collapse to only 22
// distinct (s^2,t^2) (the curve depends on (s^2,t^2,v) alone), and 332 member
// records to 68 distinct curves.
//
// This script drops MordellWeilGroup entirely.  Per base it does a naive point
// search on
//     E_{s,t} : y^2 = z(z+s^2)(z+A),   A = s^2 - t^4 + t^2,
// keeps the non-torsion points found, and builds up to NV members from them.  Every
// member so produced is an honest rational point of the family, and every statement
// about it (exact TorsionSubgroup, "D in 2*T", the x-T gate) is unconditional.  The
// only thing lost relative to the MW route is COMPLETENESS: bases whose generators
// have height above the search bound contribute nothing, and are logged as nogens.
//
// Bases are supplied one per DISTINCT (s^2,t^2) (data/claude_ov_88conic3_bases.txt),
// excluding the 22 already measured.
//
// Params: BASEFILE, NCH (children, default 8), PSB (point-search bound, default 800),
//         NV (max members per base, default 8), OUTDIR
SetColumns(0);
SetMemoryLimit(4*10^9);
Q := Rationals();  P<x> := PolynomialRing(Q);  Z := Integers();

if not assigned BASEFILE then printf "need BASEFILE\n"; quit; end if;
if not assigned NCH then NCH := 8; elif Type(NCH) eq MonStgElt then NCH := StringToInteger(NCH); end if;
if not assigned PSB then PSB := 800; elif Type(PSB) eq MonStgElt then PSB := StringToInteger(PSB); end if;
if not assigned NV then NV := 8; elif Type(NV) eq MonStgElt then NV := StringToInteger(NV); end if;
if not assigned OUTDIR then
  OUTDIR := "/tmp/claude-1000/-home-claude-torsion-jac/1c22f408-33bb-486e-870f-c5fca84bb522/scratchpad/extend";
end if;
System("mkdir -p " cat OUTDIR);

function SqFree(r)
  if r eq 0 then return 0; end if;
  m := Numerator(r)*Denominator(r);
  s := Sign(m);  m := AbsoluteValue(m);
  res := s;
  for pf in Factorization(m) do if IsOdd(pf[2]) then res *:= pf[1]; end if; end for;
  return res;
end function;

function MemberSTV(sv, tv, vv)
  if sv eq 0 or tv eq 0 or tv^2 eq 1 then return false, _; end if;
  Av := sv^2 - tv^4 + tv^2;
  if Av eq 0 then return false, _; end if;
  den := -sv^2*tv*Av*vv^2 + tv;
  if den eq 0 then return false, _; end if;
  uv := (-sv^2*Av*vv^2 - 2*Av*vv - 1)/den;
  if uv^2*sv^2 + 1 - tv^2 eq 0 then return false, _; end if;
  av := Av/(1 - tv^2);
  bv := Av/(uv^2*sv^2 + 1 - tv^2);
  cv := tv^2;
  d2 := Av*(sv^2*uv^2 + tv^4 - 2*tv^2 + 1)
          *(sv^4*uv^2 - sv^2*tv^2*uv^2 + sv^2*uv^2 - tv^6 + 3*tv^4 - 3*tv^2 + 1);
  if d2 eq 0 then return false, _; end if;
  l1 := (-av + bv + cv - 1)*x^2 + (2*av - 2*bv*cv)*x + (av*bv*cv - av*bv - av*cv + bv*cv);
  l2 := -x^2 + bv*cv;
  l3 := x^2 - av;
  f1 := l1*l2*l3;
  if Degree(f1) ne 6 then return false, _; end if;
  f := d2*f1/LeadingCoefficient(f1);
  if Discriminant(f) eq 0 then return false, _; end if;
  dd := LCM([Denominator(cc) : cc in Coefficients(f)]);
  fint := P![ cc*dd^2 : cc in Coefficients(f) ];
  cont := GCD([Z!cc : cc in Coefficients(fint)]);
  sqf := 1;
  for pf in Factorization(cont) do sqf *:= pf[1]^(2*(pf[2] div 2)); end for;
  fint := P![ (Z!cc) div sqf : cc in Coefficients(fint) ];
  return true, fint;
end function;

function Gate(uD, facs)
  lams := []; ds := [];
  for g in facs do
    if Degree(g) eq 1 then
      r := -Coefficient(g,0);
      gam := Evaluate(uD, r);
      if gam eq 0 then return "degenerate", [], [], false; end if;
      Append(~lams, SqFree(gam));  Append(~ds, 1);
    else
      K := NumberField(g);  th := K.1;
      gam := Evaluate(uD, th);
      if gam eq 0 then return "degenerate", [], [], false; end if;
      isq, ng := IsSquare(Norm(gam));
      if not isq then return "normfail", [], [], false; end if;
      del := 1 + gam/ng;
      if del eq 0 then ng := -ng; del := 1 + gam/ng; end if;
      lam := ng/Norm(del);
      if gam ne lam*del^2 then return "h90fail", [], [], false; end if;
      Append(~lams, SqFree(lam));  Append(~ds, SqFree(Discriminant(g)));
    end if;
  end for;
  k := #lams;
  for mask in [0 .. 2^k - 1] do
    cs := [ SqFree(lams[i]*(((mask div 2^(i-1)) mod 2) eq 0 select 1 else ds[i])) : i in [1..k] ];
    if #Seqset(cs) eq 1 then return "ok", lams, ds, true; end if;
  end for;
  return "ok", lams, ds, false;
end function;

function Analyse(tag, fint)
  out := [];
  C := HyperellipticCurve(fint);  J := Jacobian(C);
  T, mT := TorsionSubgroup(J);
  inv := Invariants(T);
  facs := [ g[1]/LeadingCoefficient(g[1]) : g in Factorization(fint) ];
  shape := Sort([ Degree(g) : g in facs ]);
  dbl := { mT(2*g) : g in T };
  n4 := 0; ntruth := 0; ngate := 0; nnorm := 0; ndegu := 0; ndis := 0;
  lines := [];
  for g in T do
    if Order(g) ne 4 then continue; end if;
    n4 +:= 1;
    D := mT(g);  uD := D[1];
    tr := D in dbl;
    if tr then ntruth +:= 1; end if;
    if Degree(uD) ne D[3] then ndegu +:= 1; end if;
    st, lams, ds, cp := Gate(uD, facs);
    if st eq "normfail" then nnorm +:= 1; end if;
    if cp then ngate +:= 1; end if;
    if cp ne tr then
      ndis +:= 1;
      Append(~lines, Sprintf("  DISAGREE %o u=%o st=%o gate=%o truth=%o lam=%o ds=%o", tag, uD, st, cp, tr, lams, ds));
    end if;
    if st eq "ok" and cp then
      Append(~lines, Sprintf("  D4PASS %o ds=%o lam=%o", tag, ds, lams));
    end if;
  end for;
  Append(~out, Sprintf("VERIFY %o deg=%o shape=%o tors=%o order4=%o normfail=%o degu=%o gate_pass=%o truth_pass=%o disagree=%o",
         tag, Degree(fint), shape, inv, n4, nnorm, ndegu, ngate, ntruth, ndis));
  for l in lines do Append(~out, l); end for;
  return out;
end function;

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
  gens := [];
  try
    pts := Points(Em : Bound := PSB);
    gens := [ p : p in pts | Order(p) eq 0 ];
  catch e ;
  end try;
  Append(~out, Sprintf("BASE (%o,%o) s=%o t=%o gens=%o", mstr, nstr, sv, tv, #gens));
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
    ok, fint := MemberSTV(sv, tv, vv);
    if not ok then continue; end if;
    if #Sprint(Coefficients(fint)) gt 6000 then
      Append(~out, Sprintf("VERIFY %o|%o|%o TOOBIG", sv, tv, vv)); continue; end if;
    try
      for l in Analyse(Sprintf("%o|%o|%o", sv, tv, vv), fint) do Append(~out, l); end for;
    catch e
      Append(~out, Sprintf("VERIFY %o|%o|%o ERROR %o", sv, tv, vv, e`Object));
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
printf "EXTEND_LOADED %o children=%o PSB=%o NV=%o\n", #bases, NCH, PSB, NV;

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
printf "EXTEND_DONE bases=%o lines=%o\n", #bases, nl;
quit;
