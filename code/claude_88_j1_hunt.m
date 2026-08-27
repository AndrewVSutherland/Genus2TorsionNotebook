// claude_88_j1_hunt.m — overnight (8,8) hunter on the Nicholls Lambda_334 substrate.
// For bases (m,n) of increasing height: parametrize (s,t), find rational points on the
// double-stage-1 elliptic fibration E_{s,t}: w^2 = z(z+s^2)(z+A) by direct point search,
// map each to a member (s,t,v=1/z) — on these, BOTH prescribed classes are 2-divisible
// in J2(Q) (necessary for (8,8) on J1).  Then test J1 = Jac(y^2 = d2*f1) for torsion
// beyond [4,4]: gate by (Z/8)^2 embedding in J1(F_p) at several primes, then exact
// TorsionSubgroup.  JACKPOT = two invariant factors divisible by 8 (a NEW group [8,8]).
// Params: HMAX (base height bound, default 8), PSBOUND (point search, default 10^7).
SetColumns(0);
SetMemoryLimit(8*10^9);
if not assigned HMAX then HMAX := 8; elif Type(HMAX) eq MonStgElt then HMAX := StringToInteger(HMAX); end if;
if not assigned PSBOUND then PSBOUND := 10000000; elif Type(PSBOUND) eq MonStgElt then PSBOUND := StringToInteger(PSBOUND); end if;

Q := Rationals(); Px<x> := PolynomialRing(Q);

function Member(mv, nv, vv)
  // returns ok, integral model of C1: y^2 = d2*f1 (squarefree-scaled), plus (s,t)
  al := (mv^2-1)/(2*mv);
  tv := (mv^2+1)/(2*mv);
  sv := (nv^2+al^4)/(2*nv);
  Av := sv^2 - tv^4 + tv^2;
  if Av eq 0 or sv eq 0 or tv eq 0 or tv^2 eq 1 then return false, _, _, _; end if;
  den := -sv^2*tv*Av*vv^2 + tv;
  if den eq 0 then return false, _, _, _; end if;
  uv := (-sv^2*Av*vv^2 - 2*Av*vv - 1) / den;
  if uv^2*sv^2 + 1 - tv^2 eq 0 then return false, _, _, _; end if;
  av := Av/(1 - tv^2);
  bv := Av/(uv^2*sv^2 + 1 - tv^2);
  cv := tv^2;
  d2 := Av*(sv^2*uv^2 + tv^4 - 2*tv^2 + 1)*(sv^4*uv^2 - sv^2*tv^2*uv^2 + sv^2*uv^2 - tv^6 + 3*tv^4 - 3*tv^2 + 1);
  if d2 eq 0 then return false, _, _, _; end if;
  l1 := (-av + bv + cv - 1)*x^2 + (2*av - 2*bv*cv)*x + (av*bv*cv - av*bv - av*cv + bv*cv);
  l2 := -x^2 + bv*cv;
  l3 := x^2 - av;
  f1 := l1*l2*l3;
  lc := LeadingCoefficient(f1);
  if lc eq 0 then return false, _, _, _; end if;
  f := d2/lc * f1;
  if Degree(f) ne 6 or Discriminant(f) eq 0 then return false, _, _, _; end if;
  // integral squarefree-ish scaling: multiply by den^2
  dd := LCM([Denominator(cc) : cc in Coefficients(f)]);
  fint := Px![ cc*dd^2 : cc in Coefficients(f) ];
  // strip square content to keep coefficients small
  cont := GCD([Integers()!cc : cc in Coefficients(fint)]);
  sqf := 1;
  for pf in Factorization(cont) do sqf *:= pf[1]^(2*(pf[2] div 2)); end for;
  fint := Px![ cc div sqf : cc in Coefficients(fint) ];
  return true, fint, sv, tv;
end function;

function GatePass(fint, np)
  good := 0;
  for p in PrimesInInterval(11, 200) do
    if good ge np then break; end if;
    Fp := GF(p);
    fp := PolynomialRing(Fp)!fint;
    if Degree(fp) lt 5 or Discriminant(fp) eq 0 then continue; end if;
    ok := true;
    try
      Jp := Jacobian(HyperellipticCurve(fp));
      G := AbelianGroup(Jp);
      inv := Invariants(G);
      n8 := #[ i : i in inv | i mod 8 eq 0 ];
      if n8 lt 2 then return false, p; end if;
    catch e ok := false; end try;
    if ok then good +:= 1; end if;
  end for;
  return good ge np, 0;
end function;

nb := 0; nmem := 0; ngate := 0; nexact := 0;
hts := [];
for hm in [2..HMAX] do for hn in [1..HMAX] do Append(~hts, <hm+hn, hm, hn>); end for; end for;
Sort(~hts, func<A, B | A[1] - B[1]>);

for hh in hts do
  hm := hh[2]; hn := hh[3];
  ms := [Q| ]; ns := [Q| ];
  for aa in [-hm..hm] do
    if GCD(aa, hm) eq 1 and aa ne 0 then Append(~ms, Q!aa/hm); if hm gt 1 then Append(~ms, Q!hm/aa); end if; end if;
  end for;
  for aa in [-hn..hn] do
    if GCD(aa, hn) eq 1 and aa ne 0 then Append(~ns, Q!aa/hn); if hn gt 1 then Append(~ns, Q!hn/aa); end if; end if;
  end for;
  for mv in ms do
    if mv in [Q|0,1,-1] then continue; end if;
    for nv in ns do
      if nv eq 0 then continue; end if;
      nb +:= 1;
      try
        al := (mv^2-1)/(2*mv); tv := (mv^2+1)/(2*mv); sv := (nv^2+al^4)/(2*nv);
        Av := sv^2 - tv^4 + tv^2;
        if Av eq 0 or sv eq 0 then continue; end if;
        E0 := EllipticCurve([0, sv^2+Av, 0, sv^2*Av, 0]);
        Em, mp := MinimalModel(E0);
        pts := Points(Em : Bound := PSBOUND);
        zs := {Q| };
        for P in pts do
          P0 := Inverse(mp)(P);
          zz := P0[1];
          if zz ne 0 and zz ne -sv^2 and zz ne -Av and P0[3] ne 0 then Include(~zs, zz); end if;
        end for;
        if #zs eq 0 then continue; end if;
        for zz in zs do
          vv := 1/zz;
          ok, fint, _, _ := Member(mv, nv, vv);
          if not ok then continue; end if;
          if #Sprint(Coefficients(fint)) gt 4000 then continue; end if;  // size guard
          nmem +:= 1;
          gp, killp := GatePass(fint, 6);
          if not gp then continue; end if;
          ngate +:= 1;
          printf "GATEPASS (m,n,v)=(%o,%o,%o)\n", mv, nv, vv;
          Cs := SimplifiedModel(ReducedMinimalWeierstrassModel(HyperellipticCurve(fint)));
          T := TorsionSubgroup(Jacobian(Cs));
          inv := Invariants(T);
          nexact +:= 1;
          printf "EXACT (m,n,v)=(%o,%o,%o): J1 torsion %o\n", mv, nv, vv, inv;
          n8 := #[ i : i in inv | i mod 8 eq 0 ];
          if n8 ge 2 then
            printf "JACKPOT_88 (m,n,v)=(%o,%o,%o) torsion=%o f=%o\n", mv, nv, vv, inv, fint;
          elif n8 eq 1 and #[ i : i in inv | i mod 4 eq 0 ] ge 2 then
            printf "UPGRADE_48 (m,n,v)=(%o,%o,%o) torsion=%o f=%o\n", mv, nv, vv, inv, fint;
          end if;
        end for;
      catch e
        printf "BASEERR (m,n)=(%o,%o): %o\n", mv, nv, e`Object;
      end try;
    end for;
  end for;
  printf "PROGRESS height=%o bases=%o members=%o gatepass=%o exact=%o\n", hh[1], nb, nmem, ngate, nexact;
end for;
printf "SEARCH_DONE bases=%o members=%o gatepass=%o exact=%o\n", nb, nmem, ngate, nexact;
quit;
