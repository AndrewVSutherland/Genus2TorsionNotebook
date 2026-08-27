// claude_88_worker.m — one (m,n) base of the overnight (8,8) hunt (driver: claude_88_drive.sh)
// Params: m, n (rationals as strings "a/b").  Prints BASE/MEMBER/GATEPASS/EXACT/JACKPOT_88/UPGRADE_48.
SetColumns(0);
SetMemoryLimit(2*10^9);
Q := Rationals(); Px<x> := PolynomialRing(Q);
mv := Q!eval(m);
nv := Q!eval(n);

function Member(mv, nv, vv)
  al := (mv^2-1)/(2*mv);
  tv := (mv^2+1)/(2*mv);
  sv := (nv^2+al^4)/(2*nv);
  Av := sv^2 - tv^4 + tv^2;
  if Av eq 0 or sv eq 0 or tv eq 0 or tv^2 eq 1 then return false, _; end if;
  den := -sv^2*tv*Av*vv^2 + tv;
  if den eq 0 then return false, _; end if;
  uv := (-sv^2*Av*vv^2 - 2*Av*vv - 1) / den;
  if uv^2*sv^2 + 1 - tv^2 eq 0 then return false, _; end if;
  av := Av/(1 - tv^2);
  bv := Av/(uv^2*sv^2 + 1 - tv^2);
  cv := tv^2;
  d2 := Av*(sv^2*uv^2 + tv^4 - 2*tv^2 + 1)*(sv^4*uv^2 - sv^2*tv^2*uv^2 + sv^2*uv^2 - tv^6 + 3*tv^4 - 3*tv^2 + 1);
  if d2 eq 0 then return false, _; end if;
  l1 := (-av + bv + cv - 1)*x^2 + (2*av - 2*bv*cv)*x + (av*bv*cv - av*bv - av*cv + bv*cv);
  l2 := -x^2 + bv*cv;
  l3 := x^2 - av;
  f1 := l1*l2*l3;
  if Degree(f1) ne 6 then return false, _; end if;
  f := d2*f1/LeadingCoefficient(f1);
  if Discriminant(f) eq 0 then return false, _; end if;
  dd := LCM([Denominator(cc) : cc in Coefficients(f)]);
  fint := Px![ cc*dd^2 : cc in Coefficients(f) ];
  cont := GCD([Integers()!cc : cc in Coefficients(fint)]);
  sqf := 1;
  for pf in Factorization(cont) do sqf *:= pf[1]^(2*(pf[2] div 2)); end for;
  return true, Px![ (Integers()!cc) div sqf : cc in Coefficients(fint) ];
end function;

function GatePass(fint, np)
  good := 0;
  for p in PrimesInInterval(11, 200) do
    if good ge np then break; end if;
    fp := PolynomialRing(GF(p))!fint;
    if Degree(fp) lt 5 or Discriminant(fp) eq 0 then continue; end if;
    try
      G := AbelianGroup(Jacobian(HyperellipticCurve(fp)));
      if #[ i : i in Invariants(G) | i mod 8 eq 0 ] lt 2 then return false; end if;
      good +:= 1;
    catch e
      good := good;
    end try;
  end for;
  return good ge np;
end function;

al := (mv^2-1)/(2*mv); tv := (mv^2+1)/(2*mv); sv := (nv^2+al^4)/(2*nv);
Av := sv^2 - tv^4 + tv^2;
if Av eq 0 or sv eq 0 then printf "BASE (%o,%o) degenerate\n", m, n; quit; end if;
E0 := EllipticCurve([0, sv^2+Av, 0, sv^2*Av, 0]);
Em, mp := MinimalModel(E0);
r1, r2 := RankBounds(Em);
printf "BASE (%o,%o) rankbounds [%o,%o]\n", m, n, r1, r2;
if r2 eq 0 then quit; end if;
gens := [];
try
  G, mpG := MordellWeilGroup(Em);
  gens := [ mpG(G.i) : i in [1..Ngens(G)] | Order(G.i) eq 0 ];
catch e printf "BASE (%o,%o) MW failed\n", m, n; end try;
if #gens eq 0 then quit; end if;
tor2 := [ Em!0 ] cat [ T : T in [Em![xx,0] : xx in [r[1] : r in Roots(DivisionPolynomial(Em,2))]] ];
zs := {Q| };
for g in gens do
  for k in [1,-1,2,-2,3] do
    for T in tor2 do
      P := k*g + T;
      if P eq Em!0 then continue; end if;
      P0 := Inverse(mp)(P);
      if P0[3] eq 0 then continue; end if;
      zz := P0[1];
      if zz ne 0 and zz ne -sv^2 and zz ne -Av then Include(~zs, zz); end if;
      if #zs ge 24 then break g; end if;
    end for;
  end for;
end for;
printf "BASE (%o,%o) zvals %o\n", m, n, #zs;
for zz in zs do
  ok, fint := Member(mv, nv, 1/zz);
  if not ok then continue; end if;
  if #Sprint(Coefficients(fint)) gt 6000 then continue; end if;
  printf "MEMBER (%o,%o) v=%o\n", m, n, 1/zz;
  if not GatePass(fint, 6) then continue; end if;
  printf "GATEPASS (%o,%o) v=%o\n", m, n, 1/zz;
  Cs := SimplifiedModel(ReducedMinimalWeierstrassModel(HyperellipticCurve(fint)));
  T := TorsionSubgroup(Jacobian(Cs));
  inv := Invariants(T);
  printf "EXACT (%o,%o) v=%o : J1 torsion %o\n", m, n, 1/zz, inv;
  n8 := #[ i : i in inv | i mod 8 eq 0 ];
  if n8 ge 2 then
    printf "JACKPOT_88 (%o,%o) v=%o torsion=%o f=%o\n", m, n, 1/zz, inv, fint;
  elif n8 eq 1 and #[ i : i in inv | i mod 4 eq 0 ] ge 2 then
    printf "UPGRADE_48 (%o,%o) v=%o torsion=%o f=%o\n", m, n, 1/zz, inv, fint;
  end if;
end for;
printf "BASEDONE (%o,%o)\n", m, n;
quit;
