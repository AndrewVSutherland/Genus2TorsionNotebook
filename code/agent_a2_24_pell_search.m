
//////////////////////////////////////////////////////////////////////
//  Z/2 x Z/24 via the Platonov-Petrunin fundamental-unit (Pell) method.
//
//  For a monic squarefree sextic f, the infinity divisor
//  D = infty_+ - infty_- has finite order detected by a rank drop of the
//  Hankel matrix of the Laurent coefficients of sqrt(f) (Platonov-
//  Petrunin; verified in fundamental_unit_rank_test.m).  A unit of
//  degree m signals order-m structure on Jac(y^2=f).
//
//  Scan integer monic sextics; require the mod-p Hankel rank drop for the
//  target unit degree at all filter primes (cheap necessary condition);
//  survivors get an exact TorsionSubgroup + 2-rank + simplicity check.
//  Flag TARGET_2_24 for torsion containing Z/2 x Z/24.
//
//  Usage:
//    magma -b Validate:=true agent_a2_24_pell_search.m
//    magma -b H:=3 UnitDeg:=24 NParts:=3 Part:=0 agent_a2_24_pell_search.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemGB then MemGB := 3; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

if not assigned H then H := 3; elif Type(H) eq MonStgElt then H := StringToInteger(H); end if;
if not assigned UnitDeg then UnitDeg := 24; elif Type(UnitDeg) eq MonStgElt then UnitDeg := StringToInteger(UnitDeg); end if;
if not assigned NParts then NParts := 1; elif Type(NParts) eq MonStgElt then NParts := StringToInteger(NParts); end if;
if not assigned Part then Part := 0; elif Type(Part) eq MonStgElt then Part := StringToInteger(Part); end if;
if not assigned Validate then Validate := false;
elif Type(Validate) eq MonStgElt then Validate := Validate in {"true","True","1","yes"}; end if;
if not assigned progress then progress := 200000;
elif Type(progress) eq MonStgElt then progress := StringToInteger(progress); end if;

ModPrimes := [5,7,11,13,17,19,23,29];

// ---- modular Hankel rank-drop filter ----
function SqrtCoeffsMod(coeffs, nmax, p)
    F := GF(p);
    c := [F!1];
    b := AssociativeArray(Integers());
    for n in [1..6] do b[n] := F!coeffs[n]; end for;
    for n in [1..nmax] do
        rhs := IsDefined(b, n) select b[n] else F!0;
        conv := (n eq 1) select F!0 else &+[c[i+1]*c[n-i+1] : i in [1..n-1]];
        Append(~c, (rhs - conv)/(F!2));
    end for;
    return c;
end function;

// unit of degree m <=> Hankel H_r (r = m-2) drops rank; d_{-k} = c[k+3+1]
function ModularRankDrops(coeffs, m, p)
    r := m - 2;
    c := SqrtCoeffsMod(coeffs, 2*r + 6, p);
    Hm := Matrix(GF(p), [[c[(i+j+1)+3+1] : j in [0..r-1]] : i in [0..r]]);
    return Rank(Hm) lt r;
end function;

function PassesModular(coeffs, m)
    for p in ModPrimes do
        // skip primes of bad reduction (2 | leading? monic ok; disc check cheap-ish skip)
        if not ModularRankDrops(coeffs, m, p) then return false; end if;
    end for;
    return true;
end function;

// ---- exact Laurent/Hankel over Q (from fundamental_unit_rank_test.m) ----
function SqrtLaurent(f, nmax)
    b := AssociativeArray(Integers());
    for n in [1..6] do b[n] := Q!Coefficient(f, 6-n); end for;
    c := [Q!1];
    for n in [1..nmax] do
        rhs := IsDefined(b,n) select b[n] else Q!0;
        conv := (n eq 1) select Q!0 else &+[c[i+1]*c[n-i+1] : i in [1..n-1]];
        Append(~c, (rhs - conv)/2);
    end for;
    return c;
end function;

function HasUnitDegree(f, m)
    r := m - 2;
    if r le 0 then return false; end if;
    c := SqrtLaurent(f, 2*r + 6);
    Hm := Matrix(Q, [[c[(i+j+1)+3+1] : j in [0..r-1]] : i in [0..r]]);
    if Rank(Hm) ge r then return false; end if;
    N := Nullspace(Transpose(Hm));
    for v in Basis(N) do
        if v[r] ne 0 then
            beta := &+[v[j+1]*x^j : j in [0..r-1]];
            cc := SqrtLaurent(f, Degree(beta)+8);
            alpha := &+[ (&+[Coefficient(beta,j)*cc[3-(ell-j)+1] : j in [0..Degree(beta)] | ell-j le 3]) * x^ell
                        : ell in [0..Degree(beta)+3] ];
            if Degree(alpha^2 - beta^2*f) le 0 then return true; end if;
        end if;
    end for;
    return false;
end function;

function TwoRank(fp)
    degs := [Degree(g[1]) : g in Factorization(fp)];
    k := #degs; even := 0;
    for mask in [0..2^k-1] do
        ss := 0; for i in [1..k] do if (mask div 2^(i-1)) mod 2 eq 1 then ss +:= degs[i]; end if; end for;
        if ss mod 2 eq 0 then even +:= 1; end if;
    end for;
    return Ilog2(even) - 1;
end function;

function CountCurve(fp)
    Fq := BaseRing(Parent(fp)); cnt := 0;
    for xx in Fq do vv := Evaluate(fp, xx);
        if vv eq 0 then cnt +:= 1; elif IsSquare(vv) then cnt +:= 2; end if; end for;
    if IsSquare(LeadingCoefficient(fp)) then cnt +:= 2; end if;
    return cnt;
end function;

function SimplicityCertificate(fInt)
    RT := PolynomialRing(Q); T := RT.1;
    dsc := Discriminant(fInt);
    for pp in [3,5,7,11,13,17,19,23,29,31,37,41,43,47] do
        if (Z!LeadingCoefficient(fInt)) mod pp eq 0 then continue; end if;
        if (Z!Numerator(dsc)) mod pp eq 0 then continue; end if;
        PF := PolynomialRing(GF(pp));
        fp := PF![GF(pp)!co : co in Coefficients(fInt)];
        if Degree(fp) lt 5 or not IsSquarefree(fp) then continue; end if;
        PF2 := PolynomialRing(GF(pp^2));
        fp2 := PF2![GF(pp^2)!co : co in Coefficients(fInt)];
        a1 := pp + 1 - CountCurve(fp);
        a2 := (CountCurve(fp2) - pp^2 - 1 + a1^2) div 2;
        chi := T^4 - a1*T^3 + a2*T^2 - a1*pp*T + pp^2;
        if not IsIrreducible(chi) then continue; end if;
        K := NumberField(chi); pi := K.1; drop := false;
        for nn in [2..12] do
            if Degree(MinimalPolynomial(pi^nn)) lt 4 then drop := true; break; end if;
        end for;
        if not drop then return true, pp, chi; end if;
    end for;
    return false, 0, RT!0;
end function;

function SexticFrom(cf)
    return x^6 + cf[1]*x^5 + cf[2]*x^4 + cf[3]*x^3 + cf[4]*x^2 + cf[5]*x + cf[6];
end function;

if Validate then
    print "== VALIDATION: filter recovers Platonov-Petrunin unit degrees ==";
    tests := [
        <"f14", [Q|0,2,4,6,4,4], 14>,     // (x^2+1)(x^4+5x^2+4x+4) expanded
        <"f18", [Q|-2,11,-2,-1,0,-8], 18> // placeholder; recompute below
    ];
    // recompute expansions exactly
    f14 := (x^2+1)*(x^4+5*x^2+4*x+4);
    f18 := (x^2-x+1)*(x^4-x^3+9*x^2+8*x-8);
    for pr in [<"f14",f14,14>, <"f18",f18,18>] do
        cf := [Coefficient(pr[2], 6-i) : i in [1..6]];
        okmod := PassesModular(cf, pr[3]);
        okex := HasUnitDegree(pr[2], pr[3]);
        printf "  %o: modular-pass(m=%o)=%o exact=%o (expect true)\n", pr[1], pr[3], okmod, okex;
    end for;
    print "DONE"; quit;
end if;

printf "PELL SEARCH H=%o UnitDeg=%o Part=%o/%o filterprimes=%o\n", H, UnitDeg, Part, NParts, ModPrimes;
tested := 0; modpass := 0; unitok := 0; hits := 0;
a5idx := 0;
for a5 in [-H..H] do
  a5idx +:= 1;
  if (a5idx mod NParts) ne Part then continue; end if;
  for a4 in [-H..H] do for a3 in [-H..H] do for a2 in [-H..H] do
    for a1 in [-H..H] do for a0 in [-H..H] do
      if a0 eq 0 then continue; end if;
      tested +:= 1;
      if tested mod progress eq 0 then
        printf "PROGRESS tested=%o modpass=%o unitok=%o hits=%o\n", tested, modpass, unitok, hits;
      end if;
      // skip even sextics (a5=a3=a1=0): bielliptic => split, never simple
      if a5 eq 0 and a3 eq 0 and a1 eq 0 then continue; end if;
      cf := [Q|a5,a4,a3,a2,a1,a0];
      if not PassesModular(cf, UnitDeg) then continue; end if;
      modpass +:= 1;
      f := SexticFrom(cf);
      if not IsSquarefree(f) then continue; end if;
      if not HasUnitDegree(f, UnitDeg) then continue; end if;
      unitok +:= 1;
      // exact torsion + 2-rank
      J := Jacobian(HyperellipticCurve(f));
      inv := Invariants(TorsionSubgroup(J));
      r2 := TwoRank(f);
      ord := #inv eq 0 select 1 else &*inv;
      has24 := #[n : n in inv | n mod 24 eq 0] ge 1;
      has2 := #[n : n in inv | n mod 2 eq 0] ge 2;
      // only surface interesting cases: order divisible by 8, or rank>=2
      if (ord mod 8 eq 0) or r2 ge 2 then
        printf "UNIT_HIT coeffs=%o torsion=%o 2rank=%o\n", cf, inv, r2;
      end if;
      if has24 and has2 then
        hits +:= 1;
        issimple, pp, chi := SimplicityCertificate(f);
        printf "TARGET_2_24 coeffs=%o torsion=%o simple=%o%o\n  f=%o\n",
            cf, inv, issimple, issimple select Sprintf(" (p=%o chi=%o)",pp,chi) else "", f;
      end if;
    end for; end for;
  end for; end for; end for;
end for;
printf "SEARCH_DONE tested=%o modpass=%o unitok=%o hits=%o\n", tested, modpass, unitok, hits;
print "DONE";
quit;
