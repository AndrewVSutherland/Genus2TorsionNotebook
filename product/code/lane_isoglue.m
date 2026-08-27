// lane_isoglue.m — (N,N)-gluings of ISOGENOUS pairs.  For an m-isogeny
// phi: E -> F with gcd(m,N)=1, psi = c*phi|E[N] is an anti-isometry iff
// c^2*m = -1 mod N, i.e. m mod 5 in {1,4} (N=5), m mod 7 in {3,5,6} (N=7).
// The Kani product-degeneracy would need N*(a^2+m*b^2) to be a perfect
// square for a graph-compatible (a,b) (b = c*a mod N) -- checked per pair;
// the analytic gluer's |I10|>0 guard also rejects degenerate output, so a
// surviving curve is a genuine Jacobian and T1 x T2 injects
// (gcd(N,#T1#T2)=1).  Targets (product not in KNOWN):
//   [2,6]x[2,6]  m=+-1(5): [2,2,6,6] 144      [12]x[12]: [12,12] 144
//   [2,6]x[12]:  [2,6,12] 144                 [2,8]x[2,8]: [2,2,8,8] 256
//   [2,8]x[2,6]: [2,2,2,24] 192               [2,8]x[12]: [2,4,24] 192  etc.
// Usage: cd product/code && magma -b lane_isoglue.m > ../logs/lane_isoglue.log
//   optional: H:=<int> (default 25), MaxGlue:=<int> (default 12)
SetColumns(0);
if not assigned H then H := 25; elif Type(H) eq MonStgElt then H := StringToInteger(H); end if;
if not assigned MaxGlue then MaxGlue := 12; elif Type(MaxGlue) eq MonStgElt then MaxGlue := StringToInteger(MaxGlue); end if;
SetMemoryLimit(8*10^9);
load "split_lab.m";
load "analytic_glue.m";

RQ := Rationals();

// families (as in lane_cong_sieve.m)
function E26of(tv)
    if tv in {RQ|3,-3,1,5,9} then return false, 0; end if;
    r1 := (-2*tv+10)/((tv+3)*(tv-3));
    r2 := (-tv^3+7*tv^2-11*tv+5)/(4*(tv+3)*(tv-3)^2);
    r3 := (-2*tv^2+4*tv-2)/((tv+3)^2*(tv-3));
    if r1 eq r2 or r1 eq r3 or r2 eq r3 then return false, 0; end if;
    return true, EllipticCurve((RQx.1-r1)*(RQx.1-r2)*(RQx.1-r3));
end function;
function E28of(tv)
    if tv eq 0 or tv eq -1/4 or 8*tv^2-1 eq 0 then return false, 0; end if;
    r1 := (16*tv^3+12*tv^2+2*tv)/(8*tv^2-1)^2;
    r2 := (32*tv^3+24*tv^2+8*tv+1)/(16*tv^2*(8*tv^2-1));
    r3 := (-32*tv^4-32*tv^3-12*tv^2-2*tv)/((4*tv+1)^2*(8*tv^2-1));
    if r1 eq r2 or r1 eq r3 or r2 eq r3 then return false, 0; end if;
    return true, EllipticCurve((RQx.1-r1)*(RQx.1-r2)*(RQx.1-r3));
end function;
function X1N(N, tv)
    ok := true; bp := RQ!0; cp := RQ!0;
    try
        case N:
        when 10: bp := tv^3*(tv-1)*(2*tv-1)/(tv^2-3*tv+1)^2; cp := -tv*(tv-1)*(2*tv-1)/(tv^2-3*tv+1);
        when 12: bp := tv*(2*tv-1)*(3*tv^2-3*tv+1)*(2*tv^2-2*tv+1)/(tv-1)^4; cp := -tv*(2*tv-1)*(3*tv^2-3*tv+1)/(tv-1)^3;
        end case;
    catch e
        ok := false;
    end try;
    if not ok or bp eq 0 then return false, 0; end if;
    E := 0;
    try E := EllipticCurve([1-cp, -bp, -bp, 0, 0]); catch e return false, 0; end try;
    return true, E;
end function;
function HeightRats(H)
    S := [];
    for a in [1..H], b in [1..H] do
        if GCD(a,b) eq 1 then Append(~S, a/b); Append(~S, -a/b); end if;
    end for;
    return S;
end function;

// pool of high-torsion curves (minimal models, deduped by j)
pool := [* *];
jseen := {};
for fam in [* <"26",0>, <"28",0>, <"x10",10>, <"x12",12> *] do
    cnt := 0;
    for tv in HeightRats(H) do
        ok := false; E := 0;
        if fam[1] eq "26" then ok, E := E26of(tv);
        elif fam[1] eq "28" then ok, E := E28of(tv);
        else ok, E := X1N(fam[2], tv);
        end if;
        if not ok then continue; end if;
        j := jInvariant(E);
        if j in jseen then continue; end if;
        Include(~jseen, j);
        Ti := Invariants(TorsionSubgroup(E));
        n1 := IsEmpty(Ti) select 1 else &*Ti;
        if n1 lt 10 then continue; end if;   // want >= [10]/[12]/[2,6]/[2,8]
        Append(~pool, < fam[1], tv, MinimalModel(E), Ti, j, n1 >);
        cnt +:= 1;
    end for;
    printf "POOL %o: %o curves\n", fam[1], cnt;
end for;
printf "POOL total %o\n", #pool;

// candidate (E,F,N,m) with anti-isometry available and interesting product
OKM := AssociativeArray();
OKM[5] := {1,4};  OKM[7] := {3,5,6};

function KaniDegenerate(N, m, cc)
    // graph-compatible principal subcurve exists iff N*(a^2+m*b^2) is a
    // perfect square with b = cc*a mod N (a,b) integers not both 0, small range
    for a in [-4*N..4*N] do
        for b in [-4*N..4*N] do
            if a eq 0 and b eq 0 then continue; end if;
            if (b - cc*a) mod N ne 0 then continue; end if;
            v := N*(a^2 + m*b^2);
            if IsSquare(v) then return true, <a,b>; end if;
        end for;
    end for;
    return false, <0,0>;
end function;

cands := [* *];
seencls := {};
nscan := 0;
for i in [1..#pool] do
    rec := pool[i];
    E := rec[3];
    NE := Conductor(E);
    cls := [* *];
    ok := true;
    try
        cls := IsogenousCurves(E);
    catch e
        ok := false;
    end try;
    if not ok or #cls eq 0 then continue; end if;
    clskey := Sprintf("%o", Sort([ jInvariant(cls[k]) : k in [1..#cls] ]));
    if clskey in seencls then continue; end if;   // one scan per isogeny class
    Include(~seencls, clskey);
    nscan +:= 1;
    // collect torsion of each member
    Ts := [ Invariants(TorsionSubgroup(cls[k])) : k in [1..#cls] ];
    for k1 in [1..#cls] do
        for k2 in [k1+1..#cls] do
            T1 := Ts[k1]; T2 := Ts[k2];
            n12 := (IsEmpty(T1) select 1 else &*T1) * (IsEmpty(T2) select 1 else &*T2);
            if n12 lt 100 then continue; end if;
            prod := Invariants(AbelianGroup(T1 cat T2));
            if prod in KNOWN then continue; end if;
            m := 0;
            try okm, mm := IsIsogenous(cls[k1], cls[k2]); if okm and Type(mm) eq RngIntElt then m := mm; end if; catch e; end try;
            for N in [5,7] do
                if GCD(N, n12) ne 1 then continue; end if;
                if m ne 0 and not (m mod N) in OKM[N] then continue; end if;
                Append(~cands, < cls[k1], cls[k2], N, m, T1, T2, prod, Conductor(cls[k1]) >);
            end for;
        end for;
    end for;
end for;
printf "ISOCANDS %o (classes scanned %o)\n", #cands, nscan;

// sort by product order (descending), then by conductor (ascending: small
// conductor = small invariants = fast Mestre), and glue the top ones
ords := [ &*c[7] : c in cands ];
ord_idx := Reverse(Sort([ <ords[i], -cands[i][8], i> : i in [1..#cands] ]));
nglue := 0;
for oi in ord_idx do
    c := cands[oi[3]];
    if nglue ge MaxGlue then break; end if;
    nglue +:= 1;
    E := c[1]; F := c[2]; N := c[3];
    printf "==== GLUING N=%o: %o x %o -> hope %o (cond %o) ====\n", N, c[5], c[6], c[7], c[8];
    if c[4] ne 0 then
        for cc in [1..N-1] do
            if (cc^2*c[4] + 1) mod N eq 0 then
                dg, ab := KaniDegenerate(N, c[4], cc);
                printf "  c=%o Kani-degenerate=%o %o\n", cc, dg, ab;
            end if;
        end for;
    end if;
    outs := GlueCandidates(E, F, N : prec := 200);
    for o in outs do
        C := o[1];
        okI, g := IntegralSextic(C);
        if not okI then printf "GLUEOUT: bad model\n"; continue; end if;
        // certify Lp split
        dz := Integers()!Discriminant(g);
        lc := Integers()!LeadingCoefficient(g);
        bad := 2*dz*lc*Conductor(E)*Conductor(F);
        nok := 0; nbad := 0;
        for p in PrimesUpTo(80) do
            if bad mod p eq 0 then continue; end if;
            gp := PolynomialRing(GF(p))!g;
            Lp := Numerator(ZetaFunction(HyperellipticCurve(gp)));
            PZ<T> := Parent(Lp);
            LE := 1 - TraceOfFrobenius(E,p)*T + p*T^2;
            LF := 1 - TraceOfFrobenius(F,p)*T + p*T^2;
            if Lp eq LE*LF then nok +:= 1; else nbad +:= 1; end if;
        end for;
        printf "GLUEOUT M=%o twist=%o Lpsplit=%o/%o\n", o[2], o[3], nok, nok+nbad;
        if nbad gt 0 or nok lt 4 then continue; end if;
        einv := [Integers()|-1];
        try einv := ExactTorsion(g); catch e printf "GLUEOUT: torsion failed\n"; continue; end try;
        isnew := not einv in KNOWN;
        printf "ISOGLUE_RESULT invs=%o %o (hoped %o) g=%o\n",
            einv, isnew select "*** NEW GROUP ***" else "known", c[7], g;
    end for;
end for;
printf "ISOGLUE_DONE glued=%o\n", nglue;
quit;
