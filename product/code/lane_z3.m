// lane_z3.m — sweep Frengley's Z(3,2) surface (anti-3-congruent pairs of
// elliptic curves, the gluable power for (3,3)-gluing) and glue pairs with
// 3-coprime torsion via Genus2Elliptic3 — the image is then the FULL product
// T1 x T2.  Targets among remaining gaps: [2,40] ([8]x[10]), [2,8,8]
// ([8]x[2,8]), [10,10] ([10]x[10]), [2,2,20] ([10]x[2,4]), [2,2,40]
// ([10]x[2,8]), [56] ([7]x[8]), [2,28] ([7]x[2,4]), [2,56] ([7]x[2,8]).
// Surface: z^2 = f(u,v) = 108^2 u^2 + c1(v) u + c0(v); writing
// w = 108u + c1/216 gives (z-w)(z+w) = B(v), so every (v,d) with d | B(v)
// formally (d rational, z-w=d, z+w=B/d) is a rational point: fully rational
// sweep.  Twist alignment for E2 via mod-3 trace matching.
// Usage: magma -b lane_z3.m > ../logs/z3.log
SetColumns(0);
SetMemoryLimit(4*10^9);
load "split_lab.m";  // run from product/code/

Qv<v> := FunctionField(Rationals());
c1 := -432*v^3 - 216*v^2 + 576*v - 184;
c0 := 4*v^6 - 20*v^5 + 41*v^4 - 44*v^3 + 26*v^2 - 8*v + 1;
Bv := c0 - c1^2/46656;

// WNr moduli for (3,2), from Frengley ZNr-equations.m (jj, jj_1728 as
// functions of (u,v)); we avoid loading the whole package by recomputing the
// j-pair from his WNrModuli intrinsic — so attach the spec instead.
// External dependency (beyond the Magma spec noted in split_lab.m): Sam
// Frengley's N-congruences package, expected at product/N-congruences (it is
// .gitignored here — a separate project with its own license; obtain it via
//   git clone https://github.com/SamFrengley/N-congruences product/N-congruences
// ), or point the environment variable ZNR_EQUATIONS at ZNr-equations.m.
znr := GetEnv("ZNR_EQUATIONS");
if znr eq "" then znr := "../N-congruences/ZNr-equations.m"; end if;
attached := true;
try
    Attach(znr);
catch e
    attached := false;
end try;
if not attached then
    printf "lane_z3: cannot attach %o -- clone https://github.com/SamFrengley/N-congruences into product/ or set ZNR_EQUATIONS\n", znr;
    quit;
end if;

function PotentialOrder(E0)
    // largest torsion order achievable among quadratic twists of E0, detected
    // via rational roots of division polynomials (x-coords are twist-covariant)
    best := 1;
    for m in [3,5,7,9] do
        if #Roots(DivisionPolynomial(E0, m)) gt 0 then best := Max(best, m); end if;
    end for;
    f2 := DivisionPolynomial(E0, 2);
    r2 := #Roots(f2);
    p2 := 1;
    if r2 gt 0 then
        p2 := 2;
        if #Roots(DivisionPolynomial(E0, 4) div f2) gt 0 then
            p2 := 4;
            f8 := DivisionPolynomial(E0, 8) div DivisionPolynomial(E0, 4);
            if #Roots(f8) gt 0 then p2 := 8; end if;
        end if;
    end if;
    return best*p2, r2;
end function;

function AlignTwist(E1, E2)
    // find small quadratic twist E2^d that is 3-congruent to E1 (a_p match mod 3)
    bad := &*PrimeDivisors(Integers()!(Discriminant(E1)*Discriminant(E2))*6);
    ds := Divisors(bad);
    ds := Sort(ds cat [-d : d in ds], func<x,y|AbsoluteValue(x)-AbsoluteValue(y)>);
    ps := [p : p in PrimesInInterval(5,200) | bad mod p ne 0];
    for d in ds do
        if AbsoluteValue(d) gt 1000 then break; end if;
        Ed := QuadraticTwist(E2, d);
        ok := true;
        for p in ps do
            if (TraceOfFrobenius(E1,p) - TraceOfFrobenius(Ed,p)) mod 3 ne 0 then ok := false; break; end if;
        end for;
        if ok then return true, Ed; end if;
    end for;
    return false, E2;
end function;

function HeightRats(H)
    S := [];
    for a in [1..H], b in [1..H] do
        if GCD(a,b) eq 1 then Append(~S, a/b); Append(~S, -a/b); end if;
    end for;
    return S;
end function;

seen := {};
jseen := {};
nstat := AssociativeArray();
for st in ["skip","abort","known","hit","fail"] do nstat[st] := 0; end for;
npt := 0; ncand := 0; nglue := 0;
t0 := Cputime();
vs := HeightRats(8);
dsw := HeightRats(8);
for vv in vs do
    Bval := Evaluate(Bv, vv);
    if Bval eq 0 then continue; end if;
    c1v := Evaluate(c1, vv);
    for d in dsw do
        w := (Bval/d - d)/2;
        zz := (Bval/d + d)/2;
        uu := (w - c1v/216)/108;
        if uu eq 0 then continue; end if;
        npt +:= 1;
        // moduli
        E1 := 0; E2 := 0; okm := true;
        try
            E1, E2 := ZNrModuli(3, 2, [Rationals()|uu, vv, zz] : bound := 200);
        catch e
            okm := false;
        end try;
        if not okm then continue; end if;
        jp := {jInvariant(E1), jInvariant(E2)};
        if jp in jseen then continue; end if;
        Include(~jseen, jp);
        // cheap potential filter on both components
        o1 := PotentialOrder(E1);
        o2 := PotentialOrder(E2);
        if o1*o2 lt 40 then continue; end if;
        ncand +:= 1;
        // find the twists with actual torsion: scan small twists for the max torsion
        T1best := E1; T2best := E2; n1 := #TorsionSubgroup(E1); n2 := #TorsionSubgroup(E2);
        for d1 in [-30..30] do
            if d1 eq 0 or not IsSquarefree(AbsoluteValue(d1)) then continue; end if;
            Et := QuadraticTwist(E1, d1);
            if #TorsionSubgroup(Et) gt n1 then T1best := Et; n1 := #TorsionSubgroup(Et); end if;
            Et2 := QuadraticTwist(E2, d1);
            if #TorsionSubgroup(Et2) gt n2 then T2best := Et2; n2 := #TorsionSubgroup(Et2); end if;
        end for;
        if n1*n2 lt 40 then continue; end if;
        // align E2-twist to be 3-congruent with the chosen E1-twist
        okal, E2al := AlignTwist(T1best, T2best);
        if not okal then continue; end if;
        printf "PAIR v=%o d=%o : tors %o x %o (conds %o, %o)\n", vv, d,
            Invariants(TorsionSubgroup(T1best)), Invariants(TorsionSubgroup(E2al)),
            Conductor(T1best), Conductor(E2al);
        nglue +:= 1;
        L := [];
        try L := Genus2Elliptic3(T1best, E2al); catch e L := []; end try;
        for k in [1..#L] do
            gk := "";
            try gk := Sprintf("%o", G2Invariants(L[k])); catch e gk := "bad"; end try;
            if gk eq "bad" or gk in seen then continue; end if;
            Include(~seen, gk);
            st := Funnel(L[k], Sprintf("z3|v=%o,d=%o|%o", vv, d, k));
            nstat[st] +:= 1;
        end for;
        if nglue mod 10 eq 0 then
            printf "PROGRESS z3 pts %o cands %o glues %o hits %o %o s\n", npt, ncand, nglue, nstat["hit"], Cputime()-t0;
        end if;
    end for;
end for;
printf "SEARCH_DONE z3 pts %o cands %o glues %o aborts %o known %o hits %o fails %o %o s\n",
    npt, ncand, nglue, nstat["abort"], nstat["known"], nstat["hit"], nstat["fail"], Cputime()-t0;
quit;
