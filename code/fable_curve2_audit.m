//////////////////////////////////////////////////////////////////////
// fable_curve2_audit.m  (2026-07-19, Fable session — last day!)
//
// Full Richelot-component + divisibility audit of (2,2,2,12) CURVE #2
// ((s,m,n)=(2208,-8303,-7200) on M(2,2,2,6); PR #3/#4, merged today).
// The FIRST record's audit found: exactly one of the 15 nonzero J[2]
// classes divisible by 2 (the existing 4-direction), order-12 generator
// NOT divisible, component = 18 vertices, no order>96 vertex.
// Nobody has audited curve #2.  Any second divisible 2-class here
// means [2,2,4,12] (order 192); a divisible 12-generator means
// [2,2,2,24] (order 192); any component vertex with order>=96 beyond
// the seed is a new realization.
//
// Run:  magma -b code/fable_curve2_audit.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(2*10^9);

Q := Rationals();
P<x> := PolynomialRing(Q);

function IntegralSquareScale(f)
    den := LCM([Denominator(Coefficient(f,i)) : i in [0..Degree(f)]]);
    return P!(den^2*f), den;
end function;

function NormalizePoly(D)
    fD, hD := HyperellipticPolynomials(D);
    FD := hD eq 0 select P!fD else P!(hD^2 + 4*fD);
    FI, _ := IntegralSquareScale(FD);
    return FI;
end function;

function G2Key(f)
    return Sprintf("%o", [Q!a : a in G2Invariants(HyperellipticCurve(f))]);
end function;

// curve #2
A := [1,1,1,2,2];
B := [25648128,-36568896,-52466496,-59781600,23309856];
f2 := &*[A[i] + B[i]*x : i in [1..5]];
assert Degree(f2) eq 5 and Discriminant(f2) ne 0;

C2 := HyperellipticCurve(f2);
J2 := Jacobian(C2);
time T2, mp2 := TorsionSubgroup(J2);
inv2 := Invariants(T2);
printf "SEED torsion=%o order=%o\n", inv2, #T2;
assert inv2 eq [2,2,2,12];

// ---- divisibility audit on the seed ----
// (a) order-12 generators: divisible by 2?  (b) all 15 nonzero 2-torsion
// classes; (c) all order-4 classes.
n12 := 0; n12div := 0;
n2 := 0; n2div := 0;
n4 := 0; n4div := 0;
for g in T2 do
    o := Order(g);
    D := mp2(g);
    if o eq 12 then
        n12 +:= 1;
        ok := false; try ok, h := IsDivisibleBy(D, 2); catch e ok := false; end try;
        if ok then n12div +:= 1; printf "DIV12! a 12-generator halves: [2,2,2,24] chain present\n"; end if;
    elif o eq 2 then
        n2 +:= 1;
        ok := false; try ok, h := IsDivisibleBy(D, 2); catch e ok := false; end try;
        if ok then n2div +:= 1; printf "DIV2 class %o of 15 halves\n", n2; end if;
    elif o eq 4 then
        n4 +:= 1;
        ok := false; try ok, h := IsDivisibleBy(D, 2); catch e ok := false; end try;
        if ok then n4div +:= 1; printf "DIV4! an order-4 class halves: 8-chain present\n"; end if;
    end if;
end for;
printf "DIVAUDIT ord12=%o div=%o | ord2=%o div=%o | ord4=%o div=%o\n",
    n12, n12div, n2, n2div, n4, n4div;

// ---- Richelot BFS (depth 3) with exact torsion at every vertex ----
max_depth := 3;
visited := {};
queue := [];
fI, _ := IntegralSquareScale(f2);
Include(~visited, G2Key(fI));
Append(~queue, <fI, 0, "S">);
cursor := 1;
tally := AssociativeArray();
while cursor le #queue do
    node := queue[cursor];
    f := node[1]; depth := node[2]; tag := node[3];
    J := Jacobian(HyperellipticCurve(f));
    T := TorsionSubgroup(J);
    inv := Invariants(T);
    ord := #inv eq 0 select 1 else &*inv;
    key := Sprintf("%o", inv);
    if IsDefined(tally, key) then tally[key] +:= 1; else tally[key] := 1; end if;
    printf "VERTEX %o depth=%o inv=%o\n", tag, depth, inv;
    if ord ge 96 and depth gt 0 then
        printf "BIG %o inv=%o order=%o f=%o\n", tag, inv, ord, f;
    end if;
    if inv in {[2,2,4,12],[2,2,2,24],[4,12],[2,4,12],[2,2,24],[2,24]} then
        printf "TARGETHIT %o inv=%o f=%o\n", tag, inv, f;
    end if;
    if depth lt max_depth then
        Rs := [];
        try Rs := RichelotIsogenousSurfaces(J);
        catch e printf "RICHELOT_FAIL %o\n", tag; end try;
        for i in [1..#Rs] do
            S := Rs[i]; tS := Type(S);
            if tS eq JacHyp or tS eq CrvHyp then
                D := tS eq JacHyp select Curve(S) else S;
                fN := NormalizePoly(D);
                k := G2Key(fN);
                if k notin visited then
                    Include(~visited, k);
                    Append(~queue, <fN, depth+1, Sprintf("%o.%o", tag, i)>);
                end if;
            else
                printf "NONJAC %o.%o type=%o\n", tag, i, tS;
            end if;
        end for;
    end if;
    cursor +:= 1;
end while;
printf "COMPONENT_DONE vertices=%o depth=%o inventory:\n", #queue, max_depth;
for key in Keys(tally) do printf "  %o x %o\n", tally[key], key; end for;
quit;
