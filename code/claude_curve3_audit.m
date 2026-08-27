// claude_curve3_audit.m — standing-rule depth-3 Richelot + divisibility audit of
// (2,2,2,12) CURVE #3 (surface point (5364:19661:4165), found 2026-07-20).
// Tests the T5-component-rigidity conjecture (fable_final_top10_2026_07_19):
// predicted census 1x[2,2,2,12] + 12x[2,12] + 5x[2,2,6], exactly one divisible
// 2-class, no divisible order-4 or order-12 classes.
// Adapted from code/fable_curve2_audit.m.
SetColumns(0);
SetMemoryLimit(6*10^9);
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

// curve #3 minimal model: y^2 + (x^2+x)y = fmin
fmin := 3703062294195264*x^6 - 360079374491052216*x^5 + 8901721379573296848*x^4
  - 5397945250386334945*x^3 - 86737535708373850908*x^2 + 36346694984390901540*x
  + 43035470132681030400;
f3 := (x^2+x)^2 + 4*fmin;   // simplified model y^2 = f3
assert Discriminant(f3) ne 0;

C3 := HyperellipticCurve(f3);
J3 := Jacobian(C3);
time T3, mp3 := TorsionSubgroup(J3);
inv3 := Invariants(T3);
printf "SEED torsion=%o order=%o\n", inv3, #T3;
assert inv3 eq [2,2,2,12];

n12 := 0; n12div := 0;
n2 := 0; n2div := 0;
n4 := 0; n4div := 0;
for g in T3 do
    o := Order(g);
    D := mp3(g);
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

max_depth := 3;
visited := {};
queue := [];
fI, _ := IntegralSquareScale(f3);
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
