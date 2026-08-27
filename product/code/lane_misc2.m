// lane_misc2.m — cleanup analyses for the 2026-08-13 session:
// (1) 2-descent rank bounds on the genus-3 TWISTED diagonal condition curves
//     for N=12 and N=10 (rank bound 0 would PROVE the twisted self-gluing
//     diagonals empty, completing the impossibility of [12,12]/[10,10] via
//     conjugation self-gluings);
// (2) small-point census of the five non-identity sigma-surfaces (2-of-3
//     conditions) for the [2,2,6,6] cross search, to decide which fibrations
//     deserve a structured attack.
// Usage: cd product/code && magma -b lane_misc2.m > ../logs/lane_misc2.log
SetColumns(0);
SetMemoryLimit(5*10^9);
SetClassGroupBounds("GRH");

RQ := Rationals();
QY<T> := PolynomialRing(RQ);

// ---------- (1) twisted diagonal curves ----------
printf "== twisted diagonal curves ==\n";
tw12 := -1296*T^8 + 5184*T^7 - 9072*T^6 + 9072*T^5 - 5580*T^4 + 2088*T^3 - 432*T^2 + 36*T;
tw10 := 32*T^7 - 160*T^6 + 256*T^5 - 156*T^4 + 16*T^3 + 16*T^2 - 4*T;
for entry in [* <"N=12 twisted", tw12>, <"N=10 twisted", tw10> *] do
    q := entry[2];
    // strip square content to get the squarefree kernel model
    fac := Factorization(q);
    unit := q div &*[ f[1]^f[2] : f in fac ];
    cls := QY!(Sign(RQ!Coefficient(unit,0))*SquarefreeFactorization(Integers()!Abs(RQ!Coefficient(unit,0))));
    for f in fac do if IsOdd(f[2]) then cls *:= f[1]; end if; end for;
    den := LCM([ Denominator(c) : c in Coefficients(cls) ]);
    cls := QY!(cls*den^2);
    printf "%o: kernel w^2 = %o\n", entry[1], cls;
    C := HyperellipticCurve(cls);
    printf "   genus %o\n", Genus(C);
    J := Jacobian(C);
    rb := RankBound(J);
    printf "   RankBound(J) = %o\n", rb;
    pts := Points(C : Bound := 20000);
    printf "   points to 20000: %o\n", pts;
    if rb eq 0 then
        printf "   RANK 0: torsion-only; the listed points are ALL rational points -> twisted diagonal PROVABLY empty if all degenerate\n";
    end if;
end for;

// ---------- (2) sigma-surface point census for [2,2,6,6] ----------
// conditions (from lane_2266_derive.m):
//   pair classes a12=(t+3)(t-5), a13=2(t-3), a23=-(t-1)(t-9), transposes negate.
printf "== sigma-surface census (2-of-3 conditions, H=40) ==\n";
function SFrat(x)
    n := Numerator(x)*Denominator(x);
    s := Sign(n); n := Abs(n);
    return s*SquarefreeFactorization(n);
end function;
function HeightRats(H)
    S := [];
    for a in [1..H], b in [1..H] do
        if GCD(a,b) eq 1 then Append(~S, a/b); Append(~S, -a/b); end if;
    end for;
    return S;
end function;
EXCL := {RQ|3,-3,1,5,9};
vals := [ v : v in HeightRats(40) | not v in EXCL ];
NV := #vals;
kk := [];
for v in vals do
    Append(~kk, [ SFrat((v+3)*(v-5)), SFrat(2*(v-3)), SFrat(-(v-1)*(v-9)) ]);
end for;
function KO(k, a, b)
    if a lt b then
        return a eq 1 select (b eq 2 select k[1] else k[2]) else k[3];
    else
        return -(b eq 1 select (a eq 2 select k[1] else k[2]) else k[3]);
    end if;
end function;
SIGMAS := [ [1,2,3],[2,1,3],[3,2,1],[1,3,2],[2,3,1],[3,1,2] ];
PAIRS := [ [1,2],[1,3],[2,3] ];
// for each sigma and each CHOICE of 2 of the 3 pair conditions, count matches
for si in [1..6] do
    sg := SIGMAS[si];
    for drop in [1..3] do
        keep := [ j : j in [1..3] | j ne drop ];
        M := AssociativeArray();
        for ui in [1..NV] do
            key := < KO(kk[ui],sg[PAIRS[j][1]],sg[PAIRS[j][2]]) : j in keep >;
            if IsDefined(M,key) then M[key] +:= 1; else M[key] := 1; end if;
        end for;
        nmatch := 0;
        for ti in [1..NV] do
            key := < kk[ti][j] : j in keep >;
            if IsDefined(M,key) then nmatch +:= M[key]; end if;
        end for;
        printf "sigma=%o drop=%o: %o (t,u) surface points at H=40\n", sg, drop, nmatch;
    end for;
end for;
printf "MISC2_DONE\n";
quit;
