// claude_ov_b4_orbit6members.m -- Lane 4 (route B4): the SPORADIC Richelot
// sources of Flynn's order-11 family.
//
// claude_ov_b4_orbitcovers.m proved: the locus of Flynn members carrying a
// Galois-stable (2,2)-kernel is a cover of the t-line with exactly two
// irreducible components,
//   * the "D4" component (orbit of size 9, stabiliser D4) -- GENUS 1, and it
//     is the (6,2) component E = 92.a1, rank 1: infinitely many members, all
//     already sieved for 2-rank raises to |n| <= 10^7 with ZERO survivors;
//   * the "D6" component (orbit of size 6, stabiliser D6, the TRANSVERSE
//     partitions, whose three quadratics are permuted by Galois and need not
//     be individually rational) -- GENUS 2, hence FINITELY many members.
// Its plane model has rational points at t = 0 (degenerate), -1, 65/4, -74/9.
//
// This script takes every t found on the genus-2 component and runs the exact
// test: all Richelot isogenies over Q, the 2-rank and the exact torsion of
// every codomain.  A codomain with 2-rank >= 2 would be a [2,22] curve.
//
// Run: code/claude_magma_slot.sh -b MemGB:=12 code/claude_ov_b4_orbit6members.m \
//        > results/claude_ov_b4_orbit6members.log 2>&1 &

SetColumns(0);
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();
Pol<x> := PolynomialRing(Q);

// y^2 = f  ->  an EQUIVALENT integral model, via the content:
// f = (n/d)*g with g integral primitive, so (d*y)^2 = n*d*g.
// (TorsionSubgroup refuses non-integral models; scaling x -> X/m only works
//  when f is monic, which Richelot codomains are not.)
function IntegralSextic(f)
    cs := [c : c in Coefficients(f) | c ne 0];
    d := LCM([Denominator(c) : c in cs]);
    n := GCD([Numerator(d*c) : c in cs]);
    g := (d/n)*f;
    h := (n*d)*g;
    assert &and[IsIntegral(c) : c in Coefficients(h)];
    return Parent(f)![Integers()!c : c in Coefficients(h)];
end function;

function TwoRankFromFactors(f)
    degs := [Degree(q[1]) : q in Factorization(f)];
    if &+degs eq 5 then degs := degs cat [1]; end if;
    k := #degs;
    nev := 0;
    for m in [0..2^k-1] do
        s := 0;
        for i in [1..k] do if (m div 2^(i-1)) mod 2 eq 1 then s +:= degs[i]; end if; end for;
        if s mod 2 eq 0 then nev +:= 1; end if;
    end for;
    return Ilog2(nev) - 1, Sort(degs);
end function;

TS := [Rationals() | -1, 65/4, -74/9, 0 ];
for tv in TS do
    f := x^6 + 2*x^5 + (2*tv+3)*x^4 + 2*x^3 + (tv^2+1)*x^2 + 2*tv*(1-tv)*x + tv^2;
    printf "\n=== t = %o ===\n", tv;
    if Discriminant(f) eq 0 then printf "  SINGULAR -- skipped\n"; continue; end if;
    g := IntegralSextic(f);
    C := HyperellipticCurve(g);
    J := Jacobian(C);
    tr, degs := TwoRankFromFactors(g);
    printf "  seed factor type %o, 2-rank %o\n", degs, tr;
    printf "  seed TORSION = %o\n", Invariants(TorsionSubgroup(J));
    printf "  seed G2Invariants = %o\n", G2Invariants(C);
    RS := RichelotIsogenousSurfaces(J);
    printf "  #RichelotIsogenousSurfaces = %o types %o\n", #RS, [Type(S) : S in RS];
    for i in [1..#RS] do
        S := RS[i];
        if Type(S) ne JacHyp then
            printf "  [%o] NOT a Jacobian (%o) -- the codomain is a product of elliptic curves\n", i, Type(S);
            continue;
        end if;
        D := Curve(S);
        fD := HyperellipticPolynomials(SimplifiedModel(D));
        fD := IntegralSextic(fD);
        DD := HyperellipticCurve(fD);
        trD, degsD := TwoRankFromFactors(fD);
        printf "  [%o] codomain f = %o\n", i, fD;
        printf "      factor type %o  2-RANK %o  TORSION %o\n",
            degsD, trD, Invariants(TorsionSubgroup(Jacobian(DD)));
    end for;
end for;

printf "ORBIT6MEMBERS_DONE\n";
quit;
