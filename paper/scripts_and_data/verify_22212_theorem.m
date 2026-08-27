// verify_22212_theorem.m -- Theorem 3.2 of the paper
// ("three points on the surface").
//
// For each of the three points P = [a:b:c:u:v] listed in the theorem this
// script checks:
//   (1) P lies on S: a^2+b^2+c^2 = u^2+v^2 and a^4+b^4+c^4 = u^4+v^4,
//       and P is in the open set S^o (coordinates nonzero, squares distinct);
//   (2) X_P : y^2 = x(x-a^2)(x-b^2)(x-c^2)(x-u^2)(x-v^2) has
//       Jac(X_P)(Q)_tors of invariant type exactly [2,2,2,12] (order 96);
//   (3) Jac(X_P) is geometrically simple, via the good-prime Frobenius
//       root-power certificate: an admissible prime p with chi_p irreducible
//       and Degree(MinimalPolynomial(pi^n)) = 4 for all n in [2..12]
//       (the witness prime and chi_p are printed);
//   (4) the three curves are pairwise non-isomorphic over Qbar
//       (distinct absolute Igusa invariants).
// Run from this directory:  magma -b verify_22212_theorem.m
SetColumns(0);
SetMemoryLimit(8*10^9);
Q := Rationals(); Z := Integers();
P<x> := PolynomialRing(Q);
T<t> := PolynomialRing(Q);

pts := [
  [120, 143, 266, 218, 241],
  [143, 408, 1015, 437, 1013],
  [16660, 78793, 21456, 78644, 27593]
];

// good-prime Frobenius root-power certificate for geometric simplicity
// (simplicity-certificates convention: irreducible chi_p with no degree
// drop of pi^n for n in [2..12] at one good prime certifies J geometrically
// simple).
function SimplicityCertificate(f)
    dsc := Z!Discriminant(f);
    lc := Z!LeadingCoefficient(f);
    for pp in PrimesInInterval(3, 500) do
        if dsc mod pp eq 0 or lc mod pp eq 0 then continue; end if;
        fp := PolynomialRing(GF(pp))!f;
        if Degree(fp) lt 5 or not IsSquarefree(fp) then continue; end if;
        chi := T!Reverse(Coefficients(
            EulerFactor(Jacobian(HyperellipticCurve(fp)))));
        if Degree(chi) ne 4 or not IsIrreducible(chi) then continue; end if;
        K := NumberField(chi); pi := K.1;
        drop := false;
        for nn in [2..12] do
            if Degree(MinimalPolynomial(pi^nn)) lt 4 then
                drop := true; break;
            end if;
        end for;
        if not drop then return true, pp, chi; end if;
    end for;
    return false, 0, T!0;
end function;

invs := [];
t0 := Cputime();
for i in [1..#pts] do
    a, b, c, u, v := Explode(pts[i]);
    // (1) surface membership and the open conditions
    assert a^2 + b^2 + c^2 eq u^2 + v^2;
    assert a^4 + b^4 + c^4 eq u^4 + v^4;
    sq := [a^2, b^2, c^2, u^2, v^2];
    assert 0 notin sq and #Seqset(sq) eq 5;
    printf "point %o = %o lies on S^o\n", i, pts[i];
    // (2) exact torsion
    f := x * &*[x - s : s in sq];
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    I := Invariants(AbelianGroup(TorsionSubgroup(J)));
    printf "point %o: torsion invariants %o\n", i, I;
    assert I eq [2, 2, 2, 12];
    // (3) geometric simplicity
    ok, pp, chi := SimplicityCertificate(f);
    assert ok;
    printf "point %o: geometrically simple, witness p=%o, chi=%o\n", i, pp, chi;
    Append(~invs, G2Invariants(C));
end for;
// (4) pairwise non-isomorphic (distinct absolute Igusa invariants)
assert #Seqset(invs) eq 3;
printf "the three curves are pairwise non-isomorphic\n";
printf "THEOREM 3.2 VERIFIED: 3 curves, exact [2,2,2,12], simple, distinct (%.1o s)\n", Cputime() - t0;
quit;
