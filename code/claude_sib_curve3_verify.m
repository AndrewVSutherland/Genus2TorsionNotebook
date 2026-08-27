// Independent re-verification of (2,2,2,12) curve #3 (ABC surface orbit 3).
// Fresh-session rebuild via the SYMMETRIC MODEL route -- a different code
// path from code/verify_22212_abc_hit.m (which goes (A,B,C)->(s,m,n)->
// B-vector quintic).  Here: (A,B,C) -> symmetric point [x:y:z:u:v] on
//   S: x^2+y^2+z^2 = u^2+v^2,  x^4+y^4+z^4 = u^4+v^4
// -> sextic C_P: W^2 = T(T-x^2)(T-y^2)(T-z^2)(T-u^2)(T-v^2)
// (notes/22212_symmetric_model.md), then exact torsion, strict simplicity
// certificates at FRESH primes disjoint from {37,127,131,179}, and identity
// with the recorded reduced minimal model (data/22212_abc_curve3.txt).
//
// Run: magma -b code/claude_sib_curve3_verify.m > results/claude_sib_curve3_verify.log

SetColumns(0);
SetSeed(1);
SetMemoryLimit(3*10^9);

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

// Direct-search chart of orbit 3 (smallest chart found by search_22212_abc.cpp)
A := 4165; B := 19661; Cc := 5364;

X := A^2; Y := B^2; Zc := Cc^2;
D := X - Y + Zc;
L := (X - Y)*(Y - Zc);
F := D*L;
G := D*(L + D^2);
assert D ne 0 and F gt 0 and G gt 0 and IsSquare(F) and IsSquare(G);
y0 := Isqrt(F); z0 := Isqrt(G);
assert y0^2 eq F and z0^2 eq G and G - F eq D^3;
printf "ABC = (%o,%o,%o)  D = %o  y0 = %o  z0 = %o\n", A, B, Cc, D, y0, z0;

// Birational dictionary (notes/22212_symmetric_model.md):
// [x:y:z:u:v] = [A*D : C*D : y0 : B*D : z0], reduced to a primitive vector.
raw := [A*D, Cc*D, y0, B*D, z0];
g := GCD([Abs(t) : t in raw]);
sym := [Abs(t) div g : t in raw];
printf "SYMMETRIC_POINT = %o\n", sym;

sx := sym[1]; sy := sym[2]; sz := sym[3]; su := sym[4]; sv := sym[5];
assert sx^2 + sy^2 + sz^2 eq su^2 + sv^2;
assert sx^4 + sy^4 + sz^4 eq su^4 + sv^4;
printf "SYMMETRIC_EQUATIONS_HOLD = true\n";

// Nondegeneracy: all coordinates nonzero, five squares pairwise distinct.
sq := [sx^2, sy^2, sz^2, su^2, sv^2];
assert &and[t ne 0 : t in sq];
assert #Seqset(sq) eq 5;
printf "NONDEGENERATE = true\n";

// The symmetric-model curve.
f6 := x * &*[x - t : t in sq];
assert Degree(f6) eq 6 and Discriminant(f6) ne 0;
C0 := HyperellipticCurve(f6);

Cmin := ReducedMinimalWeierstrassModel(C0);
fmin, hmin := HyperellipticPolynomials(Cmin);
printf "REDUCED_MINIMAL_POLYNOMIALS = %o, %o\n", fmin, hmin;

// Recorded reduced minimal model (data/22212_abc_curve3.txt).
frec := 3703062294195264*x^6 - 360079374491052216*x^5
      + 8901721379573296848*x^4 - 5397945250386334945*x^3
      - 86737535708373850908*x^2 + 36346694984390901540*x
      + 43035470132681030400;
hrec := x^2 + x;
Crec := HyperellipticCurve(frec, hrec);

match := (fmin eq frec) and (hmin eq hrec);
printf "MINIMAL_MODEL_MATCHES_RECORDED = %o\n", match;
if not match then
    iso := IsIsomorphic(Cmin, Crec);
    printf "ISOMORPHIC_TO_RECORDED = %o\n", iso;
    assert iso;
end if;

// Exact rational torsion (integral y^2 = f model of the minimal curve).
Cs := SimplifiedModel(Cmin);
J := Jacobian(Cs);
T := TorsionSubgroup(J);
invs := Invariants(T);
printf "TORSION_INVARIANTS = %o\n", invs;
printf "TORSION_ORDER = %o\n", #T;
assert invs eq [2,2,2,12] and #T eq 96;

// Strict geometric-simplicity witnesses at FRESH primes: chi irreducible of
// degree 4 AND deg MinPoly(pi^e) = 4 for all 2 <= e <= 12.  The four primes
// already used by verify_22212_abc_hit.m are excluded.
used := {37, 127, 131, 179};
disc := Z!Discriminant(Cmin);
nc := 0;
for pp in PrimesInInterval(41, 500) do
    if pp in used or disc mod pp eq 0 then continue; end if;
    Cp := ChangeRing(Cmin, GF(pp));
    chi := P!Reverse(Coefficients(LPolynomial(Cp)));
    if not IsIrreducible(chi) or Degree(chi) ne 4 then continue; end if;
    K<pi> := NumberField(chi);
    degs := [Degree(MinimalPolynomial(pi^e)) : e in [2..12]];
    if &and[d eq 4 : d in degs] then
        nc +:= 1;
        printf "FRESH_SIMPLICITY_WITNESS p=%o chi=%o\n", pp, chi;
        if nc ge 3 then break; end if;
    end if;
end for;
printf "FRESH_WITNESS_COUNT = %o\n", nc;
assert nc ge 2;

// Q-bar isomorphism class: distinct from curves #1 and #2, same as the
// quintic discovery model rebuilt from (s,m,n).
function CurveFromSMN(ss, mm, nn)
    avec := [1,1,1,2,2];
    bvec := [
        2*ss^2 - ss*nn,
        2*ss^2 + ss*mm - 2*ss*nn - mm*nn,
        2*ss^2 + ss*mm - ss*nn - mm*nn,
        -mm*nn,
        4*ss^2 - 4*ss*nn - mm*nn
    ];
    ff := &*[avec[i] + bvec[i]*x : i in [1..5]];
    assert Degree(ff) eq 5 and Discriminant(ff) ne 0;
    return HyperellipticCurve(ff);
end function;

g2 := G2Invariants(Cmin);
printf "G2_INVARIANTS = %o\n", g2;
assert g2 eq G2Invariants(CurveFromSMN(254097487, -10481401502, 555111200));
printf "MATCHES_QUINTIC_DISCOVERY_MODEL_G2 = true\n";
assert g2 ne G2Invariants(CurveFromSMN(336396, -689185, -166464));
assert g2 ne G2Invariants(CurveFromSMN(2208, -8303, -7200));
printf "DISTINCT_FROM_CURVES_1_AND_2 = true\n";

printf "REVERIFIED_22212_CURVE3_VIA_SYMMETRIC_MODEL\n";
quit;
