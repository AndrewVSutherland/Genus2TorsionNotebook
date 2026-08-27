// claude_ov_b4_2214_richelot.m -- Lane 4, fallback item 6:
// the project's standing rule "fan every new simple curve of torsion order >= 48
// through a depth-1..3 Richelot audit" applied, for the first time, to the
// eleven contact-7 three-root [2,2,14] curves (six from 2026-07-23, five new
// from the 2026-07-25 harvest).
//
// The quintic has factor type [1,1,1,2]; with the Weierstrass point at infinity
// the six Weierstrass points form Galois orbits {r1},{r2},{r3},{conj pair},{inf}.
// A Galois-stable partition into three pairs must pair inf with one of the three
// rational roots, then the other two rational roots together and the conjugate
// pair together: exactly THREE Galois-stable (2,2)-kernels, all defined over Q.
// A (2,2)-isogeny is injective on odd torsion (so the 7 survives) and preserves
// End^0 (so genericity survives).  Any neighbour with torsion [2,2,28] or
// [2,2,2,14] has ORDER 112 -- a new record for a geometrically simple genus-2
// Jacobian over Q.
//
//   magma -b IDX:=1 DEPTH:=3 MemGB:=3 claude_ov_b4_2214_richelot.m
//
// Markers: SEEDINFO, NODE, NEWGROUP, RECORD, SEEDDONE.

SetColumns(0);
if not assigned IDX   then IDX := 1;   elif Type(IDX) eq MonStgElt then IDX := StringToInteger(IDX); end if;
if not assigned DEPTH then DEPTH := 3; elif Type(DEPTH) eq MonStgElt then DEPTH := StringToInteger(DEPTH); end if;
if not assigned MemGB then MemGB := 3; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals(); P<x> := PolynomialRing(Q); Z := Integers();
G := func<w | -(w^5 - w^3 - w^2/2)/(w+1)^2>;

// six original orbits (2026-07-23; #1 is the RM(sqrt2) one) then the five new
trips := [
  [-3, -3/4, -3/5],            // 1  RM(sqrt2) signature
  [-10, -1/2, -10/7],          // 2
  [-5, -15/8, -15/22],         // 3
  [-1/2, -15/8, -15/19],       // 4
  [-4/9, 4/17, -4/25],         // 5
  [4/17, -5/18, -10/49],       // 6
  [-511/61, -511/625, -1/2],   // 7  NEW 2026-07-25
  [-165/41, -33/16, -165/289], // 8  NEW
  [-164/297, -1/2, 164/361],   // 9  NEW
  [-17/50, -34/189, 34/121],   // 10 NEW
  [-1/2, -13/49, 13/50]        // 11 NEW
];

function Norm2(C)
    try C := ReducedMinimalWeierstrassModel(C); catch e ; end try;
    C := SimplifiedModel(C);
    try C := ReducedModel(C); C := SimplifiedModel(C); catch e ; end try;
    return C;
end function;

T := trips[IDX];
v1 := T[1]; v2 := T[2];
c4 := (G(v1)-G(v2))/(v1^2-v2^2); c0 := G(v1) - c4*v1^2;
b := c4 - 2; a := 9/2 - c0 - c4;
h := 1 - (7/2)*x + a*x^2 + b*x^3;
num := h^2 + (x-1)^7;  f := num div x^2;
assert num eq f*x^2 and Degree(f) eq 5 and Discriminant(f) ne 0;
m := LCM([Denominator(c) : c in Coefficients(f)]);
F := P![ Z!(Coefficient(f,i)*m^(6-i)) : i in [0..5] ];
C0 := HyperellipticCurve(F);
J0 := Jacobian(C0);
tor0 := Invariants(TorsionSubgroup(J0));
printf "SEEDINFO idx=%o trip=%o ftype=%o torsion=%o\n", IDX, T,
       Sort([Degree(t[1]) : t in Factorization(F)]), tor0;

// NB: dedupe by exact Q-isomorphism, NOT by G2Invariants -- G2Invariants are
// geometric and cannot separate quadratic twists, which have DIFFERENT torsion.
seen := [* *];
queue := [* <C0, 0> *];
nnodes := 0; ngroups := {}; nrecord := 0;

while #queue gt 0 do
    it := queue[1]; Remove(~queue, 1);
    C := it[1]; d := it[2];
    dup := false;
    for CC in seen do
        if IsIsomorphic(CC, C) then dup := true; break; end if;
    end for;
    if dup then continue; end if;
    Append(~seen, C);
    nnodes +:= 1;
    J := Jacobian(C);
    ff := HyperellipticPolynomials(C);
    r := #Invariants(TwoTorsionSubgroup(J));
    tor := Invariants(TorsionSubgroup(J));
    ordt := #TorsionSubgroup(J);
    Include(~ngroups, tor);
    printf "NODE idx=%o depth=%o 2rank=%o torsion=%o order=%o ftype=%o\n",
           IDX, d, r, tor, ordt, Sort([Degree(t[1]) : t in Factorization(ff)]);
    if ordt gt 56 or (ordt eq 56 and tor ne tor0) then
        nrecord +:= 1;
        printf "NEWGROUP idx=%o depth=%o torsion=%o order=%o f=%o\n", IDX, d, tor, ordt, ff;
        if ordt ge 112 then
            printf "RECORD idx=%o depth=%o torsion=%o order=%o f=%o\n", IDX, d, tor, ordt, ff;
        end if;
    end if;
    if d ge DEPTH then continue; end if;
    RS := [];
    try
        RS := RichelotIsogenousSurfaces(J);
    catch e
        printf "RICHFAIL idx=%o depth=%o\n", IDX, d; continue;
    end try;
    printf "RSCOUNT idx=%o depth=%o nRS=%o nJac=%o\n", IDX, d, #RS,
           #[1 : S in RS | Type(S) eq JacHyp];
    for S in RS do
        if Type(S) ne JacHyp then
            printf "SPLITCODOMAIN idx=%o depth=%o\n", IDX, d; continue;
        end if;
        Append(~queue, <Norm2(Curve(S)), d+1>);
    end for;
end while;

printf "SEEDDONE idx=%o nodes=%o distinctgroups=%o newgroups=%o\n",
       IDX, nnodes, ngroups, nrecord;
quit;
