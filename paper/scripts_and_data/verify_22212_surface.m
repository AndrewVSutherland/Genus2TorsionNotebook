// verify_22212_surface.m -- the "Geometry of S" and "How we found S^o"
// claims of the paper (Section 3).
//
// S c P^4 is the complete intersection  a^2+b^2+c^2 = u^2+v^2,
// a^4+b^4+c^4 = u^4+v^4.  This script verifies:
//   (1) S is a surface whose singular subscheme is 0-dimensional of
//       degree 36, with exactly 36 rational singular points, all with
//       coordinates in {0,1,-1};
//   (2) every singular point is an ordinary double point (node): at each
//       point the two local equations have proportional nonzero linear
//       parts, and the induced quadratic form on the 3-dimensional tangent
//       space has full rank 3;
//   (3) under G = {+-1}^5/<-1> x (S3 x S2) the 36 nodes are the orbit of
//       [1:1:0:0:0] (12 points) together with the orbit of [1:1:0:1:1]
//       (24 points);
//   (4) the displayed map [a:b:c:u:v] |-> [a, u, b, a^2b^2c/(u^2-c^2),
//       v(c^2-v^2)] sends S into the surface Z:
//         x3^2 = (B^2-C^2)(B^2-A^2)(B^2-A^2-C^2),
//         x4^2 = (B^2-A^2-C^2)(B^2A^2+B^2C^2-A^4-A^2C^2-C^4)
//       with (A,B,C) = (a,u,b): both equations, after clearing the
//       denominator, lie in the ideal of S.  (This verifies the map is
//       well defined S --> Z; birationality is the discovery narrative.)
// Run from this directory:  magma -b verify_22212_surface.m
SetColumns(0);
SetMemoryLimit(8*10^9);
Q := Rationals();
t0 := Cputime();

P4<a,b,c,u,v> := ProjectiveSpace(Q, 4);
S := Scheme(P4, [a^2+b^2+c^2-u^2-v^2, a^4+b^4+c^4-u^4-v^4]);
assert Dimension(S) eq 2;

// ---- (1) the singular locus ----
Sg := SingularSubscheme(S);
assert Dimension(Sg) eq 0;
Sgr := ReducedSubscheme(Sg);
assert Degree(Sgr) eq 36;                 // 36 geometric singular points
pts := RationalPoints(Sgr);
assert #pts eq 36;                        // ... and all 36 are rational
for pt in pts do
    assert &and[ co in {Q!0, Q!1, Q!-1} : co in Coordinates(pt) ];
end for;
printf "(1) Sing(S) = 36 rational points, coordinates in {0,1,-1}\n";

// ---- (2) every singular point is a node ----
// At a singular point of the (2,4) complete intersection the two linear
// parts are proportional (rank 1); eliminating the smooth direction, the
// residual quadratic form on the tangent 3-space must be nondegenerate.
R4 := PolynomialRing(Q, 4);
function IsNodeAt(co)
    j := Min([ i : i in [1..5] | co[i] ne 0 ]);        // affine chart
    k := 0;
    sub5 := [];
    for i in [1..5] do
        if i eq j then
            Append(~sub5, R4!co[i]);
        else
            k +:= 1;
            Append(~sub5, R4!co[i] + R4.k);
        end if;
    end for;
    aa, bb, cc, uu, vv := Explode(sub5);
    F1 := aa^2 + bb^2 + cc^2 - uu^2 - vv^2;
    F2 := aa^4 + bb^4 + cc^4 - uu^4 - vv^4;
    l1 := [ MonomialCoefficient(F1, R4.i) : i in [1..4] ];
    l2 := [ MonomialCoefficient(F2, R4.i) : i in [1..4] ];
    M12 := Matrix(Q, 2, 4, l1 cat l2);
    if Rank(M12) ne 1 then return false; end if;       // 0: worse sing; 2: smooth
    if exists(i){ i : i in [1..4] | l1[i] ne 0 } then
        h := F1; lh := l1; lam := l2[i]/l1[i]; g := F2 - lam*F1;
    else
        i := [ k : k in [1..4] | l2[k] ne 0 ][1];
        h := F2; lh := l2; g := F1;
    end if;
    // quadratic part of g as a symmetric matrix
    Qg := ZeroMatrix(Q, 4, 4);
    for r in [1..4] do
        Qg[r][r] := MonomialCoefficient(g, R4.r^2);
        for s in [r+1..4] do
            Qg[r][s] := MonomialCoefficient(g, R4.r*R4.s)/2;
            Qg[s][r] := Qg[r][s];
        end for;
    end for;
    // basis of ker(lh) and the restricted form
    K := KernelMatrix(Matrix(Q, 4, 1, lh));
    assert Nrows(K) eq 3;
    return Rank(K * Qg * Transpose(K)) eq 3;
end function;
for pt in pts do
    assert IsNodeAt(Coordinates(pt));
end for;
printf "(2) all 36 singular points are ordinary double points\n";

// ---- (3) the two G-orbits ----
function Normalize(co)   // projective representative over {+-1}
    j := Min([ i : i in [1..5] | co[i] ne 0 ]);
    s := co[j] gt 0 select 1 else -1;
    return [ s*t : t in co ];
end function;
function Orbit(seed)
    S3 := [ p : p in Permutations({1,2,3}) ];
    orb := {};
    for sgn in CartesianPower({1,-1}, 5) do
        for p3 in S3, swap in [false, true] do
            im := [ sgn[1]*seed[p3[1]], sgn[2]*seed[p3[2]], sgn[3]*seed[p3[3]],
                    sgn[4]*seed[swap select 5 else 4],
                    sgn[5]*seed[swap select 4 else 5] ];
            Include(~orb, Normalize(im));
        end for;
    end for;
    return orb;
end function;
// NOTE: an earlier manuscript draft printed the 12-point orbit
// representative as [1:1:0:0:0]; that point does not lie on S
// (2 = 0 in the quadric).  The correct representative is [1:0:0:1:0].
O1 := Orbit([1,0,0,1,0]);
O2 := Orbit([1,1,0,1,1]);
assert #O1 eq 12 and #O2 eq 24 and #(O1 meet O2) eq 0;
singset := { Normalize([ Integers()!t : t in Coordinates(pt) ]) : pt in pts };
assert singset eq O1 join O2;
printf "(3) Sing(S) = orbit of [1:0:0:1:0] (12) + orbit of [1:1:0:1:1] (24)\n";

// ---- (4) the map S --> Z ----
R5<A0,B0,C0,U0,V0> := PolynomialRing(Q, 5);
IS := ideal< R5 | A0^2+B0^2+C0^2-U0^2-V0^2, A0^4+B0^4+C0^4-U0^4-V0^4 >;
// (A,B,C,x3,x4) = (a, u, b, a^2b^2c/(u^2-c^2), v(c^2-v^2))
A := A0; B := U0; C := B0;
x3num := A0^2*B0^2*C0;  x3den := U0^2 - C0^2;
x4 := V0*(C0^2 - V0^2);
E1 := x3num^2 - (B^2-C^2)*(B^2-A^2)*(B^2-A^2-C^2)*x3den^2;
E2 := x4^2 - (B^2-A^2-C^2)*(B^2*A^2+B^2*C^2-A^4-A^2*C^2-C^4);
assert NormalForm(E1, IS) eq 0;
assert NormalForm(E2, IS) eq 0;
printf "(4) the displayed map sends S into Z (both equations lie in I_S)\n";

printf "SURFACE GEOMETRY VERIFIED (%.1o s)\n", Cputime() - t0;
quit;
