// Complete enumeration of lines and conics on S:
//   Q2: x1^2+x2^2+x3^2-x4^2-x5^2 = 0,  F4: x1^4+..+x3^4-x4^4-x5^4 = 0.
// A nondegenerate rational line/conic => infinitely many [2,2,2,12] fibers.
SetColumns(0);
SetMemoryLimit(8*10^9);
Q := Rationals();

// ---------- A) LINES: charts over pivot pairs ----------
// line = row span of R1, R2 where columns (i,j) form the identity.
print "A) lines on S, full Grassmannian sweep";
npts := 0;
for i in [1..4] do for j in [i+1..5] do
    rest := [ k : k in [1..5] | k ne i and k ne j ];
    S6<a1,a2,a3,b1,b2,b3> := PolynomialRing(Q, 6, "grevlex");
    R1 := [ S6!0 : k in [1..5] ]; R2 := [ S6!0 : k in [1..5] ];
    R1[i] := 1; R2[j] := 1;
    for t in [1..3] do R1[rest[t]] := [a1,a2,a3][t]; R2[rest[t]] := [b1,b2,b3][t]; end for;
    PU<u,w> := PolynomialRing(S6, 2);
    X := [ u*R1[k] + w*R2[k] : k in [1..5] ];
    Q2 := X[1]^2+X[2]^2+X[3]^2-X[4]^2-X[5]^2;
    F4 := X[1]^4+X[2]^4+X[3]^4-X[4]^4-X[5]^4;
    eqs := Coefficients(Q2) cat Coefficients(F4);
    I := ideal<S6 | eqs>;
    d := Dimension(I);
    if d lt 0 then continue; end if;
    printf "  chart (%o,%o): dim %o", i, j, d;
    if d eq 0 then
        V := Variety(I);
        printf "  #pts %o", #V;
        for v in V do
            RR1 := [ Q!0 : k in [1..5] ]; RR2 := [ Q!0 : k in [1..5] ];
            RR1[i] := 1; RR2[j] := 1;
            for t in [1..3] do RR1[rest[t]] := v[t]; RR2[rest[t]] := v[3+t]; end for;
            gen := [ 2*RR1[k] + 3*RR2[k] : k in [1..5] ];
            sq := [ x^2 : x in gen ];
            if &and[ x ne 0 : x in gen ] and #Set(sq) eq 5 then
                npts +:= 1;
                printf "\n  NONDEG LINE gen=%o", gen;
            end if;
        end for;
    end if;
    printf "\n";
end for; end for;
printf "A) nondegenerate lines: %o\n", npts;

// ---------- B) CONICS: planes where Q2|_plane divides F4|_plane ----------
// plane = row span of R1,R2,R3, pivots (i,j,k); F4|_P = Q2|_P * G2 with
// unknown quadratic form G2 (6 coeffs, linear in the system).
print "B) conics on S, plane sweep";
nc := 0;
for i in [1..3] do for j in [i+1..4] do for k in [j+1..5] do
    rest := [ t : t in [1..5] | t ne i and t ne j and t ne k ];
    S12<a1,a2,b1,b2,c1,c2,g1,g2,g3,g4,g5,g6> := PolynomialRing(Q, 12, "grevlex");
    R1 := [ S12!0 : t in [1..5] ]; R2 := [ S12!0 : t in [1..5] ]; R3 := [ S12!0 : t in [1..5] ];
    R1[i] := 1; R2[j] := 1; R3[k] := 1;
    R1[rest[1]] := a1; R1[rest[2]] := a2;
    R2[rest[1]] := b1; R2[rest[2]] := b2;
    R3[rest[1]] := c1; R3[rest[2]] := c2;
    PU<u,v,w> := PolynomialRing(S12, 3);
    X := [ u*R1[t] + v*R2[t] + w*R3[t] : t in [1..5] ];
    Q2 := X[1]^2+X[2]^2+X[3]^2-X[4]^2-X[5]^2;
    F4 := X[1]^4+X[2]^4+X[3]^4-X[4]^4-X[5]^4;
    G2 := g1*u^2+g2*v^2+g3*w^2+g4*u*v+g5*u*w+g6*v*w;
    D := F4 - Q2*G2;
    eqs := Coefficients(D);
    I := ideal<S12 | eqs>;
    d := Dimension(I);
    printf "  chart (%o,%o,%o): dim %o\n", i, j, k, d;
    if d eq 0 then
        V := Variety(I);
        printf "    #pts %o\n", #V;
        for vv in V do
            RR := [[ Q!0 : t in [1..5]], [ Q!0 : t in [1..5]], [ Q!0 : t in [1..5]]];
            RR[1][i] := 1; RR[2][j] := 1; RR[3][k] := 1;
            RR[1][rest[1]] := vv[1]; RR[1][rest[2]] := vv[2];
            RR[2][rest[1]] := vv[3]; RR[2][rest[2]] := vv[4];
            RR[3][rest[1]] := vv[5]; RR[3][rest[2]] := vv[6];
            // conic = {Q2|_P = 0} in this plane; sample a rational point on it if any:
            printf "    CONIC plane R1=%o R2=%o R3=%o\n", RR[1], RR[2], RR[3];
            nc +:= 1;
        end for;
    end if;
end for; end for; end for;
printf "B) conic-planes found: %o\n", nc;
printf "CURVES_DONE\n";
quit;
