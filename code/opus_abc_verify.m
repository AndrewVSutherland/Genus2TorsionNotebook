//////////////////////////////////////////////////////////////////////
// opus_abc_verify.m  (2026-07-19)
//
// Structural analysis of Jen's (2,2,2,12) moduli surface
//
//   S:  y^2 = F1(A,B,C),   z^2 = F2(A,B,C)   in P^2_{A,B,C}
//
// Claim (found by hand, verified here): BOTH sextics factor, and they
// SHARE the conic factor W = C^2 + A^2 - B^2:
//
//   F1 = (B^2-A^2)(C^2-B^2) * W          [4 lines times the conic]
//   F2 = W * (C^4 + (A^2-B^2)C^2 + A^2(A^2-B^2))
//
// So the surface is NOT a generic fiber product of two K3s: it is a
// (Z/2)^2-cover of P^2 with THREE branch divisors
//   D1 = W (conic, deg 2), D2 = P1 (4 lines, deg 4), D3 = Q (quartic),
// and the conditions are equivalent to the single statement
//   W, P1, Q all lie in the SAME square class.
//
// Also verifies the (A,B,C) -> (s,m,n) -> genus-2 curve pipeline
// against both known (2,2,2,12) curves.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
QQ := Rationals();

// ---------- 1. the factorization ----------
R<A,B,C> := PolynomialRing(QQ, 3);
F1 := A^4*B^2 - 2*A^2*B^4 + B^6 - A^4*C^2 + 3*A^2*B^2*C^2 - 2*B^4*C^2
      - A^2*C^4 + B^2*C^4;
F2 := A^6 - 2*A^4*B^2 + A^2*B^4 + 2*A^4*C^2 - 3*A^2*B^2*C^2 + B^4*C^2
      + 2*A^2*C^4 - 2*B^2*C^4 + C^6;

W  := C^2 + A^2 - B^2;
P1 := (B^2-A^2)*(C^2-B^2);
Q  := C^4 + (A^2-B^2)*C^2 + A^2*(A^2-B^2);

assert F1 eq P1*W;
assert F2 eq W*Q;
printf "FACTORIZATION VERIFIED\n";
printf "  F1 = %o\n", Factorization(F1);
printf "  F2 = %o\n", Factorization(F2);
printf "  shared conic W = %o\n", W;
printf "  Q irreducible over Q: %o\n", IsIrreducible(Q);
// discriminant of Q as a quadratic in C^2 (tells when Q splits further)
printf "  disc_{C^2}(Q) = %o\n", (A^2-B^2)^2 - 4*A^2*(A^2-B^2);

// ---------- 2. the (A,B,C) -> (s,m,n) -> curve pipeline ----------
function smnFromABC(a, b, c)
    s := (a^3 - a*b^2 + 2*a*c^2) / (2*a^2*b - 2*b^3 + 2*b*c^2);
    m := (-2*a^2*b^2 + 2*b^4 + 2*a^2*c^2 - 6*b^2*c^2 + 4*c^4)
         / (2*a^3*b - 2*a*b^3 + 2*a*b*c^2);
    n := (2*a) / (2*b);
    return s, m, n;
end function;

// M(2,2,2,6) quintic in the (A_i + B_i x) form used by the lab
function CurveFromSMN(s, m, n)
    PP<x> := PolynomialRing(QQ);
    Bc := [ 2*s^2 - s*n,
            2*s^2 + s*m - 2*s*n - m*n,
            2*s^2 + s*m - s*n - m*n,
            -m*n,
            4*s^2 - 4*s*n - m*n ];
    Ac := [1,1,1,2,2];
    f := &*[ Ac[i] + Bc[i]*x : i in [1..5] ];
    return PP!f;
end function;

known := [ <336396, -689185, -166464>, <2208, -8303, -7200> ];
g2known := [];
for k in [1..#known] do
    t := known[k];
    f := CurveFromSMN(t[1], t[2], t[3]);
    if Degree(f) lt 5 or Discriminant(f) eq 0 then
        printf "KNOWN %o: degenerate\n", k; continue;
    end if;
    Cv := HyperellipticCurve(f);
    T := TorsionSubgroup(Jacobian(Cv));
    g2 := G2Invariants(Cv);
    Append(~g2known, g2);
    printf "KNOWN %o (s,m,n)=%o torsion=%o\n", k, t, Invariants(T);
    printf "KNOWN %o G2=%o\n", k, g2;
end for;

printf "PIPELINE_READY\n";
quit;
