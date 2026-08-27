// claude_222_laneB_seed_export.m — export the two known [2,22] RM witnesses
// in the [1,1,2,2] chart normalization x(x-1)*q1(x)*q2(x), together with a
// Mumford representative of the order-22 generator and its 11-part.  This is
// the anchor data for the Lane-B/B1 transversal-deformation machinery
// (notes/claude_order222_lanes_assessment_2026_07_24.md).
SetColumns(0);
SetMemoryLimit(3*10^9);
Q := Rationals(); P<x> := PolynomialRing(Q);

procedure ExportSeed(name, F)
    printf "SEED %o\n", name;
    fac := Factorization(F);
    unit := F div &*[ t[1]^t[2] : t in fac ];           // Factorization unit trap
    assert unit * &*[ t[1]^t[2] : t in fac ] eq F;
    rts := [ -Coefficient(t[1],0) : t in fac | Degree(t[1]) eq 1 ];
    assert #rts eq 2;
    r1 := rts[1]; r2 := rts[2];
    printf "  rational roots: %o, %o\n", r1, r2;
    // normalize (r1,r2) -> (0,1)
    G := Evaluate(F, r1 + (r2-r1)*x);
    gfac := Factorization(G);
    gunit := G div &*[ t[1]^t[2] : t in gfac ];
    assert gunit * &*[ t[1]^t[2] : t in gfac ] eq G;
    quads := [ t[1] : t in gfac | Degree(t[1]) eq 2 ];
    assert #quads eq 2 and IsDivisibleBy(G, x*(x-1));
    printf "  chart model: y^2 = c*x*(x-1)*q1*q2 with\n";
    printf "  c  = %o\n  q1 = %o\n  q2 = %o\n", gunit, quads[1], quads[2];
    // integral model of the chart curve and its torsion generator
    den := LCM([Denominator(c) : c in Coefficients(G)]);
    Gint := P![ c*den^2 : c in Coefficients(G) ];   // y -> den*y
    C := HyperellipticCurve(Gint);
    J := Jacobian(C);
    T, mp := TorsionSubgroup(J);
    inv := Invariants(T);
    printf "  chart torsion: %o\n", inv;
    assert inv eq [2,22];
    g22 := mp(T.2);                 // generator of the order-22 factor
    g11 := 2*g22;                   // its 11-part
    printf "  order-22 generator (Mumford, on y^2 = %o):\n", Gint;
    printf "    u22 = %o\n    v22 = %o\n", g22[1], g22[2];
    printf "  order-11 class:\n";
    printf "    u11 = %o\n    v11 = %o\n", g11[1], g11[2];
    assert 11*g11 eq J!0 and g11 ne J!0;
    printf "  checks: 22*g22 = 0: %o, 11*g22 <> 0: %o\n\n",
        22*g22 eq J!0, 11*g22 ne J!0;
end procedure;

// 19044.h.2 (completed square of y^2+(x^2+x)y = x^6-3x^5+9x^4-5x^3+12x^2-6x)
ExportSeed("19044.h.2 [RM sqrt5]", (x^2+x)^2 + 4*(x^6-3*x^5+9*x^4-5*x^3+12*x^2-6*x));
// corrected BLP C4 (addendum 2 of notes/claude_generic_222_2214_plan_2026_07_23.md)
ExportSeed("BLP C4corr [RM sqrt5]", x^6-18*x^5-4001*x^4-22524*x^3+859039*x^2-1926258*x-9043839);
print "ALL_DONE";
quit;
