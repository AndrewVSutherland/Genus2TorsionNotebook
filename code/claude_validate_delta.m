// Item 1, Step 1: validate the 2-descent (delta) framework on sample tuples.
// Curve: y^2 = f = x(x+a^2)(x+b^2)(x+c^2)(x+d^2), roots R = [0,-a^2,-b^2,-c^2,-d^2].
// delta: J(Q)/2J(Q) -> prod_e Q*/Q*2, [U,V] |-> (U(e))_e with f'(e)-correction at
// support Weierstrass points; finite-support product only (infinity ignored) -- the
// point of this script is to TEST that convention against ground truth.
// Ground truth: abstract torsion group arithmetic (halvability of each class).

QQ := Rationals();
Px<x> := PolynomialRing(QQ);

// delta vector of a Mumford class [U,V] (U monic, deg <= 2), finite components only
// component at e: prod over finite support points P of (x(P) - e),
// with f'(e) replacing any factor where x(P) = e.
// For monic U: prod(x(P) - e) = (-1)^deg(U) * U(e).
deltavec := function(f, R, U)
    vec := [];
    for e in R do
        v := Evaluate(U, e);
        if v ne 0 then
            Append(~vec, (-1)^Degree(U) * v);
        else
            U1 := U div (x - e);
            v1 := Evaluate(U1, e);
            fp := Evaluate(Derivative(f), e);
            if v1 ne 0 then
                Append(~vec, fp * (-1)^Degree(U1) * v1);
            else
                Append(~vec, QQ!1); // (x-e)^2: contributes f'(e)^2 ~ 1
            end if;
        end if;
    end for;
    return vec;
end function;

trivial := func<vec | &and[IsSquare(v) : v in vec]>;
mulvec := func<v1, v2 | [v1[i]*v2[i] : i in [1..#v1]]>;

// integer sqrt of a rational square
ratsqrt := function(q)
    ok, r := IsSquare(q);
    assert ok;
    return r;
end function;

procedure testtuple(a, b, c, d, has2244)
    A0 := a^2; B0 := b^2; C0 := c^2; D0v := d^2;
    f := x*(x+A0)*(x+B0)*(x+C0)*(x+D0v);
    R := [QQ!0, -A0, -B0, -C0, -D0v];
    Cc := HyperellipticCurve(f);
    J := Jacobian(Cc);
    Agrp, m := TorsionSubgroup(J);
    inv := Invariants(Agrp);
    printf "tuple [%o,%o,%o,%o]  torsion %o\n", a, b, c, d, inv;

    // all abstract elements and their Jacobian images
    elts := [g : g in Agrp];
    imgs := [m(g) : g in elts];

    // ground-truth halvability of a Jacobian torsion point
    isHalvable := function(P)
        for i in [1..#elts] do
            if imgs[i] eq P then
                gg := elts[i];
                return exists(t){h : h in Agrp | 2*h eq gg};
            end if;
        end for;
        error "point not in torsion list";
    end function;

    // the 16 two-torsion points with their Mumford U-polys (canonical classes)
    // e_{ij} finite pairs: U=(x-e_i)(x-e_j); e_{i,inf}: U=(x-e_i); zero: U=1
    twotors := [* *];
    Append(~twotors, <Px!1, J!0>);
    for i in [1..5] do
        U := x - R[i];
        Append(~twotors, <U, J![U, Px!0]>);
    end for;
    for i in [1..5] do for j in [i+1..5] do
        U := (x - R[i])*(x - R[j]);
        Append(~twotors, <U, J![U, Px!0]>);
    end for; end for;
    assert #twotors eq 16;

    // sanity: my delta on 2-torsion matches ground truth divisibility of e_{ij}
    nbad := 0;
    for t in twotors do
        pred := trivial(deltavec(f, R, t[1]));
        truth := isHalvable(t[2]);
        if pred ne truth then
            nbad +:= 1;
            printf "  2TORS MISMATCH U=%o pred=%o truth=%o\n", t[1], pred, truth;
        end if;
    end for;
    printf "  2-torsion delta test: %o mismatches out of 16\n", nbad;

    // D0 = Zarhin half of T0 (always rational on this chart)
    s1 := a+b+c+d; s2 := a*b+a*c+a*d+b*c+b*d+c*d;
    s3 := a*b*c+a*b*d+a*c*d+b*c*d; s4 := a*b*c*d;
    UD0 := x^2 - s2*x + s4;
    VD0 := (s1*s2 - s3)*x - s1*s4;
    D0pt := J![UD0, VD0];
    assert 2*D0pt eq J![x, Px!0];

    classlist := [* <"D0", UD0, D0pt> *];

    if has2244 then
        // H_AB from halving.m formulas (verified doubling below)
        u0 := ratsqrt((A0-C0)*(A0-D0v));
        v0 := ratsqrt((B0-C0)*(B0-D0v));
        w0 := ratsqrt((A0-C0)*(B0-C0));
        t0 := ratsqrt((A0-D0v)*(B0-D0v));
        if u0*v0 ne w0*t0 then t0 := -t0; end if;
        assert u0*v0 eq w0*t0;
        rho := a/b; sigma := (A0-C0)/w0; tau := (A0-D0v)/t0;
        s1g := 1+rho+sigma+tau;
        s2g := rho+sigma+tau+rho*sigma+rho*tau+sigma*tau;
        s3g := rho*sigma+rho*tau+sigma*tau+rho*sigma*tau;
        s4g := rho*sigma*tau;
        delta0 := 1+s2g+s4g;
        UH := x^2 + ((2*A0 + s2g*(A0+B0) + 2*B0*s4g)/delta0)*x
                  + (A0^2 + A0*B0*s2g + B0^2*s4g)/delta0;
        Lam := s1g*s2g - s1g*s4g^2 + 3*s1g*s4g + s2g*s3g*s4g + 3*s3g*s4g - s3g;
        Mm := A0*s1g*s2g + 2*A0*s1g*s4g + A0*s3g*s4g - A0*s3g
              - B0*s1g*s4g^2 + B0*s1g*s4g + B0*s2g*s3g*s4g + 2*B0*s3g*s4g;
        VH := -(b*v0/delta0^2) * (Lam*x + Mm);
        Hpt := J![UH, VH];
        assert 2*Hpt eq J![(x+A0)*(x+B0), Px!0];
        Append(~classlist, <"H_AB", UH, Hpt>);
        // D0 + H_AB: delta = product (homomorphism); ground truth on the sum point
        Append(~classlist, <"D0+H_AB", Px!0, D0pt + Hpt>); // U unused, flag below
    end if;

    // main test: for each class family and each of the 16 twists,
    // prediction [delta(class)*delta(T) trivial] vs truth [class+T halvable]
    for cl in classlist do
        if cl[1] eq "D0+H_AB" then
            base := mulvec(deltavec(f, R, UD0), deltavec(f, R, classlist[2][2]));
        else
            base := deltavec(f, R, cl[2]);
        end if;
        nbad := 0; nhalv := 0;
        for t in twotors do
            pred := trivial(mulvec(base, deltavec(f, R, t[1])));
            truth := isHalvable(cl[3] + t[2]);
            if truth then nhalv +:= 1; end if;
            if pred ne truth then
                nbad +:= 1;
                printf "  MISMATCH %o twist U=%o pred=%o truth=%o\n", cl[1], t[1], pred, truth;
            end if;
        end for;
        printf "  class %o: mismatches %o/16, halvable twists %o/16\n", cl[1], nbad, nhalv;
    end for;
end procedure;

// samples: tor2244 tuples (have (c)-conditions; H_AB exists)
tor2244samples := [[17,28,32,68],[13,57,68,112],[13,59,62,91],[5,110,133,731],
                   [33,58,63,87],[23,247,266,361]];
for t in tor2244samples do
    testtuple(t[1], t[2], t[3], t[4], true);
end for;

// tor2228 tuples (no (c)-conditions generically; D0-family only)
tor2228samples := [[2,4,23,46],[1,55,99,125],[4,11,16,44]];
for t in tor2228samples do
    testtuple(t[1], t[2], t[3], t[4], false);
end for;

print "VALIDATION COMPLETE";
quit;
