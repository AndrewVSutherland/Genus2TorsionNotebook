//////////////////////////////////////////////////////////////////////
//  Boundary component analysis modulo obstructing primes for halving Q
//  in the symmetric M_1(8,2,2) family.
//
//  For each finite field prime p, each base residue (s,p0), and each
//  monic quadratic a=X^2+A*X+B, we test whether
//
//      f_{s,p0}(X) - L*(X+1)*a(X)^2
//
//  is a square of a polynomial.  This is the coefficient-equation
//  version of the halving cover and still makes sense on singular
//  boundary fibers.
//////////////////////////////////////////////////////////////////////

prime_list := [7, 11];
max_examples := 20;

for ell in prime_list do
    F := GF(ell);
    P<X> := PolynomialRing(F);

    labels := [
        "p=0",
        "s+1=0",
        "s+2=0",
        "p-s+1=0",
        "s-p+1=0",
        "delta=0",
        "qdisc=0",
        "collision=0"
    ];

    function BoundaryValues(s, p0)
        vals := [
            p0,
            s + 1,
            s + 2,
            p0 - s + 1,
            s - p0 + 1,
            s^2 - 4*p0,
            2*s^2 + 3*s + p0 + 1,
            s^3 - s^2*p0 + s^2 - 4*s*p0 - 4*p0
        ];
        return vals;
    end function;

    function FamilyPolynomial(s, p0)
        qtilde := -X^2 + (p0*s - s^2 + 2*p0 - s - 2)*X
                   - (s^2 - p0 + s + 1);
        return ((p0-s+1)*X^2 + (2-s)*X + 1)*((s+2)*X + 1)*qtilde;
    end function;

    point_count := 0;
    open_count := 0;
    boundary_count := 0;
    split_count := 0;
    nonsplit_count := 0;
    open_split_count := 0;
    open_nonsplit_count := 0;
    boundary_split_count := 0;
    boundary_nonsplit_count := 0;
    base_hit_set := {};
    open_base_set := {};
    boundary_base_set := {};
    label_counts := AssociativeArray();
    label_base_sets := AssociativeArray();
    examples := [];

    for lab in labels do
        label_counts[lab] := 0;
        label_base_sets[lab] := {};
    end for;

    for s in F do
        for p0 in F do
            bvals := BoundaryValues(s, p0);
            is_boundary := &or [v eq 0 : v in bvals];
            delta_square := IsSquare(s^2 - 4*p0);
            base_key := Integers()!s + ell*(Integers()!p0);

            f := FamilyPolynomial(s, p0);
            L := Coefficient(f, 5);

            for A in F do
                for B in F do
                    a := X^2 + A*X + B;
                    res := f - L*(X+1)*a^2;
                    ok := IsSquare(res);
                    if not ok then
                        continue;
                    end if;

                    point_count +:= 1;
                    Include(~base_hit_set, base_key);
                    if is_boundary then
                        boundary_count +:= 1;
                        Include(~boundary_base_set, base_key);
                        if delta_square then
                            boundary_split_count +:= 1;
                        else
                            boundary_nonsplit_count +:= 1;
                        end if;
                    else
                        open_count +:= 1;
                        Include(~open_base_set, base_key);
                        if delta_square then
                            open_split_count +:= 1;
                        else
                            open_nonsplit_count +:= 1;
                        end if;
                    end if;
                    if delta_square then
                        split_count +:= 1;
                    else
                        nonsplit_count +:= 1;
                    end if;

                    for i in [1..#labels] do
                        if bvals[i] eq 0 then
                            label_counts[labels[i]] +:= 1;
                            Include(~label_base_sets[labels[i]], base_key);
                        end if;
                    end for;

                    if #examples lt max_examples then
                        active := [labels[i] : i in [1..#labels] | bvals[i] eq 0];
                        Append(~examples, <s, p0, A, B, active, delta_square, res>);
                    end if;
                end for;
            end for;
        end for;
    end for;

    print "prime", ell;
    print "halving_square_points", point_count;
    print "base_residues_with_points", #base_hit_set;
    print "open_points", open_count, "open_base_residues", #open_base_set,
          "open_split", open_split_count, "open_nonsplit", open_nonsplit_count;
    print "boundary_points", boundary_count, "boundary_base_residues", #boundary_base_set,
          "boundary_split", boundary_split_count, "boundary_nonsplit", boundary_nonsplit_count;
    print "delta_square_points", split_count, "delta_nonsquare_points", nonsplit_count;
    for lab in labels do
        print "component", lab,
              "points", label_counts[lab],
              "base_residues", #label_base_sets[lab];
    end for;
    print "examples";
    for ex in examples do
        print ex;
    end for;
end for;

quit;
