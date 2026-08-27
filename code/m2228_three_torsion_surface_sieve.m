//////////////////////////////////////////////////////////////////////
//  Full-surface M(2,2,2,8) tests for rational 3-torsion.
//
//  This complements m2228_three_torsion_equations.m.  It has two modes:
//
//    finite:
//      Over finite fields, precompute all quintic coefficient tuples
//      admitting the triple-contact 3-torsion equations, then scan all
//      good affine K3 points s2(a,b,c,d)^2 = 4abcd.
//
//    exact_file:
//      Read primitive integer K3 curve tuples [a,b,c,d], keep only those
//      with forced bad reduction at p=13, apply good-prime point-count
//      bounds for rational 3-torsion, and run exact TorsionSubgroup only
//      on survivors.
//
//  Typical runs from torsion_jac:
//
//      magma -b mode:="finite" code/m2228_three_torsion_surface_sieve.m
//
//      magma -b mode:="exact_file" \
//          tuple_file:="data/surface_tuples_B2000.txt" \
//          code/m2228_three_torsion_surface_sieve.m
//////////////////////////////////////////////////////////////////////

if not assigned mode then
    mode := "finite";
end if;

Q := Rationals();
PX<X> := PolynomialRing(Q);

finite_primes := [5,7,11,13,17,19,23,29,31,37];
bound_primes := [
    5,7,11,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,
    79,83,89,97,101,103,107,109,113,127,131,137,139,149,
    151,157,163,167,173,179,181,191,193,197,199
];

if assigned tuple_file then
    exact_tuple_file := tuple_file;
else
    exact_tuple_file := "data/surface_tuples_B2000.txt";
end if;


function ContactCoefficientKeys(p)
    F := GF(p);
    PF<x> := PolynomialRing(F);
    keys := {};
    witnesses := AssociativeArray();
    raw_count := 0;

    for L in F do
        if L eq 0 then
            continue;
        end if;
        M := L^2;

        for U in F do
            for v in F do
                if v eq 0 or U^2 - 4*v^2 eq 0 then
                    continue;
                end if;

                q := x^2 + U*x + v^2;
                for e1 in F do
                    A := 2*M*e1 + 6*(U^2+v^2) - (M+3*U)^2;
                    e2 := ((M+3*U)*A + 8*v^3 - 4*U^3 - 24*U*v^2)/(4*M);
                    e3 := (A^2 + 16*(M+3*U)*v^3
                           - 48*(U^2*v^2+v^4))/(16*M);
                    e4 := (A*v^3 - 6*U*v^4)/(2*M);
                    if e4 eq 0 then
                        continue;
                    end if;

                    f := x^5 + e1*x^4 + e2*x^3 + e3*x^2 + e4*x;
                    if Degree(GCD(q, f)) gt 0 then
                        continue;
                    end if;

                    key := <Integers()!e1, Integers()!e2,
                            Integers()!e3, Integers()!e4>;
                    Include(~keys, key);
                    raw_count +:= 1;
                    if not IsDefined(witnesses, key) then
                        witnesses[key] := <Integers()!L, Integers()!U, Integers()!v>;
                    end if;
                end for;
            end for;
        end for;
    end for;

    return keys, witnesses, raw_count;
end function;

function ScanK3FiniteField(p, keys, witnesses)
    F := GF(p);
    good := 0;
    contact := 0;
    singular := 0;
    zero_boundary := 0;
    curve_keys := {};
    contact_curve_keys := {};
    samples := [];

    for a in F do
        for b in F do
            for c in F do
                for d in F do
                    if a*b*c*d eq 0 then
                        zero_boundary +:= 1;
                        continue;
                    end if;

                    s2 := a*b + a*c + a*d + b*c + b*d + c*d;
                    if s2^2 ne 4*a*b*c*d then
                        continue;
                    end if;

                    A := a^2; B := b^2; C := c^2; D := d^2;
                    if #(Set([A,B,C,D])) ne 4 then
                        singular +:= 1;
                        continue;
                    end if;

                    e1 := A+B+C+D;
                    e2 := A*B + A*C + A*D + B*C + B*D + C*D;
                    e3 := A*B*C + A*B*D + A*C*D + B*C*D;
                    e4 := A*B*C*D;
                    key := <Integers()!e1, Integers()!e2,
                            Integers()!e3, Integers()!e4>;

                    Include(~curve_keys, key);
                    good +:= 1;

                    if key in keys then
                        contact +:= 1;
                        Include(~contact_curve_keys, key);
                        if #samples lt 5 then
                            Append(~samples, <Integers()!a, Integers()!b,
                                              Integers()!c, Integers()!d,
                                              key, witnesses[key]>);
                        end if;
                    end if;
                end for;
            end for;
        end for;
    end for;

    return good, contact, #curve_keys, #contact_curve_keys,
           singular, zero_boundary, samples;
end function;

procedure RunFiniteSurfaceSieve()
    print "Full finite-field K3 triple-contact sieve";
    for p in finite_primes do
        keys, witnesses, raw_count := ContactCoefficientKeys(p);
        good, contact, ncurves, ncontact, singular, zero_boundary, samples :=
            ScanK3FiniteField(p, keys, witnesses);

        print "p", p,
              "contact_keys", #keys,
              "raw_contact_param", raw_count,
              "good_k3_points", good,
              "contact_points", contact,
              "curvekeys", ncurves,
              "contact_curvekeys", ncontact,
              "singular", singular,
              "zero_boundary_seen", zero_boundary;
        print " samples", samples;
    end for;
end procedure;

function ReadTupleFile(filename)
    S := Read(filename);
    rows := Split(S, "\n");
    tuples := [];

    for raw in rows do
        row := raw;
        if #row gt 0 and row[#row] eq "\r" then
            row := row[1..#row-1];
        end if;
        if #row ge 2 and row[1] eq "[" and row[#row] eq "]" then
            body := row[2..#row-1];
            parts := Split(body, ",");
            tup := [ StringToInteger(part) : part in parts | #part gt 0 ];
            if #tup eq 4 then
                Append(~tuples, tup);
            end if;
        end if;
    end for;

    return tuples;
end function;

function CurvePolynomialFromTuple(tup)
    f := X;
    for a in tup do
        f *:= X + (Q!a)^2;
    end for;
    return f;
end function;

function IsBadAtPrime(tup, p)
    F := GF(p);
    squares := [];
    for a in tup do
        if F!a eq 0 then
            return true;
        end if;
        Append(~squares, (F!a)^2);
    end for;
    return #(Set(squares)) lt #squares;
end function;

function TorsionGcdBound(f, primes)
    C := HyperellipticCurve(f);
    gcd_bound := 0;
    used := [];

    for p in primes do
        try
            fp := ChangeRing(f, GF(p));
        catch e
            continue;
        end try;
        if Degree(fp) ne 5 or Discriminant(fp) eq 0 then
            continue;
        end if;

        try
            n := Integers()!#Jacobian(ChangeRing(C, GF(p)));
        catch e
            continue;
        end try;

        if gcd_bound eq 0 then
            gcd_bound := n;
        else
            gcd_bound := GCD(gcd_bound, n);
        end if;
        Append(~used, <p, n, gcd_bound>);

        if gcd_bound mod 3 ne 0 then
            return gcd_bound, used;
        end if;
    end for;

    return gcd_bound, used;
end function;

function IrreducibleFrobeniusCertificate(f)
    C := HyperellipticCurve(f);
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79] do
        try
            fp := ChangeRing(f, GF(p));
            if Discriminant(fp) eq 0 then
                continue;
            end if;
            Lp := LPolynomial(ChangeRing(C, GF(p)));
            fac := Factorization(Lp);
            if #fac eq 1 and fac[1][2] eq 1 and Degree(fac[1][1]) eq 4 then
                return true, p, Lp;
            end if;
        catch e
            continue;
        end try;
    end for;
    return false, 0, PX!0;
end function;

procedure RunExactTupleFile()
    tuples := ReadTupleFile(exact_tuple_file);
    print "Exact tuple-file 3-torsion check";
    print "tuple_file", exact_tuple_file, "tuples", #tuples;

    forced13 := 0;
    no3_bound := 0;
    survived_bound := 0;
    exact_tests := 0;
    hits := [];

    for tup in tuples do
        if not IsBadAtPrime(tup, 13) then
            continue;
        end if;
        forced13 +:= 1;

        f := CurvePolynomialFromTuple(tup);
        if Discriminant(f) eq 0 then
            continue;
        end if;

        gcd_bound, used := TorsionGcdBound(f, bound_primes);
        if gcd_bound mod 3 ne 0 then
            no3_bound +:= 1;
            continue;
        end if;

        survived_bound +:= 1;
        print "BOUND_SURVIVOR", tup, "gcd_bound", gcd_bound, "used", used;

        C := HyperellipticCurve(f);
        J := Jacobian(C);
        exact_tests +:= 1;
        try
            G, phi := TorsionSubgroup(J);
            inv := Invariants(G);
        catch e
            print " torsion_error", e`Object;
            continue;
        end try;

        if &or [ n mod 3 eq 0 : n in inv ] then
            simple, pcert, Lp := IrreducibleFrobeniusCertificate(f);
            Append(~hits, <tup, inv, simple, pcert, Lp>);
            print "HIT", tup, "torsion", inv,
                  "simple", simple, "prime", pcert, "L", Lp;
        else
            print " exact_no3", tup, "torsion", inv;
        end if;
    end for;

    print "SUMMARY forced13", forced13,
          "no3_bound", no3_bound,
          "survived_bound", survived_bound,
          "exact", exact_tests,
          "hits", #hits;
    for H in hits do
        print H;
    end for;
end procedure;

if mode eq "finite" then
    RunFiniteSurfaceSieve();
elif mode eq "exact_file" then
    RunExactTupleFile();
else
    error "Unknown mode. Use mode:=\"finite\" or mode:=\"exact_file\".";
end if;

quit;
