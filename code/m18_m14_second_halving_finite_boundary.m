
//////////////////////////////////////////////////////////////////////
//  Finite-prime and boundary diagnostics for the second-halving covers
//  inside the M_1(8,4) / [4,8] family.
//
//  Good affine chart:
//    [4,16] requires Tx divisible by 2 and P_R divisible by 2.
//    [8,8]  requires Tx divisible by 4.
//
//  Boundary chart:
//    If a prime has no good affine target residues, a rational point must
//    reduce to one of the base boundary factors listed below.  This script
//    counts those boundary residues by factor and, in boundary mode, prints
//    intersections among the main factors.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned mode then
    mode := "finite";
end if;

prime_list := [3,5,7,11,13,17,19,23,29,31,37,41,43];
if assigned primes then
    if Type(primes) eq MonStgElt then
        prime_list := [StringToInteger(s) : s in Split(primes, ",")];
    else
        prime_list := primes;
    end if;
end if;

max_samples := 12;
if assigned samples then
    if Type(samples) eq MonStgElt then
        max_samples := StringToInteger(samples);
    else
        max_samples := samples;
    end if;
end if;

Z := Integers();

function GoodHyperellipticPolynomial(f)
    return Degree(f) eq 5 and Discriminant(f) ne 0;
end function;

function DivisibleByNFinite(J, D, n)
    G, phi := AbelianGroup(J);
    a := D @@ phi;
    coords := Eltseq(a);
    invs := Invariants(G);
    for i in [1..#coords] do
        g := GCD(n, invs[i]);
        if coords[i] mod g ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;

function FamilyPolynomialFinite(F, R, w)
    P<X> := PolynomialRing(F);
    t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
    A := X^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*X + R^4;
    B := (R + 2 + 4*t)*X^2 + (R^2 + 4*R + 1 + 8*t)*X
         + (2*R^2 + R + 4*t);
    return X*A*B, t, A, B;
end function;

function YR(F, R, w)
    Qfac := R^2 - (F!1/2)*R*w^2 + (F!1/2)*R - w^2;
    return -2*R*(R-1)^2*Qfac/(w^2-1);
end function;

function BoundaryLabels(F, R, w)
    labels := [];
    if R eq 0 then Append(~labels, "R"); end if;
    if w eq 0 then Append(~labels, "w"); end if;
    if w - 1 eq 0 then Append(~labels, "w-1"); end if;
    if w + 1 eq 0 then Append(~labels, "w+1"); end if;
    if R - 1 eq 0 then Append(~labels, "R-1"); end if;
    if R + 1 eq 0 then Append(~labels, "R+1"); end if;
    if R - w eq 0 then Append(~labels, "R-w"); end if;
    if R + w eq 0 then Append(~labels, "R+w"); end if;
    if R*w - 3*R + 3*w - 1 eq 0 then Append(~labels, "Lplus"); end if;
    if R*w + 3*R + 3*w + 1 eq 0 then Append(~labels, "Lminus"); end if;
    if 2*R^2 - R*w^2 + R - 2*w^2 eq 0 then Append(~labels, "Q"); end if;
    if R^4 - 2*R^3 + R^2*w^2 - R^2 + 2*R*w^2 - w^2 eq 0 then
        Append(~labels, "Quartic");
    end if;
    return labels;
end function;

function BoundaryKey(labels)
    if #labels eq 0 then
        return "none";
    end if;
    return Join(Sort(labels), "&");
end function;

procedure Increment(~A, key)
    if IsDefined(A, key) then
        A[key] +:= 1;
    else
        A[key] := 1;
    end if;
end procedure;

procedure AddSample(~A, key, pt)
    if not IsDefined(A, key) then
        A[key] := [];
    end if;
    if #A[key] lt max_samples then
        Append(~A[key], pt);
    end if;
end procedure;

procedure FiniteOpen()
    print "FINITE_OPEN second-halving covers on good affine M_1(8,4) chart";
    print "primes", prime_list;
    print "p=2 is omitted: characteristic 2 is outside the good 2-primary reduction setup.";
    for p in prime_list do
        if p eq 2 then
            print "p", p, "skipped_char2";
            continue;
        end if;
        F := GF(p);
        P<X> := PolynomialRing(F);
        checked := 0;
        base_boundary := 0;
        good := 0;
        first := 0;
        pr2 := 0;
        tx4 := 0;
        target416 := 0;
        target88 := 0;
        both := 0;
        sample416 := [];
        sample88 := [];
        sampleFirst := [];

        for R in F do
            for w in F do
                checked +:= 1;
                labels := BoundaryLabels(F, R, w);
                if #labels ne 0 then
                    base_boundary +:= 1;
                    continue;
                end if;

                f, t, A, B := FamilyPolynomialFinite(F, R, w);
                if not GoodHyperellipticPolynomial(f) then
                    base_boundary +:= 1;
                    continue;
                end if;

                yR := YR(F, R, w);
                if yR^2 ne Evaluate(f, -R) then
                    print "WARNING P_R formula failed", p, R, w;
                    continue;
                end if;

                C := HyperellipticCurve(f);
                J := Jacobian(C);
                Tx := J![X, F!0];
                PR := J![X + R, yR];

                has_first := DivisibleByNFinite(J, Tx, 2);
                has_pr2 := DivisibleByNFinite(J, PR, 2);
                has_tx4 := DivisibleByNFinite(J, Tx, 4);

                good +:= 1;
                if has_first then
                    first +:= 1;
                    if #sampleFirst lt max_samples then
                        Append(~sampleFirst, <Z!R, Z!w>);
                    end if;
                end if;
                if has_pr2 then pr2 +:= 1; end if;
                if has_tx4 then tx4 +:= 1; end if;
                if has_first and has_pr2 then
                    target416 +:= 1;
                    if #sample416 lt max_samples then
                        Append(~sample416, <Z!R, Z!w>);
                    end if;
                end if;
                if has_tx4 then
                    target88 +:= 1;
                    if #sample88 lt max_samples then
                        Append(~sample88, <Z!R, Z!w>);
                    end if;
                end if;
                if has_first and has_pr2 and has_tx4 then
                    both +:= 1;
                end if;
            end for;
        end for;

        print "p", p,
              "checked", checked,
              "base_boundary", base_boundary,
              "good_open", good,
              "first_Tx_half", first,
              "PR_half", pr2,
              "target416", target416,
              "Tx_fourdiv_target88", target88,
              "both", both;
        print "  samples_first", sampleFirst;
        print "  samples_416", sample416;
        print "  samples_88", sample88;
    end for;
end procedure;

procedure BoundarySummary()
    print "BOUNDARY_SUMMARY for base bad locus of M_1(8,4) chart";
    print "primes", prime_list;
    for p in prime_list do
        if p eq 2 then
            print "p", p, "char2_boundary_all";
            continue;
        end if;
        F := GF(p);
        label_counts := AssociativeArray();
        exact_counts := AssociativeArray();
        exact_samples := AssociativeArray();
        boundary_total := 0;
        open_good := 0;
        open_target416 := 0;
        open_target88 := 0;

        for R in F do
            for w in F do
                labels := BoundaryLabels(F, R, w);
                is_open_good := false;
                if #labels eq 0 then
                    f, t, A, B := FamilyPolynomialFinite(F, R, w);
                    if GoodHyperellipticPolynomial(f) then
                        is_open_good := true;
                    else
                        labels := ["disc_extra"];
                    end if;
                end if;

                if is_open_good then
                    open_good +:= 1;
                    P<X> := PolynomialRing(F);
                    yR := YR(F, R, w);
                    C := HyperellipticCurve(f);
                    J := Jacobian(C);
                    Tx := J![X, F!0];
                    PR := J![X + R, yR];
                    has_first := DivisibleByNFinite(J, Tx, 2);
                    has_pr2 := DivisibleByNFinite(J, PR, 2);
                    has_tx4 := DivisibleByNFinite(J, Tx, 4);
                    if has_first and has_pr2 then open_target416 +:= 1; end if;
                    if has_tx4 then open_target88 +:= 1; end if;
                else
                    boundary_total +:= 1;
                    key := BoundaryKey(labels);
                    Increment(~exact_counts, key);
                    AddSample(~exact_samples, key, <Z!R, Z!w>);
                    for lab in labels do
                        Increment(~label_counts, lab);
                    end for;
                end if;
            end for;
        end for;

        print "p", p,
              "open_good", open_good,
              "open_target416", open_target416,
              "open_target88", open_target88,
              "boundary_total", boundary_total;
        print "  factor_counts";
        for key in Sort([k : k in Keys(label_counts)]) do
            print "   ", key, label_counts[key];
        end for;
        print "  exact_boundary_strata";
        for key in Sort([k : k in Keys(exact_counts)]) do
            print "   ", key, exact_counts[key], "samples", exact_samples[key];
        end for;
    end for;
end procedure;

if mode eq "finite" then
    FiniteOpen();
elif mode eq "boundary" then
    BoundarySummary();
elif mode eq "all" then
    FiniteOpen();
    BoundarySummary();
else
    print "Unknown mode", mode;
end if;

quit;
