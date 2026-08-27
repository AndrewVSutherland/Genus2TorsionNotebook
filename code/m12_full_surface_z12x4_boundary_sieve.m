//////////////////////////////////////////////////////////////////////
//  Boundary sieve for Z/12 x Z/4 on the full two-dimensional
//  extra-Weierstrass surface inside M(12).
//
//  We split T+1 by
//      a = (1-z^2)/(4*(r+1)).
//  An additional rational root u of Q4 = W/(T+1) gives an independent
//  rational 2-torsion candidate.  The good nonboundary mod-7 points
//  have no halves, so a rational example must reduce to the boundary
//  modulo 7.  This script records the mod-7 boundary strata and then
//  applies the same residue test at several primes before any exact
//  rational halving test.
//
//  Typical run from torsion_jac:
//      magma -b height:=30 code/m12_full_surface_z12x4_boundary_sieve.m
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 30;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;

if not assigned progress_interval then
    progress_interval := 0;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;

Q := Rationals();
P<x> := PolynomialRing(Q);
PX<X> := PolynomialRing(Q);

function IsDivisibleBy2Finite(J, D)
    G, phi := AbelianGroup(J);
    a := D @@ phi;
    coords := Eltseq(a);
    invs := Invariants(G);
    for i in [1..#coords] do
        if GCD(2, invs[i]) ne 1 and coords[i] mod 2 ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;

function PairCategoryModP(p, r, z)
    F := GF(p);
    Pp<t> := PolynomialRing(F);
    PXp<Xp> := PolynomialRing(F);

    if r eq -1 then
        return "r=-1", [];
    end if;
    if z^2 eq 1 then
        return "z=+-1", [];
    end if;
    if z eq 0 then
        return "z=0", [];
    end if;

    a := (1-z^2)/(4*(r+1));
    T := a*t^2 - t + r;
    h := (t-r)*(T+1);
    W := h^2 + 4*a*t^2*T*(T+1);
    if Degree(W) ne 6 or Discriminant(W) eq 0 then
        return "discW", [];
    end if;

    rootsT := Roots(T+1);
    if #rootsT lt 2 then
        return "badT", [];
    end if;

    Q4 := ExactQuotient(W, T+1);
    if Degree(Q4) lt 4 then
        return "u_infty", [];
    end if;

    rootsQ := [ rt[1] : rt in Roots(Q4) | rt[2] eq 1 ];
    if #rootsQ eq 0 then
        return "good_no_qroot", [];
    end if;

    div_data := [];
    nodiv_data := [];
    for wd in rootsT do
        w := wd[1];
        f5 := PXp!0;
        for i in [0..Degree(W)] do
            for j in [0..i] do
                f5 +:= Coefficient(W, i)*Binomial(i,j)*w^(i-j)*Xp^(6-j);
            end for;
        end for;
        if Degree(f5) ne 5 or Discriminant(f5) eq 0 then
            continue;
        end if;

        C := HyperellipticCurve(f5);
        J := Jacobian(C);
        for u in rootsQ do
            if u eq w then
                continue;
            end if;
            beta := 1/(u-w);
            Tbeta := J![Xp-beta, F!0];
            if Tbeta eq J!0 then
                continue;
            end if;
            if IsDivisibleBy2Finite(J, Tbeta) then
                Append(~div_data, <Integers()!w, Integers()!u,
                                   Integers()!beta, true>);
            else
                Append(~nodiv_data, <Integers()!w, Integers()!u,
                                     Integers()!beta, false>);
            end if;
        end for;
    end for;

    if #div_data gt 0 then
        return "good_qroot_halve", div_data;
    end if;
    return "good_qroot_nohalve", nodiv_data;
end function;

procedure PrintMod7Classification()
    p := 7;
    F := GF(p);
    categories := ["r=-1", "z=+-1", "z=0", "discW", "badT",
                   "u_infty", "good_no_qroot", "good_qroot_nohalve",
                   "good_qroot_halve"];
    counts := AssociativeArray();
    residues := AssociativeArray();
    details := AssociativeArray();
    for category in categories do
        counts[category] := 0;
        residues[category] := [];
        details[category] := [];
    end for;

    for r in F do
        for z in F do
            category, data := PairCategoryModP(p, r, z);
            counts[category] +:= 1;
            Append(~residues[category], <Integers()!r, Integers()!z>);
            if #data gt 0 then
                Append(~details[category], <Integers()!r, Integers()!z, data>);
            end if;
        end for;
    end for;

    print "mod 7 affine (r,z) classification";
    for category in categories do
        if counts[category] gt 0 then
            print " ", category, counts[category], residues[category];
            if #details[category] gt 0 then
                print "   details", details[category];
            end if;
        end if;
    end for;
end procedure;

function AllowedPairs(p)
    F := GF(p);
    allowed := {};
    stats := AssociativeArray();
    keys := ["boundary", "bad", "uinfty", "good_div",
             "excluded_noqroot", "excluded_nodiv"];
    for key in keys do
        stats[key] := 0;
    end for;

    for r in F do
        for z in F do
            key := <Integers()!r, Integers()!z>;
            category, data := PairCategoryModP(p, r, z);
            if category eq "r=-1" or category eq "z=+-1" or category eq "z=0" then
                Include(~allowed, key);
                stats["boundary"] +:= 1;
            elif category eq "discW" or category eq "badT" then
                Include(~allowed, key);
                stats["bad"] +:= 1;
            elif category eq "u_infty" then
                Include(~allowed, key);
                stats["uinfty"] +:= 1;
            elif category eq "good_qroot_halve" then
                Include(~allowed, key);
                stats["good_div"] +:= 1;
            elif category eq "good_no_qroot" then
                stats["excluded_noqroot"] +:= 1;
            else
                stats["excluded_nodiv"] +:= 1;
            end if;
        end for;
    end for;

    return allowed, stats;
end function;

function RationalParametersOfHeight(B)
    vals := [];
    seen := {};
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num, den) ne 1 then
                continue;
            end if;
            q := Q!num/den;
            key := Sprint(q);
            if key notin seen then
                Include(~seen, key);
                Append(~vals, q);
            end if;
        end for;
    end for;
    return vals;
end function;

function M12Data(a, r)
    T := a*x^2 - x + r;
    h := (x-r)*(T+1);
    W := h^2 + 4*a*x^2*T*(T+1);
    Q4 := ExactQuotient(W, T+1);
    return W, T, h, Q4;
end function;

function OddQuinticAtRoot(W, w)
    out := PX!0;
    for i in [0..Degree(W)] do
        for j in [0..i] do
            out +:= Coefficient(W, i)*Binomial(i,j)*w^(i-j)*X^(6-j);
        end for;
    end for;
    return out;
end function;

function IntegralModelPolynomial(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return PX!(L^2*f), L;
end function;

function IrreducibleFrobeniusCertificate(f)
    C := HyperellipticCurve(f);
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47] do
        try
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

print "M(12) full extra-Weierstrass Z/12 x Z/4 boundary sieve";
print "height", height;
PrintMod7Classification();

primes := [7,11,13,17,19,23,29,31,37,41,43];
allowed := AssociativeArray(Integers());
for p in primes do
    A, S := AllowedPairs(p);
    allowed[p] := A;
    print "p", p, "allowed", #A, "boundary", S["boundary"],
          "bad", S["bad"], "uinfty", S["uinfty"],
          "good_div", S["good_div"],
          "excluded_noqroot", S["excluded_noqroot"],
          "excluded_nodiv", S["excluded_nodiv"];
end for;

params := RationalParametersOfHeight(height);
reason := AssociativeArray();
reason_keys := ["param_boundary", "W_bad", "T_bad", "Q4_no_root",
                "odd_bad", "order12_bad", "no_independent",
                "exact_tested"];
for key in reason_keys do
    reason[key] := 0;
end for;

survivors := 0;
exact_tests := 0;
hits := [];
total := 0;

for r in params do
    for z in params do
        total +:= 1;

        pass := true;
        for p in primes do
            if Denominator(r) mod p eq 0 or Denominator(z) mod p eq 0 then
                continue;
            end if;
            rr := Integers()!((Integers(p)!Numerator(r))
                              *(Integers(p)!Denominator(r))^-1);
            zz := Integers()!((Integers(p)!Numerator(z))
                              *(Integers(p)!Denominator(z))^-1);
            if <rr, zz> notin allowed[p] then
                pass := false;
                break;
            end if;
        end for;
        if not pass then
            continue;
        end if;
        survivors +:= 1;

        if r eq -1 or z^2 eq 1 or z eq 0 then
            reason["param_boundary"] +:= 1;
            continue;
        end if;

        a := (1-z^2)/(4*(r+1));
        if a eq 0 then
            reason["param_boundary"] +:= 1;
            continue;
        end if;

        W, T, h, Q4 := M12Data(a, r);
        if Degree(W) ne 6 or Discriminant(W) eq 0 then
            reason["W_bad"] +:= 1;
            continue;
        end if;

        rootsT := Roots(T+1);
        if #rootsT lt 2 then
            reason["T_bad"] +:= 1;
            continue;
        end if;

        rootsQ := [ rt[1] : rt in Roots(Q4) | rt[2] eq 1 ];
        if #rootsQ eq 0 then
            reason["Q4_no_root"] +:= 1;
            continue;
        end if;

        local_exact := false;
        local_odd := false;
        local_order := false;
        local_independent := false;

        for wd in rootsT do
            w := wd[1];
            if w eq 0 then
                continue;
            end if;

            f5 := OddQuinticAtRoot(W, w);
            if Degree(f5) ne 5 or Discriminant(f5) eq 0 then
                continue;
            end if;
            local_odd := true;

            fI, L := IntegralModelPolynomial(f5);
            if Discriminant(fI) eq 0 then
                continue;
            end if;

            Y0 := Evaluate(h, 0);
            Xp := -1/w;
            Yp := Y0*Xp^3;
            C := HyperellipticCurve(fI);
            J := Jacobian(C);
            try
                D := J![X-Xp, L*Yp];
                ord := Order(D);
            catch e
                continue;
            end try;
            if ord ne 12 then
                continue;
            end if;
            local_order := true;
            Tdiv := 6*D;

            for u in rootsQ do
                if u eq w then
                    continue;
                end if;
                beta := 1/(u-w);
                Tbeta := J![X-beta, Q!0];
                if Tbeta eq J!0 or Tbeta eq Tdiv then
                    continue;
                end if;
                local_independent := true;
                exact_tests +:= 1;
                local_exact := true;

                if IsDivisibleBy(Tbeta, 2) then
                    simple, pcert, Lp := IrreducibleFrobeniusCertificate(f5);
                    Append(~hits, <r,z,a,w,u,beta,f5,simple,pcert,Lp>);
                    print "HIT";
                    print "  r,z,a", r, z, a;
                    print "  w,u,beta", w, u, beta;
                    print "  simple", simple, "prime", pcert, "L", Lp;
                    print "  f5", f5;
                end if;
            end for;
        end for;

        if local_exact then
            reason["exact_tested"] +:= 1;
        elif not local_odd then
            reason["odd_bad"] +:= 1;
        elif not local_order then
            reason["order12_bad"] +:= 1;
        elif not local_independent then
            reason["no_independent"] +:= 1;
        end if;

        if progress_interval gt 0 and total mod progress_interval eq 0 then
            print "progress", total, "survivors", survivors,
                  "exact_tests", exact_tests, "hits", #hits;
        end if;
    end for;
end for;

print "DONE height", height, "total", total, "survivors", survivors,
      "exact", exact_tests, "hits", #hits;
for key in reason_keys do
    print key, reason[key];
end for;
for H in hits do
    print H;
end for;
quit;
