//////////////////////////////////////////////////////////////////////
//  Quintic-contact family with rational 5-torsion.
//
//  Let
//      h(x) = 1 + a*x + b*x^2,
//      f(x) = h(x)^2 - (1+a+b)^2*x^5.
//
//  Then P=(0,1) lies on C: y^2=f(x), and
//
//      div(y-h(x)) = 5(P) - 5(infinity),
//
//  provided f is squarefree of degree 5.  Thus [P-infinity] is a
//  rational 5-torsion point.  The normalization x=1 as a root loses no
//  generality once the quintic has one rational root away from P.
//
//  The search below asks for f to split completely over Q, which adds
//  full rational 2-torsion and gives a subgroup containing
//      (Z/2)^3 x Z/10.
//
//  Typical runs from torsion_jac:
//      magma -b mode:="finite" code/m10_quintic_contact5_search.m
//      magma -b mode:="search" height:=30 code/m10_quintic_contact5_search.m
//      magma -b mode:="boundary" height:=80 code/m10_quintic_contact5_search.m
//////////////////////////////////////////////////////////////////////

if not assigned mode then
    mode := "finite";
end if;
if not assigned height then
    height := 30;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned max_hits then
    max_hits := 20;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

finite_primes := [3,5,7,11,13,17,19,23,29,31,37,41,43];
filter_primes := [3,5,7,11,13,17,19];
boundary_primes := [3,5,7,13,17];
boundary_extra_split_primes := [11,19,23,29,31,37,41,43];
boundary_filter_primes := boundary_primes cat boundary_extra_split_primes;

function RationalParametersOfHeight(B)
    vals := [];
    seen := {};
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num, den) ne 1 then
                continue;
            end if;
            r := Q!num/den;
            key := Sprint(r);
            if key notin seen then
                Include(~seen, key);
                Append(~vals, r);
            end if;
        end for;
    end for;
    return vals;
end function;

function FiveTorsionPolynomial(a, b)
    h := 1 + a*x + b*x^2;
    k := (1+a+b)^2;
    return h^2 - k*x^5, h, k;
end function;

function IsFullSplitPolynomial(f)
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        return false;
    end if;
    fac := Factorization(f);
    return #fac eq 5 and &and [ Degree(ff[1]) eq 1 and ff[2] eq 1 : ff in fac ];
end function;

function IsFullSplitFinite(f)
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        return false;
    end if;
    roots := Roots(f);
    return #roots eq 5 and &and [ rt[2] eq 1 : rt in roots ];
end function;

function ReduceRationalOrInfinity(r, p)
    den := Z!Denominator(r);
    if den mod p eq 0 then
        return -1;
    end if;
    F := GF(p);
    return Integers()!(F!(Z!Numerator(r)) / F!den);
end function;

function LocalBadAndSplitCodes(p)
    F := GF(p);
    PF<X> := PolynomialRing(F);
    bad_codes := { Integers() | };
    split_codes := { Integers() | };
    total := 0;
    smooth := 0;
    bad := 0;
    full_split := 0;

    for a in F do
        for b in F do
            total +:= 1;
            h := 1 + a*X + b*X^2;
            k := (1+a+b)^2;
            f := h^2 - k*X^5;
            code := (Integers()!a)*p + (Integers()!b);
            if k eq 0 or Degree(f) ne 5 or Discriminant(f) eq 0 then
                Include(~bad_codes, code);
                bad +:= 1;
                continue;
            end if;
            smooth +:= 1;
            if IsFullSplitFinite(f) then
                Include(~split_codes, code);
                full_split +:= 1;
            end if;
        end for;
    end for;

    return bad_codes, split_codes, total, smooth, bad, full_split;
end function;

function PairAllowedByLocalCode(ra, rb, p, bad_codes, split_codes, require_bad)
    if ra eq -1 or rb eq -1 then
        return true;
    end if;
    code := ra*p + rb;
    if code in bad_codes then
        return true;
    end if;
    if require_bad then
        return false;
    end if;
    return code in split_codes;
end function;

function PassesSplitFilter(f)
    for p in filter_primes do
        try
            fp := ChangeRing(f, GF(p));
        catch e
            continue;
        end try;
        if Degree(fp) ne 5 or Discriminant(fp) eq 0 then
            continue;
        end if;
        if not IsFullSplitFinite(fp) then
            return false;
        end if;
    end for;
    return true;
end function;

function IrreducibleFrobeniusCertificate(f)
    C := HyperellipticCurve(f);
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71] do
        try
            fp := ChangeRing(f, GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then
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
    return false, 0, P!0;
end function;

procedure VerifyFiveTorsion(a, b)
    f, h, k := FiveTorsionPolynomial(a, b);
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    D := J![x, Q!1];
    assert 5*D eq J!0;
    assert D ne J!0;
end procedure;

procedure RunFinite()
    print "Finite-field split-density check for the quintic-contact 5-torsion family";
    for p in finite_primes do
        F := GF(p);
        PF<X> := PolynomialRing(F);
        total := 0;
        smooth := 0;
        full_split := 0;
        five_ok := 0;
        samples := [];

        for a in F do
            for b in F do
                total +:= 1;
                h := 1 + a*X + b*X^2;
                k := (1+a+b)^2;
                if k eq 0 then
                    continue;
                end if;
                f := h^2 - k*X^5;
                if Degree(f) ne 5 or Discriminant(f) eq 0 then
                    continue;
                end if;
                smooth +:= 1;
                if not IsFullSplitFinite(f) then
                    continue;
                end if;
                full_split +:= 1;
                C := HyperellipticCurve(f);
                if (#Jacobian(C) mod 5) eq 0 then
                    five_ok +:= 1;
                end if;
                if #samples lt 8 then
                    Append(~samples, <Integers()!a, Integers()!b>);
                end if;
            end for;
        end for;
        print "p", p, "total", total, "smooth", smooth,
              "full_split", full_split, "five_ok", five_ok,
              "samples", samples;
    end for;
end procedure;

procedure RunSearch()
    params := RationalParametersOfHeight(height);
    checked := 0;
    smooth := 0;
    filter_survivors := 0;
    full_split := 0;
    torsion_tests := 0;
    hits := [];

    print "Quintic-contact 5-torsion rational split search";
    print "height", height, "parameters", #params;
    print "filter_primes", filter_primes;

    for a in params do
        for b in params do
            if 1+a+b eq 0 then
                continue;
            end if;
            checked +:= 1;
            f, h, k := FiveTorsionPolynomial(a, b);
            if Degree(f) ne 5 or Discriminant(f) eq 0 then
                continue;
            end if;
            smooth +:= 1;
            if not PassesSplitFilter(f) then
                continue;
            end if;
            filter_survivors +:= 1;
            if not IsFullSplitPolynomial(f) then
                continue;
            end if;
            full_split +:= 1;

            C := HyperellipticCurve(f);
            J := Jacobian(C);
            D := J![x, Q!1];
            torsion_tests +:= 1;
            G, phi := TorsionSubgroup(J);
            invs := Invariants(G);
            simple, pcert, Lp := IrreducibleFrobeniusCertificate(f);
            Append(~hits, <a,b,k,invs,simple,pcert,Lp,f>);
            print "HIT", "a", a, "b", b, "k", k,
                  "orderD", Order(D), "torsion", invs,
                  "simple", simple, "pcert", pcert, "f", f;
            if #hits ge max_hits then
                break a;
            end if;
        end for;
    end for;

    print "DONE height", height;
    print "checked", checked, "smooth", smooth,
          "filter_survivors", filter_survivors,
          "full_split", full_split,
          "torsion_tests", torsion_tests,
          "hits", #hits;
end procedure;

procedure RunBoundarySearch()
    params := RationalParametersOfHeight(height);
    residues := [ [ ReduceRationalOrInfinity(r, p) : p in boundary_filter_primes ]
                  : r in params ];

    bad_sets := [];
    split_sets := [];

    print "Boundary-focused quintic-contact 5-torsion search";
    print "height", height, "parameters", #params;
    print "boundary_primes", boundary_primes;
    print "extra_split_primes", boundary_extra_split_primes;
    print "Local residue summaries:";

    for p in boundary_filter_primes do
        bad_codes, split_codes, total, smooth, bad, full_split :=
            LocalBadAndSplitCodes(p);
        Append(~bad_sets, bad_codes);
        Append(~split_sets, split_codes);
        print "p", p, "total", total, "smooth", smooth,
              "bad_boundary", bad, "full_split_good", full_split;
    end for;

    checked := 0;
    boundary_survivors := 0;
    split_survivors := 0;
    smooth := 0;
    full_split := 0;
    torsion_tests := 0;
    hits := [];
    stop := false;

    for ia in [1..#params] do
        if stop then
            break;
        end if;
        a := params[ia];
        ares := residues[ia];
        for ib in [1..#params] do
            b := params[ib];
            if 1+a+b eq 0 then
                continue;
            end if;
            checked +:= 1;
            bres := residues[ib];

            locally_ok := true;
            for i in [1..#boundary_primes] do
                p := boundary_filter_primes[i];
                if not PairAllowedByLocalCode(ares[i], bres[i], p,
                                              bad_sets[i], split_sets[i],
                                              true) then
                    locally_ok := false;
                    break;
                end if;
            end for;
            if not locally_ok then
                continue;
            end if;
            boundary_survivors +:= 1;

            for i in [#boundary_primes+1..#boundary_filter_primes] do
                p := boundary_filter_primes[i];
                if not PairAllowedByLocalCode(ares[i], bres[i], p,
                                              bad_sets[i], split_sets[i],
                                              false) then
                    locally_ok := false;
                    break;
                end if;
            end for;
            if not locally_ok then
                continue;
            end if;
            split_survivors +:= 1;

            f, h, k := FiveTorsionPolynomial(a, b);
            if Degree(f) ne 5 or Discriminant(f) eq 0 then
                continue;
            end if;
            smooth +:= 1;
            if not IsFullSplitPolynomial(f) then
                continue;
            end if;
            full_split +:= 1;

            C := HyperellipticCurve(f);
            J := Jacobian(C);
            D := J![x, Q!1];
            torsion_tests +:= 1;
            G, phi := TorsionSubgroup(J);
            invs := Invariants(G);
            simple, pcert, Lp := IrreducibleFrobeniusCertificate(f);
            Append(~hits, <a,b,k,invs,simple,pcert,Lp,f>);
            print "HIT", "a", a, "b", b, "k", k,
                  "orderD", Order(D), "torsion", invs,
                  "simple", simple, "pcert", pcert, "f", f;
            if #hits ge max_hits then
                stop := true;
                break;
            end if;
        end for;
    end for;

    print "DONE boundary height", height;
    print "checked", checked,
          "boundary_survivors", boundary_survivors,
          "split_survivors", split_survivors,
          "smooth", smooth,
          "full_split", full_split,
          "torsion_tests", torsion_tests,
          "hits", #hits;
end procedure;

if mode eq "finite" then
    RunFinite();
elif mode eq "search" then
    RunSearch();
elif mode eq "boundary" then
    RunBoundarySearch();
else
    error "unknown mode";
end if;

quit;
