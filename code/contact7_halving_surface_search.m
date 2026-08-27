//////////////////////////////////////////////////////////////////////
//  First-halving surface inside the contact-7 rational-root family.
//
//  Start from
//
//      h = 1 - (7/2)*x + a*x^2 + b*x^3,
//      f = (h^2 + (x - 1)^7)/x^2.
//
//  Force a rational Weierstrass point by setting r = 1 - s^2 and
//  h(r) = eps*s^7, eps = +/-1.  This makes f(r)=0 and leaves b free.
//  Translate X = x-r and write
//
//      f(X+r) / X = X^4 + c3*X^3 + c2*X^2 + c1*X + c0.
//
//  Zarhin's formula says the class (r,0)-infinity is divisible by 2
//  over Q if this quartic is
//
//      prod_i (X + t_i^2),
//
//  where the elementary symmetric functions of t_i are rational.  If
//  these symmetric functions are u,v,w,z, then
//
//      c3 = u^2 - 2*v,
//      c2 = v^2 - 2*u*w + 2*z,
//      c1 = w^2 - 2*v*z,
//      c0 = z^2.
//
//  For fixed s,z, the last equation determines b linearly.  For fixed
//  s,z,u, the first two equations determine v,w, so only the c1 equation
//  remains.  The corresponding order-4 half is
//
//      Q = X^2 - v*X + z,
//      alpha = (u*v - w)*X - u*z.
//
//  Typical runs:
//      magma -b mode:="check_known" code/contact7_halving_surface_search.m
//      magma -b mode:="search" s_height:=8 u_height:=8 z_height:=20 \
//          code/contact7_halving_surface_search.m
//////////////////////////////////////////////////////////////////////

if not assigned mode then
    mode := "check_known";
end if;
if not assigned s_height then
    s_height := 8;
elif Type(s_height) eq MonStgElt then
    s_height := StringToInteger(s_height);
end if;
if not assigned u_height then
    u_height := 8;
elif Type(u_height) eq MonStgElt then
    u_height := StringToInteger(u_height);
end if;
if not assigned z_height then
    z_height := 20;
elif Type(z_height) eq MonStgElt then
    z_height := StringToInteger(z_height);
end if;
if not assigned max_exact then
    max_exact := 200;
elif Type(max_exact) eq MonStgElt then
    max_exact := StringToInteger(max_exact);
end if;
if not assigned progress_interval then
    progress_interval := 100000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;

Qq := Rationals();
P<x> := PolynomialRing(Qq);

target_primes := [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97];

function RationalParametersOfHeight(B)
    vals := [];
    seen := {};
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num, den) ne 1 then
                continue;
            end if;
            q := Qq!num/Qq!den;
            key := Sprint(q);
            if key notin seen then
                Include(~seen, key);
                Append(~vals, q);
            end if;
        end for;
    end for;
    return vals;
end function;

function IntegralModel(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return P!(L^2*f), L;
end function;

function TorsionOrder(invs)
    if #invs eq 0 then
        return 1;
    end if;
    return &*invs;
end function;

function SimpleCertificate(f)
    C := HyperellipticCurve(f);
    for p in target_primes do
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

function Contact7RootData(s, b, eps)
    r := 1 - s^2;
    if r eq 0 then
        return false, P!0, P!0, P!0, Qq!0, Qq!0;
    end if;

    a := (eps*s^7 - 1 + (Qq!7/2)*r - b*r^3)/r^2;
    h := 1 - (Qq!7/2)*x + a*x^2 + b*x^3;
    f := ExactQuotient(h^2 + (x - 1)^7, x^2);
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        return false, f, h, P!0, a, r;
    end if;
    if Evaluate(f, r) ne 0 then
        return false, f, h, P!0, a, r;
    end if;
    g := ExactQuotient(Evaluate(f, x + r), x);
    return true, f, h, g, a, r;
end function;

function BFromSZEps(s, z, eps)
    ok0, f0, h0, g0, a0, r0 := Contact7RootData(s, Qq!0, eps);
    ok1, f1, h1, g1, a1, r1 := Contact7RootData(s, Qq!1, eps);
    if not ok0 or not ok1 then
        return false, Qq!0;
    end if;
    c00 := Coefficient(g0, 0);
    c01 := Coefficient(g1, 0);
    lambda := c01 - c00;
    if lambda eq 0 then
        return false, Qq!0;
    end if;
    return true, (z^2 - c00)/lambda;
end function;

function HalvingSurfacePoint(s, u, z, eps)
    if s eq 0 or s^2 eq 1 or u eq 0 then
        return false, Qq!0, Qq!0, Qq!0, Qq!0, P!0, P!0, P!0, Qq!0, Qq!0;
    end if;

    okb, b := BFromSZEps(s, z, eps);
    if not okb then
        return false, b, Qq!0, Qq!0, Qq!0, P!0, P!0, P!0, Qq!0, Qq!0;
    end if;

    ok, f, h, g, a, r := Contact7RootData(s, b, eps);
    if not ok then
        return false, b, Qq!0, Qq!0, Qq!0, f, h, g, a, r;
    end if;

    c3 := Coefficient(g, 3);
    c2 := Coefficient(g, 2);
    c1 := Coefficient(g, 1);
    c0 := Coefficient(g, 0);
    if c0 ne z^2 then
        return false, b, Qq!0, Qq!0, Qq!0, f, h, g, a, r;
    end if;

    v := (u^2 - c3)/2;
    w := (v^2 + 2*z - c2)/(2*u);
    if c1 ne w^2 - 2*v*z then
        return false, b, v, w, Qq!0, f, h, g, a, r;
    end if;

    return true, b, v, w, c3, f, h, g, a, r;
end function;

function PassesTargetReduction(f, target)
    used := [];
    for p in target_primes do
        try
            fp := ChangeRing(f, GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then
                continue;
            end if;
            C := HyperellipticCurve(fp);
            n := Integers()!#Jacobian(C);
            Append(~used, <p,n>);
            required := target;
            while required mod p eq 0 do
                required div:= p;
            end while;
            if required gt 1 and n mod required ne 0 then
                return false, p, n, used;
            end if;
        catch e
            continue;
        end try;
    end for;
    return true, 0, 0, used;
end function;

function ExactHalvingCheck(s, u, z, eps)
    ok, b, v, w, c3, f, h, g, a, r := HalvingSurfacePoint(s, u, z, eps);
    if not ok then
        return false, [], false, false, 0, P!0, P!0, P!0, Qq!0, Qq!0, Qq!0, Qq!0;
    end if;

    C := HyperellipticCurve(f);
    J := Jacobian(C);
    D2 := J![x - r, 0];
    D7 := J![x - 1, Evaluate(h, Qq!1)];

    X := x - r;
    Qpoly := X^2 - v*X + z;
    alpha := (u*v - w)*X - u*z;
    H4 := J![Qpoly, alpha];

    ok_exact := (Order(D2) eq 2) and (Order(D7) eq 7) and
                (Order(H4) eq 4) and (2*H4 eq D2);
    if not ok_exact then
        return false, [], false, false, 0, f, h, g, a, b, r, v;
    end if;

    pass56, pbad, nbad, used := PassesTargetReduction(f, 56);
    divisible8 := false;
    ord_half := 0;
    if pass56 then
        divisible8, H8 := IsDivisibleBy(H4, 2);
        if divisible8 then
            ord_half := Order(H8);
        end if;
    end if;

    fI, L := IntegralModel(f);
    CI := HyperellipticCurve(fI);
    JI := Jacobian(CI);
    G, phi := TorsionSubgroup(JI);
    invs := Invariants(G);
    return true, invs, pass56, divisible8, ord_half, fI, h, g, a, b, r, v;
end function;

procedure CheckKnown()
    s := Qq!1/2;
    u := Qq!4;
    z := -Qq!3/16;
    eps := -Qq!1;
    ok, invs, pass56, div8, ord8, fI, h, g, a, b, r, v := ExactHalvingCheck(s,u,z,eps);
    print "CHECK_KNOWN";
    print "s", s, "u", u, "z", z, "eps", eps;
    print "ok", ok, "a", a, "b", b, "r", r, "v", v;
    print "torsion", invs, "order", TorsionOrder(invs),
          "pass56", pass56, "divisible8", div8, "ord_half", ord8;
    print "fI", fI;
    simple, pcert, Lp := SimpleCertificate(fI);
    print "simple", simple, "pcert", pcert, "Lp", Lp;
end procedure;

procedure SearchSurface()
    s_vals := RationalParametersOfHeight(s_height);
    u_vals := RationalParametersOfHeight(u_height);
    z_vals := RationalParametersOfHeight(z_height);
    checked := 0;
    surface := 0;
    exact_tests := 0;
    hits28 := 0;
    hits56 := 0;
    first_kill := AssociativeArray();

    print "CONTACT7 first-halving surface search";
    print "s_height", s_height, "#s", #s_vals,
          "u_height", u_height, "#u", #u_vals,
          "z_height", z_height, "#z", #z_vals,
          "max_exact", max_exact;

    for s in s_vals do
        for u in u_vals do
            for z in z_vals do
                for eps in [Qq!-1, Qq!1] do
                    checked +:= 1;
                    if progress_interval gt 0 and checked mod progress_interval eq 0 then
                        print "progress", "checked", checked, "surface", surface,
                              "exact_tests", exact_tests, "hits28", hits28,
                              "hits56", hits56;
                    end if;

                    ok, b, v, w, c3, f, h, g, a, r := HalvingSurfacePoint(s,u,z,eps);
                    if not ok then
                        continue;
                    end if;
                    surface +:= 1;
                    if exact_tests ge max_exact then
                        continue;
                    end if;

                    pass56, pbad, nbad, used := PassesTargetReduction(f, 56);
                    if not pass56 then
                        if IsDefined(first_kill, pbad) then
                            first_kill[pbad] +:= 1;
                        else
                            first_kill[pbad] := 1;
                        end if;
                    end if;

                    exact_ok, invs, pass56_exact, div8, ord8, fI, h0, g0, a0, b0, r0, v0 :=
                        ExactHalvingCheck(s,u,z,eps);
                    exact_tests +:= 1;
                    if not exact_ok then
                        print "EXACT_FAIL", "s", s, "u", u, "z", z, "eps", eps,
                              "a", a, "b", b, "r", r, "v", v, "w", w;
                        continue;
                    end if;

                    if TorsionOrder(invs) mod 28 eq 0 then
                        hits28 +:= 1;
                    end if;
                    if TorsionOrder(invs) mod 56 eq 0 or div8 then
                        hits56 +:= 1;
                    end if;

                    print "SURFACE_HIT", "s", s, "u", u, "z", z, "eps", eps,
                          "a", a, "b", b, "r", r, "v", v, "w", w,
                          "torsion", invs, "pass56", pass56_exact,
                          "divisible8", div8, "ord_half", ord8;
                    print "  f =", fI;
                end for;
            end for;
        end for;
    end for;

    print "DONE contact7 first-halving surface search";
    print "checked", checked, "surface", surface, "exact_tests", exact_tests,
          "hits28", hits28, "hits56", hits56;
    print "FIRST_KILL_56";
    for p in Sort([ k : k in Keys(first_kill) ]) do
        print " ", p, first_kill[p];
    end for;
end procedure;

if mode eq "check_known" then
    CheckKnown();
elif mode eq "search" then
    SearchSurface();
else
    print "unknown mode", mode;
end if;

quit;
