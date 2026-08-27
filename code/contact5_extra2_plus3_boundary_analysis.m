//////////////////////////////////////////////////////////////////////
//  Component-wise boundary analysis for trying to add rational
//  3-torsion to the contact-5 [2,20] families.
//
//  Obstruction-prime summary from the good-reduction sieve:
//    - linear branch: p=7 and p=13 have no good residues with 3 | #J(F_p)
//    - qq branch:     p=11 has no good residues with 3 | #J(F_p)
//
//  This script classifies the bad residues at those primes, solves the
//  cubic-contact 3-torsion equations on the bad fibers, and lifts the
//  finite t=-3 contact points to p^2 to see whether first-order branches
//  move off the singular divisor.
//
//  Typical run:
//      magma code/contact5_extra2_plus3_boundary_analysis.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

function TLinearMod(z, R)
    den := z^4 + 4*z^3 + 8*z^2 + 8*z + 4;
    if not IsUnit(den) then
        return false, R!0;
    end if;
    return true, -(z^4 + 4*z + 4)/den;
end function;

function TQQMod(r, R)
    den := (r^2 - 2)^2*(r^2 - 2*r + 2);
    if not IsUnit(den) then
        return false, R!0;
    end if;
    num := r^6 - 2*r^5 + 2*r^4 - 4*r^3 + 4*r^2 - 8*r + 8;
    return true, -num/den;
end function;

function FamilyCoeffsMod(t, R)
    P<x> := PolynomialRing(R);
    inv2 := (R!2)^-1;
    inv4 := (R!4)^-1;
    b := (t^2 - 1)*inv2;
    h := 1 + t*x + b*x^2;
    f := h^2 - ((t + 1)^4)*inv4*x^5;
    return f, [ Coefficient(f,i) : i in [0..5] ];
end function;

function ComponentLabel(label, p, a)
    F := GF(p);
    z := F!a;
    if label eq "linear" then
        ok, t := TLinearMod(z, F);
    else
        ok, t := TQQMod(z, F);
    end if;
    if not ok then
        return "pole", F!0;
    end if;
    if t eq -F!1 then
        return "t=-1", t;
    end if;
    if t eq -F!3 then
        return "t=-3", t;
    end if;
    return "disc-cubic", t;
end function;

function ContactSolutionsFromCoeffs(coeffs, R)
    sols := [];
    for m in R do
        if not IsUnit(m) then
            continue;
        end if;
        inv2m := (2*m)^-1;
        for U in R do
            for V in R do
                N := (3*m^2*U + coeffs[6])*inv2m;
                Rr := (3*m^2*(U^2 + V) + coeffs[5] - N^2)*inv2m;
                S := (m^2*(U^3 + 6*U*V) + coeffs[4] - 2*N*Rr)*inv2m;
                e2 := Rr^2 + 2*N*S - coeffs[3] - 3*m^2*(U^2*V + V^2);
                e1 := 2*Rr*S - coeffs[2] - 3*m^2*U*V^2;
                e0 := S^2 - coeffs[1] - m^2*V^3;
                if e2 eq 0 and e1 eq 0 and e0 eq 0 then
                    Append(~sols, <m,U,V,N,Rr,S>);
                end if;
            end for;
        end for;
    end for;
    return sols;
end function;

function ContactStats(f)
    F := BaseRing(Parent(f));
    P<x> := Parent(f);
    coeffs := [ Coefficient(f,i) : i in [0..5] ];
    sols := ContactSolutionsFromCoeffs(coeffs, F);
    nd_q := 0;
    gcd1 := 0;
    for sol in sols do
        U := sol[2];
        V := sol[3];
        q := x^2 + U*x + V;
        if Discriminant(q) ne 0 then
            nd_q +:= 1;
        end if;
        if GCD(q, f) eq 1 then
            gcd1 +:= 1;
        end if;
    end for;
    return sols, nd_q, gcd1;
end function;

function LiftCountP2(label, p, a, sol)
    n := p^2;
    R := Integers(n);
    m0 := Integers()!sol[1];
    U0 := Integers()!sol[2];
    V0 := Integers()!sol[3];
    total := 0;
    on_tminus3 := 0;
    off_tminus3 := 0;

    for za in [0..p-1] do
        z := R!(a + p*za);
        if label eq "linear" then
            ok, t := TLinearMod(z, R);
        else
            ok, t := TQQMod(z, R);
        end if;
        if not ok then
            continue;
        end if;
        f, coeffs := FamilyCoeffsMod(t, R);
        for ma in [0..p-1] do
            for Ua in [0..p-1] do
                for Va in [0..p-1] do
                    m := R!(m0 + p*ma);
                    U := R!(U0 + p*Ua);
                    V := R!(V0 + p*Va);
                    if not IsUnit(m) then
                        continue;
                    end if;
                    inv2m := (2*m)^-1;
                    N := (3*m^2*U + coeffs[6])*inv2m;
                    Rr := (3*m^2*(U^2 + V) + coeffs[5] - N^2)*inv2m;
                    S := (m^2*(U^3 + 6*U*V) + coeffs[4] - 2*N*Rr)*inv2m;
                    e2 := Rr^2 + 2*N*S - coeffs[3] - 3*m^2*(U^2*V + V^2);
                    e1 := 2*Rr*S - coeffs[2] - 3*m^2*U*V^2;
                    e0 := S^2 - coeffs[1] - m^2*V^3;
                    if e2 eq 0 and e1 eq 0 and e0 eq 0 then
                        total +:= 1;
                        if t + 3 eq 0 then
                            on_tminus3 +:= 1;
                        else
                            off_tminus3 +:= 1;
                        end if;
                    end if;
                end for;
            end for;
        end for;
    end for;
    return total, on_tminus3, off_tminus3;
end function;

bad_cases := [ <"linear",7>, <"linear",13>, <"qq",11> ];

print "Boundary residue classification and mod-p contact cover";
for cc in bad_cases do
    label := cc[1];
    p := cc[2];
    F := GF(p);
    print "CASE", label, "p", p;
    for a in [0..p-1] do
        comp, t := ComponentLabel(label, p, a);
        is_bad := false;
        if comp eq "pole" or comp eq "t=-1" then
            is_bad := true;
        else
            f, coeffs := FamilyCoeffsMod(t, F);
            if Degree(f) ne 5 or Discriminant(f) eq 0 then
                is_bad := true;
            end if;
        end if;
        if not is_bad then
            continue;
        end if;
        if comp eq "pole" then
            print " residue", a, "component", comp, "modp_contact", "use infinity chart";
            continue;
        end if;
        f, coeffs := FamilyCoeffsMod(t, F);
        sols, nd_q, gcd1 := ContactStats(f);
        print " residue", a, "component", comp, "t", t, "degree", Degree(f),
              "disc", Discriminant(f), "contact_total", #sols,
              "nondeg_q", nd_q, "gcd1", gcd1;
    end for;
end for;

print "";
print "First-order p^2 lifts of finite t=-3 contact points";
lift_cases := [ <"linear",7,[3,5]>, <"linear",13,[5,11]>, <"qq",11,[3,5,7,8]> ];
for cc in lift_cases do
    label := cc[1];
    p := cc[2];
    residues := cc[3];
    F := GF(p);
    print "CASE", label, "p", p;
    for a in residues do
        comp, t := ComponentLabel(label, p, a);
        f, coeffs := FamilyCoeffsMod(t, F);
        sols, nd_q, gcd1 := ContactStats(f);
        total_lifts := 0;
        total_on := 0;
        total_off := 0;
        for sol in sols do
            total, on, off := LiftCountP2(label, p, a, sol);
            total_lifts +:= total;
            total_on +:= on;
            total_off +:= off;
        end for;
        print " residue", a, "component", comp, "modp_solutions", #sols,
              "p2_lifts", total_lifts, "on_t=-3", total_on,
              "off_t=-3", total_off;
    end for;
end for;

print "";
print "Infinity chart for pole components";
for p in [13] do
    F := GF(p);
    P<x> := PolynomialRing(F);
    // s=1/t, y' = 2*s^2*y gives boundary fiber y'^2 = x^4 - x^5.
    f0 := x^4 - x^5;
    sols, nd_q, gcd1 := ContactStats(f0);
    print "p", p, "pole fiber x^4-x^5", "contact_total", #sols,
          "nondeg_q", nd_q, "gcd1", gcd1;
end for;

quit;
