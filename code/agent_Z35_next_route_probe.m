//////////////////////////////////////////////////////////////////////
//  Z/35 next-route probe: structural analysis of the interrupted
//  simultaneous contact-7 / point-contact-5 equations.
//
//  The old search used
//
//      c2 - b = d,     c2 + b = e,
//      r = (d*e + 7)/5, b = (e-d)/2, c2 = (e+d)/2,
//
//  then solved two derivative equations for a and c0, leaving two
//  residual equations N0(d,e,c1)=N1(d,e,c1)=0.
//
//  This script eliminates c1, records the degenerate factors, and checks
//  the small-prime residue structure.  It is intentionally not another
//  broad height search.
//
//  Typical run:
//      magma code/agent_Z35_next_route_probe.m
//
//  Optional:
//      magma -b print_phi:=1 code/agent_Z35_next_route_probe.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned print_phi then
    print_phi := 0;
elif Type(print_phi) eq MonStgElt then
    print_phi := StringToInteger(print_phi);
end if;

Q := Rationals();
Z := Integers();

R<D,E,C> := PolynomialRing(Q, 3);
K := FieldOfFractions(R);
d := K!D;
e := K!E;
c1 := K!C;

r := (d*e + 7)/5;
b := (e - d)/2;
c2 := (e + d)/2;

a := (2*c1*c2 + 4*c2^2*r - 4*b^2*r - 10*r^2 + 28*r - 21)/(2*b);
c0 := (a^2 + 6*a*b*r + 6*b^2*r^2 - 7*b - c1^2
       - 6*c1*c2*r - 6*c2^2*r^2 + 10*r^3 - 42*r^2
       + 63*r - 35)/(2*c2);

// The two remaining equations from contact7_contact5_point_enum.py.
eq0 := -(4*a^2*r^2 + 8*a*b*r^3 - 28*a*r + 8*a
         + 4*b^2*r^4 - 28*b*r^2 + 8*b*r
         - 4*c0^2 - 8*c0*c1*r - 8*c0*c2*r^2
         - 4*c1^2*r^2 - 8*c1*c2*r^3 - 4*c2^2*r^4
         + 4*r^5 - 28*r^4 + 84*r^3 - 140*r^2
         + 140*r - 35)/4;
eq1 := -2*a^2*r - 6*a*b*r^2 + 7*a - 4*b^2*r^3
       + 14*b*r - 2*b + 2*c0*c1 + 4*c0*c2*r
       + 2*c1^2*r + 6*c1*c2*r^2 + 4*c2^2*r^3
       - 5*r^4 + 28*r^3 - 63*r^2 + 70*r - 35;

N0 := R!Numerator(eq0);
N1 := R!Numerator(eq1);

function IntegralPoly(poly)
    L := LCM([ Denominator(c) : c in Coefficients(poly) ]);
    ZR := PolynomialRing(Integers(), 3);
    return ZR!(L*poly), L;
end function;

function Contact7F(F, aa, bb)
    P<x> := PolynomialRing(F);
    h := 1 - (F!7/F!2)*x + aa*x^2 + bb*x^3;
    f := ExactQuotient(h^2 + (x - 1)^7, x^2);
    return f, h;
end function;

function ResidueTag(dd, ee, cc)
    F := Parent(dd);
    tags := [];
    if dd eq 0 then Append(~tags, "D=0"); end if;
    if ee eq 0 then Append(~tags, "E=0"); end if;
    if dd eq ee then Append(~tags, "D=E"); end if;
    if dd + ee eq 0 then Append(~tags, "D=-E"); end if;
    if F!(dd*ee + 2) eq 0 then Append(~tags, "r=1"); end if;
    if #tags eq 0 then
        return "nondegenerate";
    end if;
    return Join(tags, ",");
end function;

function NormalizedProjectiveTriple(uu, vv, ww)
    F := Parent(uu);
    if uu eq 0 and vv eq 0 and ww eq 0 then
        return "zero";
    end if;
    arr := [uu, vv, ww];
    idx := Min([ i : i in [1..3] | arr[i] ne 0 ]);
    scale := arr[idx]^-1;
    return Sprint(<scale*arr[1], scale*arr[2], scale*arr[3]>);
end function;

S<U,V,W> := PolynomialRing(GF(3), 3);

function LowPartAt(poly, ctr)
    RF := Parent(ChangeRing(poly, GF(3)));
    phi := hom<RF -> S | S!(ctr[1] + U), S!(ctr[2] + V), S!(ctr[3] + W)>;
    g := phi(ChangeRing(poly, GF(3)));
    mons := Monomials(g);
    coeffs := Coefficients(g);
    if #mons eq 0 then
        return -1, S!0;
    end if;
    mindeg := Minimum([ Degree(m) : m in mons ]);
    low := S!0;
    for i in [1..#mons] do
        if Degree(mons[i]) eq mindeg then
            low +:= coeffs[i]*mons[i];
        end if;
    end for;
    return mindeg, low;
end function;

function EvalMod(poly, vals, m)
    return Evaluate(poly, [ Z!v : v in vals ]) mod m;
end function;

procedure PointFiniteTable(primes)
    print "POINT_CONTACT_FINITE_TABLE";
    print "columns: p residual_solutions smooth_open good_pass_5 bad_or_boundary";
    for p in primes do
        F := GF(p);
        N0p := ChangeRing(N0, F);
        N1p := ChangeRing(N1, F);
        total := 0;
        smooth := 0;
        pass5 := 0;
        bad := 0;

        for dd in F do
            for ee in F do
                for cc in F do
                    if dd eq 0 or ee eq 0 or dd eq ee or dd + ee eq 0 then
                        continue;
                    end if;
                    if F!(dd*ee + 2) eq 0 then
                        continue;
                    end if;
                    if Evaluate(N0p, [dd, ee, cc]) ne 0
                       or Evaluate(N1p, [dd, ee, cc]) ne 0 then
                        continue;
                    end if;

                    total +:= 1;
                    rr := (dd*ee + 7)/5;
                    bb := (ee - dd)/2;
                    c2f := (ee + dd)/2;
                    aa := (2*cc*c2f + 4*c2f^2*rr - 4*bb^2*rr
                           - 10*rr^2 + 28*rr - 21)/(2*bb);
                    f, h := Contact7F(F, aa, bb);
                    if Degree(f) eq 5 and Discriminant(f) ne 0
                       and Evaluate(h, F!1) ne 0 then
                        smooth +:= 1;
                        n := Z!#Jacobian(HyperellipticCurve(f));
                        if n mod 5 eq 0 then
                            pass5 +:= 1;
                        end if;
                    else
                        bad +:= 1;
                    end if;
                end for;
            end for;
        end for;
        print p, total, smooth, pass5, bad;
    end for;
end procedure;

procedure PrintF3Centers()
    F := GF(3);
    N0p := ChangeRing(N0, F);
    N1p := ChangeRing(N1, F);
    centers := [];
    cats := AssociativeArray();

    for dd in F do
        for ee in F do
            for cc in F do
                if Evaluate(N0p, [dd, ee, cc]) eq 0
                   and Evaluate(N1p, [dd, ee, cc]) eq 0 then
                    tag := ResidueTag(dd, ee, cc);
                    Append(~centers, <Z!dd, Z!ee, Z!cc, tag>);
                    if not IsDefined(cats, tag) then
                        cats[tag] := 0;
                    end if;
                    cats[tag] +:= 1;
                end if;
            end for;
        end for;
    end for;

    print "F3_RESIDUAL_CENTERS", #centers;
    for item in centers do
        print item;
    end for;
    print "F3_CENTER_TAG_COUNTS";
    for key in Sort([ k : k in Keys(cats) ]) do
        print key, cats[key];
    end for;
end procedure;

procedure PrintCenterGeometry()
    centers := [ <0,1,1>, <2,0,2>, <2,1,0>, <1,1,1>, <2,2,2> ];
    F := GF(3);

    print "F3_CENTER_LOW_PARTS_AND_TANGENTS";
    for ctr in centers do
        deg0, low0 := LowPartAt(N0, ctr);
        deg1, low1 := LowPartAt(N1, ctr);
        dirs := [];
        for uu in F do
            for vv in F do
                for ww in F do
                    if uu eq 0 and vv eq 0 and ww eq 0 then
                        continue;
                    end if;
                    if Evaluate(low0, [uu, vv, ww]) eq 0
                       and Evaluate(low1, [uu, vv, ww]) eq 0 then
                        Append(~dirs, NormalizedProjectiveTriple(uu, vv, ww));
                    end if;
                end for;
            end for;
        end for;
        dirs := Sort(SetToSequence(SequenceToSet(dirs)));
        print "CENTER", ctr;
        print "  N0_low_degree", deg0, "N0_low_factor", Factorization(low0);
        print "  N1_low_degree", deg1, "N1_low_factor", Factorization(low1);
        print "  projective_tangent_directions", #dirs, dirs;

        dd := F!ctr[1];
        ee := F!ctr[2];
        cc := F!ctr[3];
        if dd ne ee and dd + ee ne 0 then
            rr := (dd*ee + 7)/5;
            bb := (ee - dd)/2;
            c2f := (ee + dd)/2;
            aa := (2*cc*c2f + 4*c2f^2*rr - 4*bb^2*rr
                   - 10*rr^2 + 28*rr - 21)/(2*bb);
            f, h := Contact7F(F, aa, bb);
            print "  contact7_reduction", "a", Z!aa, "b", Z!bb,
                  "h1", Z!Evaluate(h, F!1), "disc", Z!Discriminant(f),
                  "factor", Factorization(f);
        else
            print "  contact7_reduction formula has b=0 or c2=0 pole";
        end if;
    end for;
end procedure;

procedure PrintLiftCounts()
    centers := [ <0,1,1>, <2,0,2>, <2,1,0>, <1,1,1>, <2,2,2> ];
    N0Z, L0 := IntegralPoly(N0);
    N1Z, L1 := IntegralPoly(N1);
    print "INTEGRAL_CLEARING_DENOMINATORS", L0, L1;
    print "THREE_ADIC_LIFT_COUNTS";
    print "Each center is lifted as center + 3*(u,v,w) modulo 3^k.";

    for m in [9, 27] do
        print "MODULUS", m;
        step_bound := m div 3 - 1;
        for ctr in centers do
            count := 0;
            dirs := AssociativeArray();
            for uu in [0..step_bound] do
                for vv in [0..step_bound] do
                    for ww in [0..step_bound] do
                        vals := [Z!ctr[1] + 3*uu, Z!ctr[2] + 3*vv, Z!ctr[3] + 3*ww];
                        if EvalMod(N0Z, vals, m) eq 0
                           and EvalMod(N1Z, vals, m) eq 0 then
                            count +:= 1;
                            if m eq 9 then
                                key := NormalizedProjectiveTriple(GF(3)!uu, GF(3)!vv, GF(3)!ww);
                                if not IsDefined(dirs, key) then
                                    dirs[key] := 0;
                                end if;
                                dirs[key] +:= 1;
                            end if;
                        end if;
                    end for;
                end for;
            end for;
            print "CENTER", ctr, "lifts", count;
            if m eq 9 then
                for key in Sort([ k : k in Keys(dirs) ]) do
                    print " ", key, dirs[key];
                end for;
            end if;
        end for;
    end for;
end procedure;

print "Z35 next-route probe: simultaneous contact7/contact5 point structure";
print "N0 degree/terms", Degree(N0), #Terms(N0);
print "N1 degree/terms", Degree(N1), #Terms(N1);

Res := Resultant(N0, N1, C);
fac := Factorization(Res);
print "RESULTANT_IN_D_E";
print "degree", Degree(Res), "terms", #Terms(Res), "number_of_factors", #fac;
print "factor_summary", [ <Degree(item[1]), item[2], #Terms(item[1])> : item in fac ];
print "linear_degenerate_factors";
for item in fac do
    if Degree(item[1]) le 1 then
        print item;
    end if;
end for;
Phi := fac[#fac][1];
print "nondegenerate_factor_degree_terms", Degree(Phi), #Terms(Phi);
if print_phi ne 0 then
    print "PHI38";
    print Phi;
end if;

PointFiniteTable([3,7,11,13,17,19,23,29,31]);
PrintF3Centers();
PrintCenterGeometry();
PrintLiftCounts();

quit;
