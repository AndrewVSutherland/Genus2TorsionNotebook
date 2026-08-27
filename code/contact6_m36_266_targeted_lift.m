//////////////////////////////////////////////////////////////////////
//  Bounded component/lifting diagnostic for two [2,6,6] fibers.
//
//  The first fiber (eps,r,b)=(1,2,3) is the p=19 residue previously
//  singled out by an irreducible Frobenius quartic.  The second fiber
//  (1,4,5) passes the explicit non-automorphism test modulo 19.
//
//  Work with M=L^2.  For each fixed fiber the three cubic-contact
//  equations are zero-dimensional in (M,U,v).  The script:
//    * saturates and computes a lexicographic elimination basis over Q;
//    * enumerates smooth, nonboundary roots at a short list of primes;
//    * Hensel-lifts every such root to p^lift_precision;
//    * tries every requested two-prime CRT pairing and validates each
//      rational reconstruction in the exact equations.
//
//  Typical bounded run:
//      magma -b lift_precision:=4 \
//          code/contact6_m36_266_targeted_lift.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned lift_precision then
    lift_precision := 4;
elif Type(lift_precision) eq MonStgElt then
    lift_precision := StringToInteger(lift_precision);
end if;
if not assigned do_primary then
    do_primary := false;
elif Type(do_primary) eq MonStgElt then
    do_primary := do_primary in {"true", "True", "1", "yes"};
end if;
if not assigned output_file then
    output_file := "data/contact6_m36_266_targeted_lift.txt";
end if;

Q := Rationals();
Z := Integers();
RQ<M,U,v> := PolynomialRing(Q, 3, "lex");

function PrimitiveIntegralPolynomial(g)
    den := LCM([Denominator(c) : c in Coefficients(g)]);
    h := RQ!(den*g);
    nums := [Z!c : c in Coefficients(h)];
    nonzero := [Abs(n) : n in nums | n ne 0];
    if #nonzero gt 0 then
        cont := GCD(nonzero);
        if cont gt 1 then
            h := RQ!(h/cont);
        end if;
    end if;
    return h;
end function;

function FiberPolynomials(eps, r, b)
    a := (eps*(r-1)^3 - 1 - r^3 - b*r^2)/r;
    c1 := 2*a + 6;
    c2 := a^2 + 2*b - 15;
    c3 := 2*a*b + 22;
    c4 := 2*a + b^2 - 15;
    c5 := 2*b + 6;
    B := c5*M + 3*U;
    Delta := 4*c4*M + 12*(U^2 + v^2) - B^2;
    F3 := B*Delta + 16*v^3 - 8*c3*M - 8*U^3 - 48*U*v^2;
    F2 := Delta^2 + 64*B*v^3 - 64*c2*M
          - 192*(U^2*v^2 + v^4);
    F1 := Delta*v^3 - 4*c1*M - 12*U*v^4;
    return [PrimitiveIntegralPolynomial(g) : g in [F1,F2,F3]], a;
end function;

function AutoA(b,r)
    return [
        2*b^3*r^4 + 30*b^2*r^4 - 36*b^2*r^3 + 8*b^2*r^2
          + 126*b*r^4 - 296*b*r^3 + 216*b*r^2 - 48*b*r
          + 162*r^4 - 612*r^3 + 816*r^2 - 464*r + 96,
        -b^4*r^4 - 12*b^3*r^4 + 12*b^3*r^3 - 2*b^3*r^2
          - 54*b^2*r^4 + 108*b^2*r^3 - 66*b^2*r^2 + 12*b^2*r
          - 108*b*r^4 + 324*b*r^3 - 342*b*r^2 + 152*b*r - 24*b
          - 81*r^4 + 324*r^3 - 470*r^2 + 300*r - 72,
        b^4*r^3 + 6*b^3*r^3 - 6*b^3*r^2 - 14*b^2*r^2
          + 12*b^2*r - 54*b*r^3 + 102*b*r^2 - 48*b*r
          - 81*r^3 + 270*r^2 - 284*r + 96
    ];
end function;

function AutoB(b,r)
    return [
        12*b^3*r^4 - 12*b^3*r^3 + 4*b^3*r^2 + 108*b^2*r^4
          - 188*b^2*r^3 + 120*b^2*r^2 - 36*b^2*r + 4*b^2
          + 324*b*r^4 - 804*b*r^3 + 732*b*r^2 - 296*b*r + 48*b
          + 324*r^4 - 1044*r^3 + 1224*r^2 - 612*r + 108,
        -6*b^3*r^4 + 6*b^3*r^3 - 2*b^3*r^2 - 54*b^2*r^4
          + 94*b^2*r^3 - 66*b^2*r^2 + 24*b^2*r - 4*b^2
          - 162*b*r^4 + 378*b*r^3 - 342*b*r^2 + 144*b*r - 24*b
          - 162*r^4 + 450*r^3 - 470*r^2 + 216*r - 36,
        6*b^2*r^2 - 6*b^2*r + 2*b^2 + 24*b*r^3 - 24*b*r^2
          + 4*b*r + 72*r^3 - 142*r^2 + 90*r - 18
    ];
end function;

function ExactAutoLabel(b,r)
    inA := &and[e eq 0 : e in AutoA(b,r)];
    inB := &and[e eq 0 : e in AutoB(b,r)];
    return inA select "A" else (inB select "B" else "none");
end function;

function ModValue(g, z, F)
    return F!Evaluate(g, <Z!z[1],Z!z[2],Z!z[3]>);
end function;

function GoodAuxiliaryRoot(eps, r, b, z, p)
    F := GF(p);
    if (Denominator(r) mod p) eq 0 or (Denominator(b) mod p) eq 0 then
        return false;
    end if;
    rr := F!r;
    bb := F!b;
    ee := F!eps;
    mm := F!z[1];
    uu := F!z[2];
    vv := F!z[3];
    if mm eq 0 or not IsSquare(mm) or vv eq 0 or uu^2-4*vv^2 eq 0 then
        return false;
    end if;
    P<X> := PolynomialRing(F);
    aa := (ee*(rr-1)^3 - 1 - rr^3 - bb*rr^2)/rr;
    h := 1 + aa*X + bb*X^2 + X^3;
    f := h^2 - (X-1)^6;
    q := X^2 + uu*X + vv^2;
    return Degree(f) eq 5 and Discriminant(f) ne 0
           and Evaluate(h,F!1) ne 0 and Evaluate(f,rr) eq 0
           and Discriminant(q) ne 0 and Degree(GCD(q,f)) eq 0;
end function;

function SmoothRoots(polys, eps, r, b, p)
    F := GF(p);
    derivs := [[Derivative(polys[i],j) : j in [1..3]] : i in [1..3]];
    roots := [];
    for mm in [1..p-1] do
        if not IsSquare(F!mm) then
            continue;
        end if;
        for uu in [0..p-1] do
            for vv in [1..p-1] do
                z := [mm,uu,vv];
                if &or[ModValue(g,z,F) ne 0 : g in polys] then
                    continue;
                end if;
                if not GoodAuxiliaryRoot(eps,r,b,z,p) then
                    continue;
                end if;
                J := Matrix(F,3,3,
                    [ModValue(derivs[i][j],z,F) : i,j in [1..3]]);
                if Determinant(J) ne 0 then
                    Append(~roots,z);
                end if;
            end for;
        end for;
    end for;
    return roots;
end function;

function HenselLift(polys, z0, p, precision)
    z := [Z!c : c in z0];
    derivs := [[Derivative(polys[i],j) : j in [1..3]] : i in [1..3]];
    F := GF(p);
    modulus := p;
    for e in [2..precision] do
        vals := [Z!Evaluate(g,<z[1],z[2],z[3]>) : g in polys];
        if &or[val mod modulus ne 0 : val in vals] then
            return false, z;
        end if;
        J := Matrix(F,3,3,
            [F!Evaluate(derivs[i][j],<z[1],z[2],z[3]>)
             : i,j in [1..3]]);
        if Determinant(J) eq 0 then
            return false, z;
        end if;
        rhs := Matrix(F,3,1,[F!(-(val div modulus)) : val in vals]);
        delta := J^-1*rhs;
        newmodulus := modulus*p;
        z := [(z[j] + modulus*(Z!delta[j,1])) mod newmodulus
              : j in [1..3]];
        modulus := newmodulus;
    end for;
    ok := &and[(Z!Evaluate(g,<z[1],z[2],z[3]>)) mod modulus eq 0
               : g in polys];
    return ok, z;
end function;

function BalancedRationalReconstruction(a, modulus)
    aa := a mod modulus;
    if aa eq 0 then
        return true, Q!0;
    end if;
    bound := Isqrt(modulus div 2);
    r0 := modulus;
    r1 := aa;
    t0 := Z!0;
    t1 := Z!1;
    while Abs(r1) gt bound do
        quot := r0 div r1;
        r0, r1 := r1, r0-quot*r1;
        t0, t1 := t1, t0-quot*t1;
    end while;
    if t1 eq 0 or Abs(t1) gt bound or Abs(r1) gt bound then
        return false, Q!0;
    end if;
    num := r1;
    den := t1;
    if den lt 0 then
        num := -num;
        den := -den;
    end if;
    if den eq 0 or GCD(num,den) ne 1
       or (num-aa*den) mod modulus ne 0 then
        return false, Q!0;
    end if;
    return true, Q!num/den;
end function;

procedure ExactFiberAnalysis(label, polys)
    print "EXACT_FIBER", label;
    boundary := M*v*(U^2-4*v^2);
    I := ideal<RQ | polys>;
    Isat := Saturation(I, ideal<RQ | boundary>);
    dim, degs := Dimension(Isat);
    print " saturated_dimension", dim, "component_degrees", degs;
    G := GroebnerBasis(Isat);
    print " groebner_length", #G,
          "degree_term_summary",
          [<TotalDegree(g),#Terms(g)> : g in G];
    univariate_v := [g : g in G |
                     Degree(g,1) eq 0 and Degree(g,2) eq 0];
    print " univariate_v_count", #univariate_v;
    for g in univariate_v do
        print "  univariate_v_degree", Degree(g,3),
              "factor_degrees",
              [<Degree(fe[1],3),fe[2]> : fe in Factorization(g)];
        print "  univariate_v", g;
    end for;
    if do_primary then
        comps := PrimaryDecomposition(Isat);
        print " primary_components", #comps;
        for i in [1..#comps] do
            cdim, cdegs := Dimension(comps[i]);
            print "  component", i, "dimension", cdim,
                  "degrees", cdegs, "prime", IsPrime(comps[i]),
                  "basis_summary",
                  [<TotalDegree(g),#Terms(g)> : g in Basis(comps[i])];
        end for;
    end if;
end procedure;

procedure TryCRT(label, polys, p1, lifts1, p2, lifts2, precision)
    m1 := p1^precision;
    m2 := p2^precision;
    modulus := m1*m2;
    pairs := 0;
    reconstructed := 0;
    exact := 0;
    seen := { };
    for z1 in lifts1 do
        for z2 in lifts2 do
            pairs +:= 1;
            vals := [];
            ok := true;
            for j in [1..3] do
                residue := CRT([z1[j],z2[j]],[m1,m2]);
                flag, value := BalancedRationalReconstruction(residue,modulus);
                if not flag then
                    ok := false;
                    break;
                end if;
                Append(~vals,value);
            end for;
            if not ok then
                continue;
            end if;
            reconstructed +:= 1;
            key := Sprint(vals);
            if key in seen then
                continue;
            end if;
            Include(~seen,key);
            if &and[Evaluate(g,<vals[1],vals[2],vals[3]>) eq 0
                    : g in polys] then
                exact +:= 1;
                squareM, sqrtM := IsSquare(vals[1]);
                print " EXACT_CRT_POINT", label, vals,
                      "M_square", squareM,
                      "L", squareM select sqrtM else Q!0;
            end if;
        end for;
    end for;
    print "CRT_SUMMARY", label, "primes", p1, p2,
          "precision", precision, "modulus", modulus,
          "balanced_height_bound", Isqrt(modulus div 2),
          "pairs", pairs, "reconstructed", reconstructed,
          "exact_points", exact;
end procedure;

SetOutputFile(output_file);
print "Bounded [2,6,6] component/Hensel/CRT diagnostic";
print "lift_precision", lift_precision, "do_primary", do_primary;
print "go_stop", "all smooth roots; fixed primes <=31; no height box; no exact torsion unless an exact CRT point appears";

fibers := [
    <"r2_b3",Q!1,Q!2,Q!3,[13,19]>,
    <"r4_b5",Q!1,Q!4,Q!5,[11,13,17,19,23,29,31]>
];

for fib in fibers do
    label := fib[1];
    eps := fib[2];
    r := fib[3];
    b := fib[4];
    primes := fib[5];
    polys, a := FiberPolynomials(eps,r,b);
    print "\nFIBER", label, "eps", eps, "r", r, "a", a, "b", b;
    print " exact_auto_label", ExactAutoLabel(b,r),
          "AutoA", AutoA(b,r), "AutoB", AutoB(b,r);
    print " polynomial_summaries",
          [<TotalDegree(g),#Terms(g)> : g in polys];
    ExactFiberAnalysis(label,polys);

    records := [];
    for p in primes do
        if (Denominator(a) mod p) eq 0 then
            continue;
        end if;
        roots := SmoothRoots(polys,eps,r,b,p);
        lifts := [];
        for z in roots do
            ok, zp := HenselLift(polys,z,p,lift_precision);
            if ok then
                Append(~lifts,zp);
            end if;
        end for;
        autoAmod := &and[(Numerator(e) mod p) eq 0 : e in AutoA(b,r)];
        autoBmod := &and[(Numerator(e) mod p) eq 0 : e in AutoB(b,r)];
        print " PRIME", p, "smooth_good_roots", #roots,
              "lifted", #lifts,
              "auto_mod_p", autoAmod or autoBmod,
              "root_samples", roots[1..Min(#roots,4)];
        Append(~records,<p,lifts>);
    end for;

    if label eq "r2_b3" then
        rec13 := [rec : rec in records | rec[1] eq 13][1];
        rec19 := [rec : rec in records | rec[1] eq 19][1];
        TryCRT(label,polys,13,rec13[2],19,rec19[2],lift_precision);
    else
        rec19 := [rec : rec in records | rec[1] eq 19][1];
        others := [rec : rec in records | rec[1] ne 19 and #rec[2] gt 0];
        if #others gt 0 then
            TryCRT(label,polys,others[1][1],others[1][2],
                   19,rec19[2],lift_precision);
        else
            print "CRT_SUMMARY", label,
                  "no second useful prime in bounded scan";
        end if;
    end if;
end for;

print "\nDONE";
UnsetOutputFile();
print "wrote", output_file;
quit;
