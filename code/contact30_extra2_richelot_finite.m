//////////////////////////////////////////////////////////////////////
//  Finite diagnostic for a Richelot route from order 30 to order 60.
//
//  In the simultaneous contact-5/contact-6 family,
//
//      f = Q2*C3,  deg(Q2)=2, deg(C3)=3.
//
//  If C3 has a rational root rho, then
//
//      f = Q2 * Q2b * (x-rho)
//
//  gives a rational maximal-isotropic Richelot kernel (the linear
//  factor is paired with infinity).  The odd order-15 subgroup survives
//  the 2-isogeny.  Hence a Richelot codomain with 2-primary exponent 4
//  has an element of order 60.
//
//  This script checks, over small finite fields, whether that precise
//  distinguished Richelot neighbor can have exponent divisible by 60.
//  It is only a local viability diagnostic, not a rational search.
//
//  Typical run:
//      magma -b primes:="7,11,13,17,19,23,29,31" \
//          code/contact30_extra2_richelot_finite.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(4*10^9);

Z := Integers();

if assigned primes then
    if Type(primes) eq MonStgElt then
        prime_list := [StringToInteger(s) : s in Split(primes, ",") | #s gt 0];
    else
        prime_list := primes;
    end if;
else
    prime_list := [7,11,13,17,19,23,29,31];
end if;
prime_list := [Z!p : p in prime_list | p notin {2,3,5}];

function TorsionExponent(invs)
    e := 1;
    for n in invs do
        e := LCM(e, Z!n);
    end for;
    return e;
end function;

function FamilyFp(a, branch, F)
    P<x> := PolynomialRing(F);
    denR := a^2-F!5;
    if denR eq 0 then
        return false, P!0, P!0, P!0;
    end if;
    t := (F!5*a^2-F!20*a+F!19)/denR;
    Y := -F!2*(F!5*a^2-F!22*a+F!25)/denR;
    u := t^3;
    if u eq 0 then
        return false, P!0, P!0, P!0;
    end if;
    s := t^5+t^4+(F!5/F!2)*t^3+(F!1/F!2)*t
       + F!branch*t*(t-F!1/F!2)*(t+1)*Y;
    C := (u^2+1)/(F!2*u);
    c := (u^2-1)/(F!2*u);
    if c eq 0 then
        return false, P!0, P!0, P!0;
    end if;
    denq := u^6+F!6*u^4*s-F!2*u^4+F!15*u^3*s-u*s^3+u^2;
    if denq eq 0 then
        return false, P!0, P!0, P!0;
    end if;
    numq := F!15*u^5+F!90*u^4+F!20*u^3*s-F!6*u^2*s^2
          + F!231*u^3+F!2*u^2*s-F!15*u*s^2+F!90*u^2
          - F!20*u*s+F!15*u-F!2*s;
    q := numq/denq;
    A := (s+q)/F!2;
    B := (F!15-s*q)/F!2;
    h6 := x^3+A*x^2+B*x+C;
    Q2 := h6-(x-1)^3;
    C3 := h6+(x-1)^3;
    f := Q2*C3;
    return true, f, Q2, C3;
end function;

function GoodGenus2Polynomial(f)
    return Degree(f) in {5,6} and Discriminant(f) ne 0;
end function;

function Bracket(A, B)
    return Derivative(A)*B-A*Derivative(B);
end function;

function CoefficientMatrixDeterminant(A, B, C)
    F := BaseRing(Parent(A));
    M := Matrix(F, 3, 3,
        [Coefficient(A,i) : i in [0..2]] cat
        [Coefficient(B,i) : i in [0..2]] cat
        [Coefficient(C,i) : i in [0..2]]);
    return Determinant(M);
end function;

function DistinguishedRichelotPolynomial(A, B, C, source_order)
    // One of +/-Delta*product is the chosen Richelot codomain; the
    // correct sign is recognized by equality of the finite Jacobian
    // orders, which every isogenous codomain must satisfy.
    Delta := CoefficientMatrixDeterminant(A, B, C);
    if Delta eq 0 then
        return false, Parent(A)!0, "zero_determinant";
    end if;
    g0 := Delta*Bracket(B,C)*Bracket(C,A)*Bracket(A,B);
    for sign in [1,-1] do
        g := sign*g0;
        if not GoodGenus2Polynomial(g) then
            continue;
        end if;
        try
            if #Jacobian(HyperellipticCurve(g)) eq source_order then
                return true, g, sign eq 1 select "plus" else "minus";
            end if;
        catch e
            continue;
        end try;
    end for;
    return false, Parent(A)!0, "order_mismatch";
end function;

print "FINITE contact30 extra-2 distinguished Richelot route";
print "primes", prime_list;

for p in prime_list do
    F := GF(p);
    P<x> := PolynomialRing(F);
    for branch in [-1,1] do
        total := p;
        good_sources := 0;
        source_parameters := 0;
        kernels := 0;
        dual_good := 0;
        dual_exp60 := 0;
        dual_2exp4 := 0;
        formula_fail := 0;
        samples := [];

        for aa in [0..p-1] do
            ok, f, Q2, C3 := FamilyFp(F!aa, branch, F);
            if not ok or not GoodGenus2Polynomial(f) or
               Degree(Q2) ne 2 or Degree(C3) ne 3 then
                continue;
            end if;
            good_sources +:= 1;
            roots := Roots(C3);
            if #roots eq 0 then
                continue;
            end if;
            source_parameters +:= 1;
            J := Jacobian(HyperellipticCurve(f));
            source_order := #J;

            for rr in roots do
                rho := rr[1];
                L := x-rho;
                Q2b := ExactQuotient(C3, L);
                if Degree(Q2b) ne 2 then
                    continue;
                end if;
                kernels +:= 1;
                dual_ok, g, sign_label :=
                    DistinguishedRichelotPolynomial(Q2, Q2b, L, source_order);
                if not dual_ok then
                    formula_fail +:= 1;
                    continue;
                end if;
                dual_good +:= 1;
                G, phi := AbelianGroup(Jacobian(HyperellipticCurve(g)));
                invs := Invariants(G);
                exponent := TorsionExponent(invs);
                if Valuation(exponent, 2) ge 2 then
                    dual_2exp4 +:= 1;
                end if;
                if exponent mod 60 eq 0 then
                    dual_exp60 +:= 1;
                    if #samples lt 8 then
                        Append(~samples,
                            <aa,Z!rho,invs,sign_label,
                             [Degree(fe[1]) : fe in Factorization(g)]>);
                    end if;
                end if;
            end for;
        end for;

        print "p", p, "branch", branch,
              "total", total,
              "good_sources", good_sources,
              "source_parameters_with_C3_root", source_parameters,
              "distinguished_kernels", kernels,
              "dual_good", dual_good,
              "dual_2exp_at_least4", dual_2exp4,
              "dual_exp60", dual_exp60,
              "formula_fail", formula_fail,
              "samples", samples;
    end for;
end for;

quit;
