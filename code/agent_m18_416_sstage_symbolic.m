
//////////////////////////////////////////////////////////////////////
//  Symbolic second-stage tower for P_R halving in M_1(8,4), on the
//  Pell-parametrized C1 locus.
//
//  Setup (validated in agent_m18_416_descent_conditions.m):
//    K = -2R(R^2-1),  w = (m^2+K)/(m^2-K)  (free C1 parametrization),
//    C2 <=> G(R,m) = 2(R^2-1)*(R(2R+1) - W(R+2)) = square,  W = w^2.
//  Second stage at the A-component (quadratic At, monic):
//    u_A = X_R - theta_A,  alpha_A = X_R + g1_A/2,  N_A = At(X_R),
//    S_A <=> 2*(alpha_A + n_A) in {1, d_A} mod squares, n_A^2 = N_A.
//  On the double cover g^2 = G:  N_A = G*h_A^2  ==>  n_A = +-g*h_A.
//
//  This script computes and FACTORS over Q(R,m):
//    d_A, d_B (component discriminants),
//    alpha_A, alpha_B,
//    N_A, N_B, and the exact cofactors h_A^2 = N_A/G, h_B^2 = N_B/(C3class)
//  making the whole tower explicit.  Every factorization is printed with
//  multiplicities so the squarefree cores (the actual squareclass
//  content) are visible.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();

K2<R,m> := RationalFunctionField(Q, 2);
PX<x> := PolynomialRing(K2);

K := -2*R*(R^2-1);
w := (m^2 + K)/(m^2 - K);
W := w^2;

t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
c4 := R + 2 + 4*t;
Apol := x^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*x + R^4;
Bpol := c4*x^2 + (R^2 + 4*R + 1 + 8*t)*x + (2*R^2 + R + 4*t);

// sanity: c4 = 2(R^2-1)/(w^2-1)
assert c4 eq 2*(R^2-1)/(w^2-1);

XR := -c4*R;
// sanity: C1 on the parametrization is a square
sQ := 2*m/(m^2 - K);        // Pell: w^2 - K*s^2 = 1 with s = 2m/(m^2-K)
assert w^2 - K*sQ^2 eq 1;
c1val := -c4*R;
// -c4*R = 2R(R^2-1)/(1-w^2) * ... check squareness directly:
issq1, rt1 := IsSquare(c1val);
printf "C1 = -c4*R is a square in Q(R,m): %o\n", issq1;

// monic transforms
At := c4^2*Evaluate(Apol, x/c4);
Bt := c4*Evaluate(Bpol, x/c4);
At := PX![K2!co : co in Coefficients(At)];
Bt := PX![K2!co : co in Coefficients(Bt)];
assert LeadingCoefficient(At) eq 1 and LeadingCoefficient(Bt) eq 1;

G := 2*(R^2-1)*(R*(2*R+1) - W*(R+2));

function FactorRat(val)   // factor numerator and denominator with mults
    num := Numerator(val); den := Denominator(val);
    return Factorization(num), Factorization(den);
end function;

procedure ShowFactored(name, val)
    fn, fd := FactorRat(val);
    printf "%o:\n  num: %o\n  den: %o\n", name, fn, fd;
    // squarefree core
    core := K2!1;
    for pr in fn do if IsOdd(pr[2]) then core *:= K2!pr[1]; end if; end for;
    for pr in fd do if IsOdd(pr[2]) then core *:= K2!pr[1]; end if; end for;
    printf "  squarefree core: %o\n", core;
end procedure;

print "\n== G (the C2 conic-bundle condition) ==";
ShowFactored("G", G);

for pair in [<"A", At>, <"B", Bt>] do
    nm := pair[1]; gq := pair[2];
    printf "\n== component %o ==\n", nm;
    g1 := Coefficient(gq, 1);
    dsc := Discriminant(gq);
    alpha := XR + g1/2;
    Nval := Evaluate(gq, XR);
    ShowFactored(Sprintf("d_%o = disc", nm), dsc);
    ShowFactored(Sprintf("alpha_%o", nm), alpha);
    ShowFactored(Sprintf("N_%o", nm), Nval);
    // cofactor: N / G should be a perfect square for A (C2 <=> G sq),
    // and N_B / (C1-class * G) or similar for B; test both
    for div0 in [<"N/G", Nval/G>, <"N/(G*C1)", Nval/(G*c1val)>] do
        lbl := div0[1]; ratio := div0[2];
        issq, rt := IsSquare(ratio);
        printf "  %o is a perfect square in Q(R,m): %o\n", lbl, issq;
        if issq then
            printf "    sqrt = h_%o with:\n", nm;
            ShowFactored(Sprintf("    h_%o", nm), rt);
        end if;
    end for;
end for;

print "\n== the second-stage objects on the cover g^2 = G ==";
print "S_A <=> 2*(alpha_A + g*h_A) in {1, d_A} mod squares   (n_A = g*h_A)";
print "S_B <=> 2*(alpha_B + n_B)   in {1, d_B} mod squares";
print "with the exact identities 2(a+n)*2(a-n) = 4*d*beta^2 == d mod sq.";
print "DONE";
quit;
