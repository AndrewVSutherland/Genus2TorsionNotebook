//////////////////////////////////////////////////////////////////////
//  Generic s=m^2 condition for [4,12] on the full M(2,12) surface.
//
//  The earlier equations use
//      g(X) - (X-beta)(mX+n)^2 = g4 (X^2 + A X + B)^2.
//
//  Put s=m^2, t=mn, N=n^2.  After solving the top coefficients for
//  A and B, the remaining two equations are linear in t and N, and
//  the only nonlinear relation left is t^2=sN.  This script derives
//  the resulting single generic polynomial F_s(r,z,u,s)=0.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned mode then
    mode := "summary";
end if;
if not assigned output_file then
    output_file := "data/m12_full_surface_z12x4_s_condition.txt";
end if;

Q := Rationals();
R<r,z,u,s> := PolynomialRing(Q, 4, "grevlex");
K := FieldOfFractions(R);
Px<x> := PolynomialRing(K);
PX<X> := PolynomialRing(K);
Rbase<rb,zb,sb> := PolynomialRing(Q, 3, "grevlex");
Kbase := FieldOfFractions(Rbase);
PU<U> := PolynomialRing(Kbase);

function PrimitiveNumerator(f)
    num := R!Numerator(f);
    coeffs := Coefficients(num);
    if #coeffs eq 0 then
        return num;
    end if;
    dens := [Denominator(c) : c in coeffs | c ne 0];
    if #dens gt 0 then
        L := LCM(dens);
        num := R!(L*num);
    end if;
    coeffs := Coefficients(num);
    nums := [Integers()!c : c in coeffs | c ne 0];
    if #nums eq 0 then
        return num;
    end if;
    content := GCD([Abs(c) : c in nums]);
    if content gt 1 then
        num := R!(num/content);
    end if;
    return num;
end function;

function ReduceInU(f, q)
    to_u := hom<R -> PU | PU!(Kbase!rb), PU!(Kbase!zb), U, PU!(Kbase!sb)>;
    _, rem_u := Quotrem(to_u(f), to_u(q));
    dens := [Denominator(Coefficient(rem_u, i)) : i in [0..Degree(rem_u)] |
             Coefficient(rem_u, i) ne 0];
    if #dens eq 0 then
        return R!0;
    end if;
    den := LCM(dens);
    back := hom<Rbase -> R | r,z,s>;
    out := R!0;
    for i in [0..Degree(rem_u)] do
        ci := Coefficient(rem_u, i);
        if ci ne 0 then
            out +:= back(Rbase!(den*ci))*u^i;
        end if;
    end for;
    return PrimitiveNumerator(out);
end function;

procedure PrintPolyInfo(label, f)
    print label,
          "total_degree", TotalDegree(f),
          "degree_r", Degree(f, 1),
          "degree_z", Degree(f, 2),
          "degree_u", Degree(f, 3),
          "degree_s", Degree(f, 4),
          "terms", #Terms(f);
end procedure;

a := (1-z^2)/(4*(r+1));
T := a*x^2 - x + r;
h := (x-r)*(T+1);
W := h^2 + 4*a*x^2*T*(T+1);
Q4 := ExactQuotient(W, T+1);
w := 2*(r+1)/(1+z);
assert Evaluate(T+1, w) eq 0;

f5 := PX!0;
for i in [0..Degree(W)] do
    for j in [0..i] do
        f5 +:= Coefficient(W, i)*Binomial(i,j)*w^(i-j)*X^(6-j);
    end for;
end for;

beta := 1/(u-w);
g, rem := Quotrem(f5, X-beta);

g4 := Coefficient(g, 4);
g3 := Coefficient(g, 3);
g2 := Coefficient(g, 2);
g1 := Coefficient(g, 1);
g0 := Coefficient(g, 0);

A := (g3 - s)/(2*g4);
D := g2 + beta*s - g4*A^2;
P := 2*(beta + A);
Q0 := g1 - A*D;
L := 8*(g4*beta - s)*(beta + A) + 4*D;
C := 4*(g4*beta - s)*Q0 + 4*g4*g0 - D^2;

// Generic elimination relation after t=-C/L and N=P*t+Q0.
F := C^2 + s*P*C*L - s*Q0*L^2;

qeq := PrimitiveNumerator(Evaluate(Q4, u));
remeq := PrimitiveNumerator(Coefficient(rem, 0));
Feq := PrimitiveNumerator(F);
Leq := PrimitiveNumerator(L);
Ceq := PrimitiveNumerator(C);
Peq := PrimitiveNumerator(P);
Q0eq := PrimitiveNumerator(Q0);
Fred := ReduceInU(Feq, qeq);
Gmain := ExactQuotient(Fred, (r+1)^3);

procedure Summary()
    print "Full M(2,12) [4,12] generic s=m^2 condition";
    print "variables: r,z,u,s";
    print "a =", a;
    print "w =", w;
    print "beta =", beta;
    PrintPolyInfo("Q4(u)", qeq);
    PrintPolyInfo("remainder", remeq);
    print "factor remainder:";
    print Factorization(remeq);
    PrintPolyInfo("P", Peq);
    PrintPolyInfo("Q0", Q0eq);
    PrintPolyInfo("L", Leq);
    PrintPolyInfo("C", Ceq);
    PrintPolyInfo("F_s", Feq);
    PrintPolyInfo("F_s_mod_Q4", Fred);
    PrintPolyInfo("G_s_main", Gmain);
    facF := Factorization(Feq);
    print "F_s_factor_count", #facF;
    for i in [1..#facF] do
        PrintPolyInfo(Sprintf("F_s_factor_%o", i), facF[i][1]);
        if #Terms(facF[i][1]) le 10 then
            print "F_s_factor_poly", i, facF[i][1];
        end if;
        print "F_s_factor_exp", i, facF[i][2];
    end for;
    facFred := Factorization(Fred);
    print "F_s_mod_Q4_factor_count", #facFred;
    for i in [1..#facFred] do
        PrintPolyInfo(Sprintf("F_s_mod_Q4_factor_%o", i), facFred[i][1]);
        if #Terms(facFred[i][1]) le 10 then
            print "F_s_mod_Q4_factor_poly", i, facFred[i][1];
        end if;
        print "F_s_mod_Q4_factor_exp", i, facFred[i][2];
    end for;
end procedure;

procedure WriteOutput(filename)
    out := Open(filename, "w");
    fprintf out, "# Full M(2,12) [4,12] generic s=m^2 condition\n";
    fprintf out, "# variables r,z,u,s\n";
    fprintf out, "# Q4(u)=0 and F_s=0 are the generic conditions; then t=-C/L and N=P*t+Q0.\n";
    fprintf out, "# Special branch L=C=0 must be checked separately.\n";
    fprintf out, "a = %o\n", a;
    fprintf out, "w = %o\n", w;
    fprintf out, "beta = %o\n", beta;
    fprintf out, "Q4u = %o\n", qeq;
    fprintf out, "remainder = %o\n", remeq;
    fprintf out, "P = %o\n", Peq;
    fprintf out, "Q0 = %o\n", Q0eq;
    fprintf out, "L = %o\n", Leq;
    fprintf out, "C = %o\n", Ceq;
    fprintf out, "F_s = %o\n", Feq;
    fprintf out, "F_s_mod_Q4 = %o\n", Fred;
    fprintf out, "G_s_main = %o\n", Gmain;
    delete out;
    print "wrote", filename;
end procedure;

if mode eq "summary" then
    Summary();
elif mode eq "write" then
    WriteOutput(output_file);
else
    error "unknown mode";
end if;

quit;
