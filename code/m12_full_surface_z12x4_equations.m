//////////////////////////////////////////////////////////////////////
//  Symbolic [4,12] conditions on the full M(2,12) surface.
//
//  M(2,12) chart:
//      a = (1-z^2)/(4*(r+1)),  T = a*x^2 - x + r,
//      W = ((x-r)(T+1))^2 + 4*a*x^2*T*(T+1).
//
//  A rational root u of Q4 = W/(T+1) gives an independent rational
//  2-torsion class.  Choose w = 2*(r+1)/(1+z), a root of T+1, and move
//  it to infinity: x = w + 1/X.  Then beta = 1/(u-w).
//
//  For f(X)=(X-beta)g(X), a half of [beta-infinity] is equivalent to
//      g(X) - (X-beta)*(mX+n)^2 = c*(X^2 + A X + B)^2,
//  where c is the leading coefficient of g.
//
//  The equations are therefore:
//      Q4(u)=0,
//      coeff_i(g - (X-beta)*(mX+n)^2 - c*(X^2+A*X+B)^2)=0, i=0..3.
//
//  Modes:
//      summary      print degrees/term counts and selected factorizations
//      write        write explicit equations to output_file
//      eliminate    compute a Groebner/elimination basis in a bounded way
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned mode then
    mode := "summary";
end if;
if not assigned output_file then
    output_file := "data/m12_full_surface_z12x4_equations.txt";
end if;

Q := Rationals();
R<r,z,u,A,B,m,n> := PolynomialRing(Q, 7, "grevlex");
K := FieldOfFractions(R);
Px<x> := PolynomialRing(K);
PX<X> := PolynomialRing(K);

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
c := Coefficient(g, 4);
U := X^2 + A*X + B;
E := g - (X-beta)*(m*X+n)^2 - c*U^2;

qeq := PrimitiveNumerator(Evaluate(Q4, u));
remeq := PrimitiveNumerator(Coefficient(rem, 0));
halving_eqs := [PrimitiveNumerator(Coefficient(E, i)) : i in [0..3]];
eqs := [qeq, remeq] cat halving_eqs;

// Reduced halving equations: solve the X^3 and X^2 coefficients for A and B.
g4 := Coefficient(g, 4);
g3 := Coefficient(g, 3);
g2 := Coefficient(g, 2);
g1 := Coefficient(g, 1);
g0 := Coefficient(g, 0);
Ared := (g3 - m^2)/(2*g4);
Bred := (g2 - 2*m*n + beta*m^2 - g4*Ared^2)/(2*g4);
R1 := g1 - n^2 + 2*beta*m*n - 2*g4*Ared*Bred;
R0 := g0 + beta*n^2 - g4*Bred^2;
red_eqs := [PrimitiveNumerator(R1), PrimitiveNumerator(R0)];


procedure PrintPolyInfo(label, f)
    print label,
          "total_degree", TotalDegree(f),
          "degree_r", Degree(f, 1),
          "degree_z", Degree(f, 2),
          "degree_u", Degree(f, 3),
          "degree_A", Degree(f, 4),
          "degree_B", Degree(f, 5),
          "degree_m", Degree(f, 6),
          "degree_n", Degree(f, 7),
          "terms", #Terms(f);
end procedure;

procedure Summary()
    print "Full M(2,12) [4,12] symbolic equations";
    print "variables: r,z,u,A,B,m,n";
    print "a =", a;
    print "w =", w;
    print "beta =", beta;
    PrintPolyInfo("Q4(u)", qeq);
    PrintPolyInfo("remainder", remeq);
    print "factor remainder:";
    print Factorization(remeq);
    for i in [1..#halving_eqs] do
        PrintPolyInfo(Sprintf("H%o", i-1), halving_eqs[i]);
    end for;
    PrintPolyInfo("R1_reduced", red_eqs[1]);
    PrintPolyInfo("R0_reduced", red_eqs[2]);
    print "factor R1_reduced:";
    print Factorization(red_eqs[1]);
    print "factor R0_reduced:";
    print Factorization(red_eqs[2]);
    print "factor Q4(u):";
    print Factorization(qeq);
    print "factor H0:";
    print Factorization(halving_eqs[1]);
    print "factor H3:";
    print Factorization(halving_eqs[4]);
end procedure;

procedure WriteEquations(filename)
    out := Open(filename, "w");
    fprintf out, "# Full M(2,12) [4,12] equations\n";
    fprintf out, "# variables r,z,u,A,B,m,n\n";
    fprintf out, "a = %o\n", a;
    fprintf out, "w = %o\n", w;
    fprintf out, "beta = %o\n", beta;
    fprintf out, "Q4u = %o\n", qeq;
    fprintf out, "remainder = %o\n", remeq;
    for i in [1..#halving_eqs] do
        fprintf out, "H%o = %o\n", i-1, halving_eqs[i];
    end for;
    fprintf out, "Ared = %o\n", Ared;
    fprintf out, "Bred = %o\n", Bred;
    fprintf out, "R1 = %o\n", red_eqs[1];
    fprintf out, "R0 = %o\n", red_eqs[2];
    delete out;
    print "wrote", filename;
end procedure;

procedure EliminateLight()
    print "Computing Groebner basis for the five [4,12] equations.";
    print "This may be expensive; use as a diagnostic, not the default path.";
    I := ideal<R | eqs>;
    G := GroebnerBasis(I);
    print "basis_length", #G;
    for i in [1..#G] do
        PrintPolyInfo(Sprintf("G%o", i), G[i]);
    end for;
end procedure;

if mode eq "summary" then
    Summary();
elif mode eq "write" then
    WriteEquations(output_file);
elif mode eq "eliminate" then
    EliminateLight();
else
    error "unknown mode";
end if;

quit;
