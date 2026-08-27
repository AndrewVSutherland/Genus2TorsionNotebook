//////////////////////////////////////////////////////////////////////
//  Symbolic [3,6] and [6,6] conditions on the contact-6 M(6) family.
//
//  The base family is
//
//      h6 = 1 + a*x + b*x^2 + x^3,
//      f  = h6^2 - (x - 1)^6.
//
//  Then f(0)=0 and deg(f)<=5.  For smooth f and h6(1) != 0, the point
//  P=(1,h6(1)) gives a marked class D=P-infinity of order dividing 6,
//  generically exact order 6.
//
//  The independent 3-torsion class is imposed by cubic contact
//
//      h3(x)^2 - f(x) = m^2*q(x)^3,
//      q = x^2 + U*x + v^2.
//
//  We set L=1/m and eliminate N,R,S from
//
//      h3 = m*x^3 + N*x^2 + R*x + S,    S=m*v^3.
//
//  The halving cover for [6,6] introduces
//
//      u = x^2 + alpha*x + beta,
//      ell = l3*x^3 + l2*x^2 + l1*x + l0,
//
//  and imposes
//
//      ell^2 - f = l3^2*u^2*q,
//      ell + h3 == 0 mod q.
//
//  The sign in the last congruence just chooses one of E or -E; since
//  E has order 3, existence of a rational half is unchanged by negation.
//
//  Typical runs:
//      magma code/contact6_m36_symbolic.m
//      magma -b mode:=write code/contact6_m36_symbolic.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned mode then
    mode := "summary";
end if;
if not assigned output_file then
    output_file := "data/contact6_m36_symbolic_equations.txt";
end if;

Q := Rationals();
R<a,b,L,U,v,alpha,beta,l0,l1,l2,l3> := PolynomialRing(Q, 11);
P<x> := PolynomialRing(R);

h6 := 1 + a*x + b*x^2 + x^3;
f := h6^2 - (x - 1)^6;

c1 := Coefficient(f, 1);
c2 := Coefficient(f, 2);
c3 := Coefficient(f, 3);
c4 := Coefficient(f, 4);
c5 := Coefficient(f, 5);

// Cubic-contact elimination.  With B=c5*L^2+3U and
// Delta=4*c4*L^2+12*(U^2+v^2)-B^2, one has
// N=B/(2L), R=Delta/(8L), S=v^3/L.
B := c5*L^2 + 3*U;
Delta := 4*c4*L^2 + 12*(U^2 + v^2) - B^2;

F3 := B*Delta + 16*v^3 - 8*c3*L^2 - 8*U^3 - 48*U*v^2;
F2 := Delta^2 + 64*B*v^3 - 64*c2*L^2
      - 192*(U^2*v^2 + v^4);
F1 := Delta*v^3 - 4*c1*L^2 - 12*U*v^4;

q := x^2 + U*x + v^2;
u := x^2 + alpha*x + beta;
ell := l3*x^3 + l2*x^2 + l1*x + l0;

halve_identity := ell^2 - f - l3^2*u^2*q;
assert Coefficient(halve_identity, 6) eq 0;
halve_eqs := [Coefficient(halve_identity, i) : i in [0..5]];

// Clear the common denominator 8L in ell + h3 mod q.
h3_clear := 8*x^3 + 4*B*x^2 + Delta*x + 8*v^3;
halve_congruence := (8*L*ell + h3_clear) mod q;
halve_rem_eqs := [Coefficient(halve_congruence, i) : i in [0..1]];

function SummaryLine(label, F)
    return Sprintf("%o: total_degree %o; degrees [a,b,L,U,v,alpha,beta,l0,l1,l2,l3] %o; terms %o",
                   label, TotalDegree(F),
                   [Degree(F, i) : i in [1..11]], #Terms(F));
end function;

procedure PrintSummary()
    print "Contact-6 M(6) base family";
    print "h6 =", h6;
    print "f  =", f;
    print "coefficients c5..c1 =", [c5,c4,c3,c2,c1];
    print "";
    print "[3,6] cubic-contact equations";
    print "B     =", B;
    print "Delta =", Delta;
    print SummaryLine("F3", F3);
    print SummaryLine("F2", F2);
    print SummaryLine("F1", F1);
    print "";
    print "[6,6] halving equations";
    for i in [1..#halve_eqs] do
        print SummaryLine(Sprintf("H%o", i), halve_eqs[i]);
    end for;
    for i in [1..#halve_rem_eqs] do
        print SummaryLine(Sprintf("R%o", i), halve_rem_eqs[i]);
    end for;
    print "";
    print "Nonboundary conditions:";
    print "  discriminant(f) != 0";
    print "  (b+3)*(a+b+2)*L*v*(U^2-4*v^2)*l3 != 0";
    print "  gcd(q,f)=1, gcd(u,f)=1";
    print "  the new 3-torsion class is independent from 2D";
    print "  exact [6,6] should be checked by TorsionSubgroup on specializations";
end procedure;

procedure WriteEquations(file)
    out := Open(file, "w");
    fprintf out, "Contact-6 M(6) base family\n";
    fprintf out, "h6 = %o\n", h6;
    fprintf out, "f = %o\n", f;
    fprintf out, "c5 = %o\nc4 = %o\nc3 = %o\nc2 = %o\nc1 = %o\n\n",
            c5, c4, c3, c2, c1;

    fprintf out, "[3,6] cubic-contact cover\n";
    fprintf out, "q = %o\n", q;
    fprintf out, "B = %o\n", B;
    fprintf out, "Delta = %o\n", Delta;
    fprintf out, "N = B/(2L), R = Delta/(8L), S = v^3/L\n";
    fprintf out, "F3 = %o\n", F3;
    fprintf out, "F2 = %o\n", F2;
    fprintf out, "F1 = %o\n\n", F1;

    fprintf out, "[6,6] halving cover\n";
    fprintf out, "u = %o\n", u;
    fprintf out, "ell = %o\n", ell;
    fprintf out, "Identity equations: coeffs of ell^2 - f - l3^2*u^2*q in degrees 0..5\n";
    for i in [1..#halve_eqs] do
        fprintf out, "H%o = %o\n", i, halve_eqs[i];
    end for;
    fprintf out, "Congruence equations: coeffs of 8L*ell + (8*x^3+4*B*x^2+Delta*x+8*v^3) mod q\n";
    for i in [1..#halve_rem_eqs] do
        fprintf out, "R%o = %o\n", i, halve_rem_eqs[i];
    end for;
    fprintf out, "\nNonboundary conditions:\n";
    fprintf out, "discriminant(f) != 0\n";
    fprintf out, "(b+3)*(a+b+2)*L*v*(U^2-4*v^2)*l3 != 0\n";
    fprintf out, "gcd(q,f)=1, gcd(u,f)=1\n";
    fprintf out, "E independent from 2D for [3,6]\n";
    fprintf out, "3H independent from 3D, or directly TorsionSubgroup has two 6-divisible invariants, for [6,6]\n";
    delete out;
    print "Wrote", file;
end procedure;

if mode eq "summary" then
    PrintSummary();
elif mode eq "write" then
    PrintSummary();
    WriteEquations(output_file);
else
    error "unknown mode";
end if;

quit;
