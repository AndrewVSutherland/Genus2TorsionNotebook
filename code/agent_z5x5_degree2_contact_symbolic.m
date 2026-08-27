//////////////////////////////////////////////////////////////////////
//  One rational contact-5 plus one degree-2 contact-5 class.
//
//  Work on an odd quintic model C: y^2=f(x).  A rational 5-contact
//  at x=0 has
//
//      f = h^2 - K*x^5,          deg(h) <= 2.
//
//  A degree-2 Mumford class D=[U,V], U=x^2+s*x+t, has order dividing
//  5 if there is a function y-H(x) with
//
//      H^2 - f = U^5.
//
//  Equivalently, after scaling y, seek
//
//      f = H^2 - U^5 = h^2 - K*x^5,
//
//  with H monic of degree 5.  The coefficients of H are forced through
//  degree 6 by the condition deg(H^2-U^5) <= 5; one constant parameter
//  remains.  The rational contact at 0 is then exactly the condition
//  that the quartic tail of H^2-U^5 is a square.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Q := Rationals();
R<s,t,m> := PolynomialRing(Q, 3, "grevlex");
K := FieldOfFractions(R);
P<x> := PolynomialRing(K);

U := x^2 + s*x + t;
U5 := U^5;

// H=x^5+a4*x^4+...+a1*x+m.  Cancel degrees 9,8,7,6
// in H^2-U^5 recursively.  The constant term m is free.
a4 := Coefficient(U5, 9)/2;
a3 := (Coefficient(U5, 8) - a4^2)/2;
a2 := (Coefficient(U5, 7) - 2*a4*a3)/2;
a1 := (Coefficient(U5, 6) - 2*a4*a2 - a3^2)/2;
H := x^5 + a4*x^4 + a3*x^3 + a2*x^2 + a1*x + m;
F := H^2 - U5;

assert Degree(F) le 5;

q := [Coefficient(F, i) : i in [0..5]];
q0 := q[1]; q1 := q[2]; q2 := q[3];
q3 := q[4]; q4 := q[5]; q5 := q[6];

// On the open chart q4 != 0, the quartic
// q4*x^4+q3*x^3+q2*x^2+q1*x+q0 is a square iff C1=C0=0.
C1 := R!Numerator(8*q4^2*q1 - 4*q4*q3*q2 + q3^3);
C0 := R!Numerator(64*q4^3*q0 - (4*q4*q2 - q3^2)^2);

print "# degree-2/contact-5 symbolic normal form";
print "U =", U;
print "H =", H;
print "";
print "# F = H^2-U^5 has degree", Degree(F);
for i in [0..5] do
    printf "F%o = %o\n", i, q[i+1];
end for;

print "";
print "# On q4 != 0 the quartic tail is a square iff:";
print "C1 = 8*q4^2*q1 - 4*q4*q3*q2 + q3^3 = 0";
print "C0 = 64*q4^3*q0 - (4*q4*q2 - q3^2)^2 = 0";
print "";
printf "q4 = %o\n", q4;
printf "C1 degree=%o total_degree=%o terms=%o\n", Degree(C1), TotalDegree(C1), #Terms(C1);
printf "C0 degree=%o total_degree=%o terms=%o\n", Degree(C0), TotalDegree(C0), #Terms(C0);

print "";
print "# factorization of the two square-tail equations";
print "Factorization(C1):";
for pair in Factorization(C1) do
    printf "  exponent %o degree %o total_degree %o terms %o: %o\n",
        pair[2], Degree(pair[1]), TotalDegree(pair[1]), #Terms(pair[1]), pair[1];
end for;
print "Factorization(C0):";
for pair in Factorization(C0) do
    printf "  exponent %o degree %o total_degree %o terms %o: %o\n",
        pair[2], Degree(pair[1]), TotalDegree(pair[1]), #Terms(pair[1]), pair[1];
end for;

print "";
print "# gcd and a small elimination diagnostic";
G := GCD(C1, C0);
print "GCD(C1,C0) =", G;
print "Eliminating m by a resultant in the open chart can be expensive,";
print "so this script records only the degrees unless do_resultant:=true.";

if assigned do_resultant then
    if Type(do_resultant) eq MonStgElt then
        do_resultant := (do_resultant eq "true") or (do_resultant eq "True");
    end if;
else
    do_resultant := false;
end if;

if do_resultant then
    time Resm := Resultant(C1, C0, m);
    printf "Resultant_m degree=%o total_degree=%o terms=%o\n",
        Degree(Resm), TotalDegree(Resm), #Terms(Resm);
    print "Factorization(Resultant_m):";
    for pair in Factorization(Resm) do
        printf "  exponent %o degree %o total_degree %o terms %o: %o\n",
            pair[2], Degree(pair[1]), TotalDegree(pair[1]), #Terms(pair[1]), pair[1];
    end for;

    print "";
    print "# branch diagnostics for rational resultant factors";

    procedure PrintFactorizationOrZero(F, indent)
        if F eq 0 then
            printf "%oidentically zero\n", indent;
        else
            for pair in Factorization(F) do
                printf "%oexponent %o: %o\n", indent, pair[2], pair[1];
            end for;
        end if;
    end procedure;

    Rb<S,M> := PolynomialRing(Q, 2, "grevlex");
    phi_double := hom<R -> Rb | S, S^2/4, M>;
    phi_54 := hom<R -> Rb | S, 5*S^2/4, M>;

    print "Branch t=s^2/4 (U has a double root):";
    print "  C1 factors:";
    PrintFactorizationOrZero(phi_double(C1), "    ");
    print "  C0 factors:";
    PrintFactorizationOrZero(phi_double(C0), "    ");

    print "Branch t=5*s^2/4:";
    print "  C1 factors:";
    PrintFactorizationOrZero(phi_54(C1), "    ");
    print "  C0 factors:";
    PrintFactorizationOrZero(phi_54(C0), "    ");

    Rbeta<S2,B> := PolynomialRing(Q, 2, "grevlex");
    psi_54 := hom<R -> Rbeta | S2, 5*S2^2/4, B*S2^5>;
    print "  with m=B*s^5, C1 factors:";
    PrintFactorizationOrZero(psi_54(C1), "    ");
    print "  with m=B*s^5, C0 factors:";
    PrintFactorizationOrZero(psi_54(C0), "    ");

    Rz<Z> := PolynomialRing(Q);
    P4 := 1 - 8*Z + 16/5*Z^2;
    P16 := 1 - 144/5*Z + 8256/25*Z^2 - 240384/125*Z^3
        + 3708416/625*Z^4 - 5419008/625*Z^5
        + 5586944/3125*Z^6 + 19070976/3125*Z^7
        + 23658496/15625*Z^8;
    print "Dehomogenized non-linear ratio factors in Z=t/s^2:";
    print "  P4 factorization:", Factorization(P4);
    print "  P16 factorization:", Factorization(P16);
end if;

print "";
print "# Recover h=e*x^2+d*x+c on q4 != 0 by";
print "e^2=q4, 2*e*d=q3, 2*c*d=q1, c^2=q0.";
print "The search script checks the square condition directly, including";
print "q4=0 boundary cases, before constructing Jacobian classes.";

quit;
