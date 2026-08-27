//////////////////////////////////////////////////////////////////////
// Exact cubic-contact cover over the simple M(2,12) line.
//
// The M(2,12) chart is
//      a=(1-z^2)/(4*(r+1)).
// The simple [2,12] component a=(1-r)/4 is therefore z^2=r^2.
// We use z=r, for which the chosen root of T+1 is x=2.
//
// Besides x=2, the sextic has the rational root
//      u=2*r/(2-r).
// The Mobius coordinate X=(x-u)/(x-2), equivalently
//      x=(2*X-u)/(X-1),
// sends u to 0 and 2 to infinity.  Thus
//      F_r(X)=(X-1)^6 W_r((2X-u)/(X-1))
// is an odd quintic with F_r(0)=0.
//
// A nonzero rational order-3 direction has cubic contact
//      h3^2-F_r=m^2*q^3,
//      q=X^2+U*X+v^2,
//      h3=m*X^3+N*X^2+R*X+m*v^3.
// Put L=1/m and eliminate N,R.  The resulting exact cover is cut out
// by E3=E2=E1=0 below in (r,L,U,v).  It contains the built-in direction
// 4D; an independent direction must have a different q-support.
//
// This file only derives/writes equations.  It does not launch a search.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned mode then mode := "summary"; end if;
if not assigned output_file then
    output_file := "data/m612_m212_line_contact_cover.txt";
end if;

Q := Rationals();
K<r,L,U,v> := RationalFunctionField(Q,4);
P<X> := PolynomialRing(K);
FP := FieldOfFractions(P);

a := (1-r)/4;
T := a*X^2-X+r;
h := (X-r)*(T+1);
W := h^2+4*a*X^2*T*(T+1);

assert Evaluate(T+1,2) eq 0;
u := 2*r/(2-r);
assert Evaluate(W,u) eq 0;

xold := FP!((2*X-u)/(X-1));
Ffrac := FP!((X-1)^6)*Evaluate(W,xold);
assert Denominator(Ffrac) eq 1;
F := P!Numerator(Ffrac);
assert Degree(F) eq 5;
assert Coefficient(F,0) eq 0;

c1 := Coefficient(F,1);
c2 := Coefficient(F,2);
c3 := Coefficient(F,3);
c4 := Coefficient(F,4);
c5 := Coefficient(F,5);

// Cubic-contact triangular elimination, identical to the validated
// contact6_m36 and m2228_three_torsion formulas.
B := c5*L^2+3*U;
Delta := 4*c4*L^2+12*(U^2+v^2)-B^2;
E3 := B*Delta+16*v^3-8*c3*L^2-8*U^3-48*U*v^2;
E2 := Delta^2+64*B*v^3-64*c2*L^2
      -192*(U^2*v^2+v^4);
E1 := Delta*v^3-4*c1*L^2-12*U*v^4;

// The q-support of the built-in class 4D.  It follows by transforming
// q0=X0^2+r*X0+(r-1)/4 under X=1+(2-u)X0.
U0 := -2*(2*r^2-3*r+2)/(2-r);
v0 := r/(2-r);

procedure PrintSummary()
    print "simple M(2,12) line a=(1-r)/4, equivalently z^2=r^2";
    print "u =",u;
    print "F_r(X) =",F;
    print "built-in q support: U0 =",U0,"v0^2 =",v0^2;
    print "contact auxiliaries B =",B;
    print "contact auxiliaries Delta =",Delta;
    print "cleared equation summaries";
    for tup in [<"E3",E3>,<"E2",E2>,<"E1",E1>] do
        num := Numerator(tup[2]);
        print tup[1],"numerator total degree",TotalDegree(num),
              "terms",#Terms(num),"denominator",Denominator(tup[2]);
    end for;
    print "nonboundary: (r-1)*(r-2)*L*v != 0, Disc(F) != 0";
    print "independence: (U,v^2) != (U0,v0^2)";
end procedure;

procedure WriteEquations(file)
    out := Open(file,"w");
    fprintf out,"a = %o\n",a;
    fprintf out,"u = %o\n",u;
    fprintf out,"F = %o\n",F;
    fprintf out,"c5,c4,c3,c2,c1 = %o\n",[c5,c4,c3,c2,c1];
    fprintf out,"B = %o\nDelta = %o\n",B,Delta;
    fprintf out,"E3_num = %o\n",Numerator(E3);
    fprintf out,"E2_num = %o\n",Numerator(E2);
    fprintf out,"E1_num = %o\n",Numerator(E1);
    fprintf out,"U0 = %o\nv0 = %o\n",U0,v0;
    fprintf out,"Nonboundary: (r-1)*(r-2)*L*v*Disc(F) != 0\n";
    fprintf out,"Independent support: (U,v^2) != (U0,v0^2)\n";
    delete out;
end procedure;

PrintSummary();
if mode eq "write" then
    WriteEquations(output_file);
elif mode ne "summary" then
    error "mode must be summary or write";
end if;
quit;
