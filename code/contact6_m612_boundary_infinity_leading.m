//////////////////////////////////////////////////////////////////////
// Exact leading and first-order equations at base-projective infinity
// for the contact-6 cubic-contact core.
//
// In the B=1 chart [A:B:T]=[t:1:e] use a=t/e,b=1/e; in the
// A=1 chart use a=1/e,b=t/e.  The base degrees of (F1,F2,F3) are
// (2,4,3).  On the nondegenerate contact chart the leading equations
// force M=L^2=1.  We print the two compatibility polynomials for a
// first-order transverse deformation at L=+1.  L=-1 gives the same
// conditions because the core equations depend only on L^2.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q:=Rationals();
R<e,t,L,U,v>:=PolynomialRing(Q,5,"grevlex");
K:=FieldOfFractions(R);

function CoreEquations(a,b,L,U,v)
    M:=L^2;
    c1:=2*a+6; c2:=a^2+2*b-15; c3:=2*a*b+22;
    c4:=2*a+b^2-15; c5:=2*b+6;
    B3:=c5*M+3*U;
    Delta3:=4*c4*M+12*(U^2+v^2)-B3^2;
    F3:=B3*Delta3+16*v^3-8*c3*M-8*U^3-48*U*v^2;
    F2:=Delta3^2+64*B3*v^3-64*c2*M
        -192*(U^2*v^2+v^4);
    F1:=Delta3*v^3-4*c1*M-12*U*v^4;
    return [F1,F2,F3];
end function;

S<tt,UU,vv>:=PolynomialRing(Q,3,"grevlex");
atL1:=hom<R -> S | S!0,tt,S!1,UU,vv>;

print "CONTACT6_M612_BOUNDARY_INFINITY_LEADING";
procedure PrintFactorOrZero(label,g)
    if g eq 0 then print label,"0";
    else print label,Factorization(g); end if;
end procedure;
for label in ["infB","infA"] do
    if label eq "infB" then
        aa:=(K!t)/(K!e); bb:=(K!1)/(K!e);
        desc:="[A:B:T]=[t:1:e]";
    else
        aa:=(K!1)/(K!e); bb:=(K!t)/(K!e);
        desc:="[A:B:T]=[1:t:e]";
    end if;
    raw:=CoreEquations(aa,bb,K!L,K!U,K!v);
    degs:=[2,4,3]; G:=[];
    for i in [1..3] do
        gg:=(K!e)^degs[i]*raw[i];
        assert Denominator(gg) eq 1;
        Append(~G,R!Numerator(gg));
    end for;
    leading:=[Coefficient(g,e,0):g in G];
    print "CHART",label,desc;
    print " LEADING_FACTORS",[Factorization(g):g in leading];
    assert &and[Evaluate(g,[Q!0,Q!0,Q!1,Q!0,Q!0]) eq 0
                :g in leading];

    h:=[S!atL1(Coefficient(g,e,1)):g in G];
    dL:=[S!atL1(Derivative(g,L)):g in G];
    // The middle leading equation has a double zero at L^2=1, so its
    // first-order equation is h[2]=0.  The other two share one correction
    // L1; eliminate it by h1*dL3-h3*dL1=0.
    assert dL[2] eq 0;
    compat2:=h[2];
    compat13:=h[1]*dL[3]-h[3]*dL[1];
    PrintFactorOrZero(" FIRST_ORDER_MIDDLE",compat2);
    PrintFactorOrZero(" FIRST_ORDER_COMPAT_13",compat13);
    PrintFactorOrZero(" FIRST_ORDER_IDEAL_GCD",GCD(compat2,compat13));
    if label eq "infA" then
        S0<LL,UU0,vv0>:=PolynomialRing(Q,3,"grevlex");
        atEndpoint:=hom<R -> S0 | S0!0,S0!0,LL,UU0,vv0>;
        endpoint_h:=[S0!atEndpoint(Coefficient(g,e,1)):g in G];
        print " ENDPOINT_[1:0:0]_FIRST_ORDER";
        for i in [1..3] do
            PrintFactorOrZero(Sprintf("  H%o",i),endpoint_h[i]);
        end for;
    end if;
end for;
print "CONTACT6_M612_BOUNDARY_INFINITY_LEADING_DONE";
quit;
