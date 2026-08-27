//////////////////////////////////////////////////////////////////////
// Exact lift audit for the three low (s,r)-resultant factors in the
// fixed-Weierstrass order-20 + rational-3 cover.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned mode then mode:="quadratics"; end if;
if mode notin {"quadratics","c5"} then error "mode quadratics or c5"; end if;
if not assigned MemGB then MemGB:=8; end if;
if Type(MemGB) eq MonStgElt then MemGB:=StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q:=Rationals();

procedure ReportIdeal(label,I)
    bs:=Basis(I);
    unit:=&or[g ne 0 and TotalDegree(g) eq 0:g in bs];
    print label,"basis",#bs,"unit",unit;
    if not unit then
        degout:=-1; try degout:=Degree(I); catch e degout:=-2; end try;
        print label,"dimension",Dimension(I),"degree",degout,
              "basis_degrees",[TotalDegree(g):g in bs],
              "basis_terms",[#Terms(g):g in bs];
        if Dimension(I) eq 0 then
            try print label,"Q_VARIETY",Variety(I); catch e
                print label,"Q_VARIETY_FAILED",e`Object;
            end try;
        end if;
    end if;
end procedure;

if mode eq "quadratics" then
    for eps in [1,-1] do
        // eps=1: s=r^2+2, factor r^2-s+2.
        // eps=-1: s=2-r^2, factor r^2+s-2.
        R<r,a,j>:=PolynomialRing(Q,3,"grevlex");
        s:=R!2+eps*r^2;
        F4 := r^5*j^4 - 4*s^2*r^3*j^2 - 6*r^3*a*j^2 + 20*r^3*j^2
              -4*s^2*j^2 -12*r^3+12*r^2*a-3*r*a^2
              +16*s*j^2-16*j^2;
        F3 := 2*s^3*r^3*j^2+s^2*r^2*j^4-8*s^2*r^3*j^2
              -3*r^4*a*j^2-4*s*r^2*j^4-2*s*r^3*j^2+4*r^2*j^4
              +20*r^3*j^2-3*s^2*a*j^2+4*r^4-12*r^3*a+9*r^2*a^2
              -2*r*a^3+12*s*a*j^2-12*a*j^2;
        F2 := -s^4*r^4*j^2+8*s^3*r^4*j^2-4*r^7*j^2-14*s^2*r^4*j^2
              +s^4*j^4-16*s*r^4*j^2-6*s^2*r^2*a*j^2-8*s^3*j^4
              -12*r^6+12*r^5*a-3*r^4*a^2+39*r^4*j^2
              +24*s*r^2*a*j^2+24*s^2*j^4-24*r^2*a*j^2
              -32*s*j^4+16*j^4;
        K:=FieldOfFractions(R); KZ<z>:=PolynomialRing(K);
        x:=z+1;
        Fs:=((1-s)^2+(1-s^2)*x+2*s*x^2)^2-4*x^5;
        q:=z^2+a*z+r^2;
        resqF:=R!Numerator(Resultant(q,Fs));
        Ds:=8*s^3-59*s^2-18*s+197;
        boundary:=r*j*(s-1)*(s-2)*Ds*(a^2-4*r^2)*resqF;
        I:=ideal<R|F4,F3,F2>;
        ReportIdeal(Sprintf("QUAD_%o_RAW",eps),I);
        print "QUAD",eps,"SATURATION_BEGIN";
        time Isat:=Saturation(I,ideal<R|boundary>);
        ReportIdeal(Sprintf("QUAD_%o_OPEN",eps),Isat);
    end for;
    print "FIXED_LOWFACTOR_QUADRATICS_DONE";
    quit;
end if;

R<s,r,a,j>:=PolynomialRing(Q,4,"grevlex");
F4 := r^5*j^4 - 4*s^2*r^3*j^2 - 6*r^3*a*j^2 + 20*r^3*j^2
      -4*s^2*j^2 -12*r^3+12*r^2*a-3*r*a^2+16*s*j^2-16*j^2;
F3 := 2*s^3*r^3*j^2+s^2*r^2*j^4-8*s^2*r^3*j^2
      -3*r^4*a*j^2-4*s*r^2*j^4-2*s*r^3*j^2+4*r^2*j^4
      +20*r^3*j^2-3*s^2*a*j^2+4*r^4-12*r^3*a+9*r^2*a^2
      -2*r*a^3+12*s*a*j^2-12*a*j^2;
F2 := -s^4*r^4*j^2+8*s^3*r^4*j^2-4*r^7*j^2-14*s^2*r^4*j^2
      +s^4*j^4-16*s*r^4*j^2-6*s^2*r^2*a*j^2-8*s^3*j^4
      -12*r^6+12*r^5*a-3*r^4*a^2+39*r^4*j^2+24*s*r^2*a*j^2
      +24*s^2*j^4-24*r^2*a*j^2-32*s*j^4+16*j^4;
C5:=s^4*r-4*s^2*r^3-8*s^3*r-8*r^4+14*s^2*r+20*r^3
    +8*s^2+16*s*r-32*s-39*r+32;
K:=FieldOfFractions(R); KZ<z>:=PolynomialRing(K); x:=z+1;
Fs:=((1-s)^2+(1-s^2)*x+2*s*x^2)^2-4*x^5;
q:=z^2+a*z+r^2;
resqF:=R!Numerator(Resultant(q,Fs));
Ds:=8*s^3-59*s^2-18*s+197;
boundary:=r*j*(s-1)*(s-2)*Ds*(a^2-4*r^2)*resqF;
I:=ideal<R|F4,F3,F2,C5>;
ReportIdeal("C5_RAW",I);
print "C5_SATURATION_BEGIN";
time Isat:=Saturation(I,ideal<R|boundary>);
ReportIdeal("C5_OPEN",Isat);
print "FIXED_LOWFACTOR_C5_DONE";
quit;
