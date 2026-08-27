//////////////////////////////////////////////////////////////////////
//  Low-degree exact cubic-contact cover for the full contact-5/order-20
//  family plus rational 3-torsion.
//
//  Starting from t, put
//
//      s=(t-1)/(t+1),       t=(1+s)/(1-s),
//      x_old=(1-s)*z.
//
//  Then, without any y-scaling,
//
//      g_s(z)=f_t((1-s)z)
//            =[1+(1+s)z+2s z^2]^2-4(1-s)z^5.
//
//  Its coefficients are
//
//      c0=1, c1=2(s+1), c2=s^2+6s+1,
//      c3=4s(s+1), c4=4s^2, c5=4(s-1).
//
//  Thus every coefficient has degree at most two in s, while the constant
//  square equation remains E^2-V^3=M=L^2.  The finite singular/boundary
//  fibers are s=1 (t=infinity/degree drop), s=2 (t=-3), and the cubic
//
//      Ds=8s^3-59s^2-18s+197.
//
//  The missing t=-1 fiber is s=infinity and is singular; it is recorded
//  explicitly rather than silently treated as part of the affine chart.
//
//  Cubic contact:
//
//      q=z^2+U z+V,
//      H=z^3+A z^2+B z+E,
//      H^2-q^3=M*g_s,       M=L^2 != 0.
//
//  The repeated-support branch q=(z-rho)^2 is a legitimate 3-class, but
//  cannot give exact cyclic Z/60: the identity factors g_s into rational
//  degrees 2+3, and the fixed rational Weierstrass root splits one factor
//  again, forcing rational 2-rank at least two.  It is retained as a
//  [2,60]-type side diagnostic; the exact-cyclic generic cover may then be
//  saturated by disc(q).
//
//  Modes:
//      summary    exact symbolic assertions and degrees (default)
//      vzero      exact V=0 factorization
//      repeat     q=(z-rho)^2 with M=L^2
//      quotient   saturated generic quotient in (s,M,U,V)
//      eliminate  eliminate U,V to (s,M), then pull back M=L^2
//      param       directly impose M=L^2 using the rational constant conic
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned mode then mode := "summary"; end if;
if mode notin {"summary","vzero","repeat","quotient","eliminate","param"} then
    error "mode must be summary, vzero, repeat, quotient, eliminate, or param";
end if;
if not assigned MemGB then
    MemGB := 8;
elif Type(MemGB) eq MonStgElt then
    MemGB := StringToInteger(MemGB);
end if;
if not assigned PointBound then
    PointBound := 100;
elif Type(PointBound) eq MonStgElt then
    PointBound := StringToInteger(PointBound);
end if;
if not assigned MaxPlaneDegree then
    MaxPlaneDegree := 16;
elif Type(MaxPlaneDegree) eq MonStgElt then
    MaxPlaneDegree := StringToInteger(MaxPlaneDegree);
end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();
Z := Integers();
R<s,mu,U,V> := PolynomialRing(Q,4,"grevlex");
PZ<z> := PolynomialRing(R);

h := 1+(1+s)*z+2*s*z^2;
g := h^2-4*(1-s)*z^5;
c0 := R!Coefficient(g,0);
c1 := R!Coefficient(g,1);
c2 := R!Coefficient(g,2);
c3 := R!Coefficient(g,3);
c4 := R!Coefficient(g,4);
c5 := R!Coefficient(g,5);
assert [c0,c1,c2,c3,c4,c5] eq
       [1,2*(s+1),s^2+6*s+1,4*s*(s+1),4*s^2,4*(s-1)];

Ds := 8*s^3-59*s^2-18*s+197;
assert R!Discriminant(g) eq 256*(1-s)^2*(2-s)^4*Ds;

q := z^2+U*z+V;
A := (mu*c5+3*U)/2;
B := (mu*c4+3*(U^2+V)-A^2)/2;
E := (mu*c3+U^3+6*U*V-2*A*B)/2;
H := z^3+A*z^2+B*z+E;
identity := H^2-q^3-mu*g;
assert &and[Coefficient(identity,i) eq 0 : i in [3..6]];

G2raw := B^2+2*A*E-3*(U^2*V+V^2)-mu*c2;
G1raw := 2*B*E-3*U*V^2-mu*c1;
G0raw := E^2-V^3-mu;
assert [Coefficient(identity,i) : i in [0..2]] eq [G0raw,G1raw,G2raw];

function PrimitiveR(poly)
    if poly eq 0 then return poly; end if;
    den := LCM([Denominator(cc) : cc in Coefficients(poly)]);
    ints := [Z!(den*cc) : cc in Coefficients(poly)];
    cont := GCD(ints);
    if cont eq 0 then cont := 1; end if;
    return R!((Q!den/Q!Abs(cont))*poly);
end function;

G2 := PrimitiveR(G2raw);
G1 := PrimitiveR(G1raw);
G0 := PrimitiveR(G0raw);
discq := U^2-4*V;
resqg := PrimitiveR(R!Resultant(q,g));
family_boundary := (s-1)*(s-2)*Ds;
generic_boundary := mu*discq*resqg*family_boundary;
Iraw := ideal<R | G2,G1,G0>;

procedure PrintIdeal(label,I,components)
    bb := Basis(I);
    print label,"dimension",Dimension(I),"basis_size",#bb,
          "basis_degrees",[TotalDegree(f) : f in bb];
    if not components then return; end if;
    time comps := PrimaryDecomposition(I);
    print label,"components",#comps;
    for i in [1..#comps] do
        C := comps[i];
        bc := Basis(C);
        print " COMPONENT",i,"dimension",Dimension(C),
              "prime",IsPrime(C),"basis_size",#bc,
              "basis_degrees",[TotalDegree(f) : f in bc];
        if #bc le 8 then print bc; end if;
    end for;
end procedure;

procedure Summary()
    print "CONTACT5/ORDER20 +3 LOW-DEGREE s,x CHART";
    print "s=(t-1)/(t+1), t=(1+s)/(1-s), x_old=(1-s)z";
    print "s=infinity is t=-1 and is singular";
    print "g_s =",g;
    print "coefficients",[c0,c1,c2,c3,c4,c5];
    print "discriminant",R!Discriminant(g);
    print "A =",A;
    print "B =",B;
    print "E =",E;
    print "G2 degree",TotalDegree(G2),"terms",#Terms(G2),G2;
    print "G1 degree",TotalDegree(G1),"terms",#Terms(G1),G1;
    print "G0 degree",TotalDegree(G0),"terms",#Terms(G0),G0;
    print "resultant(q,g) degree",TotalDegree(resqg),"terms",#Terms(resqg);
    if #Terms(resqg) le 80 then print resqg; end if;
end procedure;

//////////////////////////////////////////////////////////////////////
// V=0 with distinct split support q=z(z+U).
//////////////////////////////////////////////////////////////////////

procedure VZeroBranch()
    // With E=eps*L and L*U!=0, the equations collapse to
    //
    //   B=eps*L*(s+1), A=2*eps*L*s, U^3=2*eps*L,
    //   s=3/U-1,
    //   (U-1)^3*(2U^2+3U+3)=0.
    //
    // The quadratic discriminant is -15; U=1 gives s=2, the singular
    // t=-3 fiber.  The q=z^2 sub-slice U=0 is handled in repeat mode.
    QU<W> := PolynomialRing(Q);
    splitpoly := 2*W^5-3*W^4-2*W^2+6*W-3;
    assert splitpoly eq (W-1)^3*(2*W^2+3*W+3);
    print "VZERO relations s=3/U-1, L=eps*U^3/2";
    print "VZERO polynomial",splitpoly;
    print "VZERO factorization",Factorization(splitpoly);
    print "VZERO roots",Roots(splitpoly),"quadratic_discriminant",-15;
    singular_pt := [Q!2,Q!1/4,Q!1,Q!0];
    assert Evaluate(G2raw,singular_pt) eq 0;
    assert Evaluate(G1raw,singular_pt) eq 0;
    assert Evaluate(G0raw,singular_pt) eq 0;
    print "VZERO singular_control (s,M,U,V)",singular_pt;
    print "VZERO_CONCLUSION only rational U=1 gives s=2 (t=-3), singular";
end procedure;

//////////////////////////////////////////////////////////////////////
// Repeated-support q=(z-rho)^2.
//////////////////////////////////////////////////////////////////////

procedure RepeatedSupportBranch()
    // rho=0 gives q=z^2.  If M=L^2 and E=eps*L, then G1 and G2 give
    // B=eps*L*(s+1), A=2*eps*L*s.  The recursive formula gives E=0,
    // contradicting L!=0.  Hence rho=0 is empty and may be saturated.
    print "REPEATED_EXACT_CYCLIC_OBSTRUCTION";
    print "M*g=(H-(z-rho)^3)*(H+(z-rho)^3) has degrees 2+3;";
    print "the fixed rational root forces >=3 branch factors and 2-rank >=2";
    S<SS,MS,Rho,LS> := PolynomialRing(Q,4,"grevlex");
    mp := hom<R -> S | SS,MS,-2*Rho,Rho^2>;
    Jraw := ideal<S | mp(G2),mp(G1),mp(G0),LS^2-MS>;
    qres := S!mp(resqg);
    sboundary := MS*Rho*qres*(SS-1)*(SS-2)*
                 (8*SS^3-59*SS^2-18*SS+197);
    print "REPEATED q=(z-rho)^2; rho=0 eliminated by E=0 contradiction";
    print "raw degrees",[TotalDegree(f) : f in
          [mp(G2),mp(G1),mp(G0),LS^2-MS]];
    PrintIdeal("REPEATED_RAW",Jraw,false);
    time J := Saturation(Jraw,ideal<S | sboundary>);
    PrintIdeal("REPEATED_SATURATED",J,true);
    if Dimension(J) eq 0 then
        try
            pts := Variety(J);
            print "REPEATED_RATIONAL_POINTS",#pts,pts;
        catch err
            print "REPEATED_VARIETY_FAILED",err`Object;
        end try;
    end if;
end procedure;

//////////////////////////////////////////////////////////////////////
// Generic quotient and elimination to (s,M).
//////////////////////////////////////////////////////////////////////

procedure QuotientComponents()
    print "GENERIC QUOTIENT: repeated support already split off";
    time J := Saturation(Iraw,ideal<R | generic_boundary>);
    PrintIdeal("GENERIC_QUOTIENT",J,true);
end procedure;

procedure EliminateToSM()
    // U,V first in lex order; U,V-free basis elements generate the exact
    // elimination ideal in Q[s,M].
    RL<UL,VL,SL,ML> := PolynomialRing(Q,4,"lex");
    mp := hom<R -> RL | SL,ML,UL,VL>;
    Jraw := ideal<RL | mp(G2),mp(G1),mp(G0)>;
    time J := Saturation(Jraw,ideal<RL | mp(generic_boundary)>);
    gb := Basis(J);
    elim := [f : f in gb | Degree(f,1) eq 0 and Degree(f,2) eq 0];
    print "SM_LEX dimension",Dimension(J),"basis_size",#gb,
          "elim_generators",#elim;

    S<SP,MP> := PolynomialRing(Q,2);
    toS := hom<RL -> S | S!0,S!0,SP,MP>;
    plane := [S!toS(f) : f in elim | f ne 0];
    if #plane eq 0 then print "SM_ELIM_EMPTY"; return; end if;
    core := plane[1];
    for i in [2..#plane] do core := GCD(core,plane[i]); end for;
    print "SM_CORE degree",TotalDegree(core),"terms",#Terms(core),
          "factors",#Factorization(core);

    SW<SQ,LQ> := PolynomialRing(Q,2);
    squaremap := hom<S -> SW | SQ,LQ^2>;
    for fe in Factorization(core) do
        base := fe[1];
        print " SM_FACTOR multiplicity",fe[2],"degree",TotalDegree(base),
              "terms",#Terms(base);
        if #Terms(base) le 100 then print base; end if;
        pull := SW!squaremap(base);
        for pe in Factorization(pull) do
            pc := pe[1];
            print "  SQUARE_FACTOR multiplicity",pe[2],
                  "degree",TotalDegree(pc),"terms",#Terms(pc);
            if TotalDegree(pc) gt MaxPlaneDegree then continue; end if;
            try
                Caff := Curve(AffineSpace(SW),pc);
                Cp := ProjectiveClosure(Caff);
                print "   projective_degree",Degree(Cp),
                      "genus",Genus(Cp),"nonsingular",IsNonsingular(Cp);
                pts := Points(Cp : Bound := PointBound);
                print "   points_bound",PointBound,"count",#pts;
                if #pts le 40 then print pts; end if;
            catch err
                print "   geometry_failed",err`Object;
            end try;
        end for;
    end for;
end procedure;

//////////////////////////////////////////////////////////////////////
// Direct square cover via E^2-L^2=V^3.
//////////////////////////////////////////////////////////////////////

procedure ParametrizedSquareCover()
    // On V!=0 put w=(E+L)/V.  Then
    //
    //   L=V(w^2-V)/(2w), E=V(w^2+V)/(2w), M=L^2.
    //
    // Impose E_recursive=E_param and G2=G1=0; G0 is automatic.
    S<SS,US,VS,WS> := PolynomialRing(Q,4,"grevlex");
    K := FieldOfFractions(S);
    Lpar := K!(VS*(WS^2-VS)/(2*WS));
    Epar := K!(VS*(WS^2+VS)/(2*WS));
    mupar := Lpar^2;
    mp := hom<R -> K | K!SS,mupar,K!US,K!VS>;
    P3 := S!Numerator(mp(E)-Epar);
    P2 := S!Numerator(mp(G2raw));
    P1 := S!Numerator(mp(G1raw));
    assert S!Numerator(Epar^2-VS^3-mupar) eq 0;
    qres := S!Numerator(mp(resqg));
    sboundary := VS*WS*(WS^2-VS)*(US^2-4*VS)*qres*
                 (SS-1)*(SS-2)*(8*SS^3-59*SS^2-18*SS+197);
    Jraw := ideal<S | P3,P2,P1>;
    print "PARAM_SQUARE equation_degrees",
          [TotalDegree(f) : f in [P3,P2,P1]],
          "terms",[#Terms(f) : f in [P3,P2,P1]];
    PrintIdeal("PARAM_RAW",Jraw,false);
    time J := Saturation(Jraw,ideal<S | sboundary>);
    PrintIdeal("PARAM_SATURATED",J,true);
end procedure;

Summary();
if mode eq "vzero" then
    VZeroBranch();
elif mode eq "repeat" then
    RepeatedSupportBranch();
elif mode eq "quotient" then
    QuotientComponents();
elif mode eq "eliminate" then
    EliminateToSM();
elif mode eq "param" then
    ParametrizedSquareCover();
end if;

quit;
