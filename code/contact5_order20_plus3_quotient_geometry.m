//////////////////////////////////////////////////////////////////////
//  Exact cubic-contact quotient cover for the full one-parameter
//  contact-5/order-20 family, targeting exact cyclic Z/60.
//
//  This direct-t file is retained as an independent algebraic cross-check.
//  For actual component work prefer
//      code/contact5_order20_plus3_sx_geometry.m,
//  whose s=(t-1)/(t+1), x=(1-s)z normalization has coefficient degree <=2.
//
//      h_t = 1+t*x+(t^2-1)*x^2/2,
//      f_t = h_t^2-(t+1)^4*x^5/4.
//
//  Every smooth specialization has a marked order-20 point.  Generically
//  f_t=(x-1)*quartic with irreducible quartic, so there is exactly one
//  rational 2-direction and an independent 3-class is compatible with
//  exact cyclic Z/60.
//
//  A rational 3-class has a monic cubic-contact identity
//
//      H^2-q^3 = M*f_t,
//      q=x^2+U*x+V,
//      H=x^3+A*x^2+B*x+E,
//
//  where M=L^2 is a nonzero rational square.  Matching x^5,x^4,x^3
//  determines A,B,E recursively.  The remaining x^2,x^1,x^0
//  coefficients give three equations in (t,M,U,V), hence a curve.
//
//  This file keeps the quotient curve (M arbitrary) separate from the
//  final double cover M=L^2.  The repeated-support locus disc(q)=0 is a
//  legitimate 3-torsion possibility (Mumford support 2P), so it is split
//  off before the generic quotient is saturated by disc(q).  It cannot
//  yield exact cyclic Z/60: the identity factors f into degrees 2+3, and
//  the family's fixed rational root forces a third branch factor, hence
//  rational 2-rank at least two.  Repeat mode is only a [2,60] side lane.
//  The generic open also removes the singular t-boundary, M=0, and
//  resultant(q,f)=0.  Its modes are:
//
//      summary    print/assert the exact equations (default; cheap)
//      vzero      exact V=0 split-support derivation (cheap)
//      repeat     q=(x-rho)^2 with M=L^2, including rho=0
//      decompose  saturated Q-component decomposition
//      eliminate  lex elimination to the (t,M)-plane, then M=L^2
//      param       impose M=L^2 on V!=0 by rationally parametrizing
//                  E^2-L^2=V^3 and decompose the resulting curve
//
//  The V=0 branch is already solved exactly in this file.  With E=eps*L,
//  eps=+/-1, saturation by L*U gives
//
//      t=3/U,  L=eps*U^3/2,
//      (U+1)^3*(U^2+9U+24)=0.
//
//  The quadratic has discriminant -15.  The only rational solution is
//  U=-1, t=-3, L=-eps/2, which lies on the singular family boundary.
//  Thus V=0 contains no smooth rational Z/60 candidate.
//
//  Audited finite baseline (2026-07-10):
//
//    * exact F_p/F_{p^2} point counts are controlled by the known
//      20 | #J(F_p) at every good residue;
//    * every tested prime through 101 has good residues with 3|#J(F_p),
//      so the full t-line has no local emptiness obstruction;
//    * masks through p=47 leave 55 height-1000 residues, 53 smooth;
//    * the stronger 19-prime mask set p=7,...,79 checks all 1,216,767
//      reduced rationals of naive height <=1000 and leaves only t=-1,-3,
//      both singular.  No smooth point is omitted in that height box.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned mode then
    mode := "summary";
end if;
if mode notin {"summary","vzero","repeat","decompose","eliminate","param"} then
    error "mode must be summary, vzero, repeat, decompose, eliminate, or param";
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
R<t,mu,U,V> := PolynomialRing(Q,4,"grevlex");
PX<x> := PolynomialRing(R);

b := (t^2-1)/2;
h5 := 1+t*x+b*x^2;
f := h5^2-(t+1)^4*x^5/4;
c0 := R!Coefficient(f,0);
c1 := R!Coefficient(f,1);
c2 := R!Coefficient(f,2);
c3 := R!Coefficient(f,3);
c4 := R!Coefficient(f,4);
c5 := R!Coefficient(f,5);

assert c0 eq 1;
assert c1 eq 2*t;
assert c2 eq 2*t^2-1;
assert c3 eq t*(t^2-1);
assert c4 eq (t^2-1)^2/4;
assert c5 eq -(t+1)^4/4;

q3 := x^2+U*x+V;
A := (mu*c5+3*U)/2;
B := (mu*c4+3*(U^2+V)-A^2)/2;
E := (mu*c3+U^3+6*U*V-2*A*B)/2;
H := x^3+A*x^2+B*x+E;
identity := H^2-q3^3-mu*f;

// The top three coefficients vanish by construction.
assert Coefficient(identity,6) eq 0;
assert Coefficient(identity,5) eq 0;
assert Coefficient(identity,4) eq 0;
assert Coefficient(identity,3) eq 0;

G2raw := B^2+2*A*E-3*(U^2*V+V^2)-mu*c2;
G1raw := 2*B*E-3*U*V^2-mu*c1;
G0raw := E^2-V^3-mu*c0;
assert Coefficient(identity,2) eq G2raw;
assert Coefficient(identity,1) eq G1raw;
assert Coefficient(identity,0) eq G0raw;

function PrimitiveR(g)
    if g eq 0 then return g; end if;
    den := LCM([Denominator(cc) : cc in Coefficients(g)]);
    ints := [Z!(den*cc) : cc in Coefficients(g)];
    cont := GCD(ints);
    if cont eq 0 then cont := 1; end if;
    return R!((Q!den/Q!Abs(cont))*g);
end function;

G2 := PrimitiveR(G2raw);
G1 := PrimitiveR(G1raw);
G0 := PrimitiveR(G0raw);

discq := U^2-4*V;
resqf := PrimitiveR(R!Resultant(q3,f));
disc3 := 32*t^3+152*t^2+173*t+37;
family_boundary := (t+1)*(t+3)*disc3;
boundary := mu*discq*resqf*family_boundary;
Iraw := ideal<R | G2,G1,G0>;

procedure PrintIdeal(label,I,components)
    print label,"dimension",Dimension(I),"basis_size",#Basis(I),
          "basis_degrees",[TotalDegree(g) : g in Basis(I)];
    if components then
        time comps := PrimaryDecomposition(I);
        print label,"primary_components",#comps;
        for i in [1..#comps] do
            C := comps[i];
            print " COMPONENT",i,"dimension",Dimension(C),
                  "prime",IsPrime(C),"basis_size",#Basis(C),
                  "basis_degrees",[TotalDegree(g) : g in Basis(C)];
            if #Basis(C) le 8 then print Basis(C); end if;
        end for;
    end if;
end procedure;

procedure Summary()
    print "FULL CONTACT5/ORDER20 + 3 CUBIC-CONTACT QUOTIENT";
    print "f_t =",f;
    print "disc(f_t) proportional to (t+1)^11*(t+3)^4*disc3";
    print "disc3 =",disc3;
    print "q =",q3;
    print "A =",A;
    print "B =",B;
    print "E =",E;
    print "G2 degree",TotalDegree(G2),"terms",#Terms(G2),G2;
    print "G1 degree",TotalDegree(G1),"terms",#Terms(G1),G1;
    print "G0 degree",TotalDegree(G0),"terms",#Terms(G0),G0;
    print "disc(q) =",discq;
    print "resultant(q,f) degree",TotalDegree(resqf),"terms",#Terms(resqf);
    if #Terms(resqf) le 80 then print "resultant(q,f) =",resqf; end if;
    print "boundary degree",TotalDegree(boundary);
    print "RAW_QUOTIENT generator_degrees",
          [TotalDegree(g) : g in [G2,G1,G0]],
          "generator_terms",[#Terms(g) : g in [G2,G1,G0]];
end procedure;

//////////////////////////////////////////////////////////////////////
// Exact V=0 branch.
//////////////////////////////////////////////////////////////////////

procedure VZeroBranch()
    S<T,L,W> := PolynomialRing(Q,3,"grevlex");
    QU<Zeta> := PolynomialRing(Q);
    print "V=0 SPLIT-SUPPORT BRANCH";
    splitpoly := W^5+12*W^4+54*W^3+100*W^2+81*W+24;
    assert splitpoly eq (W+1)^3*(W^2+9*W+24);
    splitpolyQ := Zeta^5+12*Zeta^4+54*Zeta^3+100*Zeta^2+81*Zeta+24;
    print "derived relations t=3/U, L=eps*U^3/2";
    print "U polynomial =",splitpoly;
    print "factorization",Factorization(splitpolyQ);
    print "quadratic discriminant",9^2-4*24;
    print "rational roots",Roots(splitpolyQ);

    for eps in [-1,1] do
        mp := hom<R -> S | T,L^2,W,S!0>;
        H3 := S!(mp(E)-eps*L);
        H2 := S!mp(G2raw);
        H1 := S!mp(G1raw);
        qres := S!mp(resqf);
        sboundary := L*W*qres*(T+1)*(T+3)*
                     (32*T^3+152*T^2+173*T+37);
        Jraw := ideal<S | H3,H2,H1>;
        print " eps",eps,"raw";
        PrintIdeal("  V0_RAW",Jraw,false);
        J := Saturation(Jraw,ideal<S | sboundary>);
        PrintIdeal("  V0_SATURATED",J,true);

        // The only Q-root of the eliminated polynomial is the singular
        // point U=-1,t=-3,L=-eps/2; verify it in the unsaturated equations.
        pt := [Q!-3,Q!(-eps)/2,Q!-1];
        assert Evaluate(H3,pt) eq 0;
        assert Evaluate(H2,pt) eq 0;
        assert Evaluate(H1,pt) eq 0;
        print "  singular_control",pt,
              "identity_multiplier_M",(Q!(-eps)/2)^2;
    end for;
    print "V0_CONCLUSION only rational point has t=-3, hence singular";
end procedure;

//////////////////////////////////////////////////////////////////////
// Repeated-support branch disc(q)=0.
//////////////////////////////////////////////////////////////////////

procedure RepeatedSupportBranch()
    // A non-squarefree Mumford u-polynomial is not a degeneracy: it
    // represents the divisor 2P-2*infinity.  Put
    //
    //     q=(x-rho)^2, U=-2*rho, V=rho^2,
    //
    // and impose M=L^2 exactly.  The rho=0 sub-slice q=x^2 is first
    // eliminated algebraically: with E=eps*L, G1 gives B=eps*L*t and G2
    // gives A=eps*L*(t^2-1)/2, after which the recursive degree-3 formula
    // gives E=0, contradicting L!=0.  We may therefore saturate rho here.
    S<TS,MS,Rho,LS> := PolynomialRing(Q,4,"grevlex");
    mp := hom<R -> S | TS,MS,-2*Rho,Rho^2>;
    Jraw := ideal<S | mp(G2),mp(G1),mp(G0),LS^2-MS>;
    qres := S!mp(resqf);
    sboundary := MS*Rho*qres*(TS+1)*(TS+3)*
                 (32*TS^3+152*TS^2+173*TS+37);
    print "REPEATED_SUPPORT q=(x-rho)^2, M=L^2";
    print "rho=0 sub-slice eliminated: E_recursive=0 contradicts E=eps*L";
    print "raw generator degrees",
          [TotalDegree(g) : g in [mp(G2),mp(G1),mp(G0),LS^2-MS]];
    PrintIdeal("REPEATED_RAW",Jraw,false);
    time J := Saturation(Jraw,ideal<S | sboundary>);
    PrintIdeal("REPEATED_SATURATED",J,true);
    if Dimension(J) eq 0 then
        try
            pts := Variety(J);
            print "REPEATED_RATIONAL_POINTS",#pts,pts;
            for pt in pts do
                tv := Q!pt[1];
                mv := Q!pt[2];
                rv := Q!pt[3];
                lv := Q!pt[4];
                print " REPEATED_POINT","t",tv,"M",mv,
                      "rho",rv,"L",lv,"square_check",lv^2 eq mv;
            end for;
        catch err
            print "REPEATED_VARIETY_FAILED",err`Object;
        end try;
    end if;
end procedure;

//////////////////////////////////////////////////////////////////////
// Quotient decomposition and (t,M) elimination.
//////////////////////////////////////////////////////////////////////

procedure DecomposeQuotient()
    print "SATURATING QUOTIENT CURVE";
    time Isat := Saturation(Iraw,ideal<R | boundary>);
    PrintIdeal("SATURATED_QUOTIENT",Isat,true);
end procedure;

procedure EliminateToTM()
    // Put U,V first in lex order, so the polynomials free of U,V form the
    // exact elimination ideal in Q[t,M].
    RL<UL,VL,TL,ML> := PolynomialRing(Q,4,"lex");
    mp := hom<R -> RL | TL,ML,UL,VL>;
    Jraw := ideal<RL | mp(G2),mp(G1),mp(G0)>;
    J := Saturation(Jraw,ideal<RL | mp(boundary)>);
    print "LEX_SATURATED dimension",Dimension(J);
    time gb := Basis(J);
    print "LEX_BASIS size",#gb,"degrees",[TotalDegree(g) : g in gb];
    elim := [g : g in gb | Degree(g,1) eq 0 and Degree(g,2) eq 0];
    print "TM_ELIM generators",#elim;

    S<TM,MM> := PolynomialRing(Q,2);
    toS := hom<RL -> S | S!0,S!0,TM,MM>;
    plane := [S!toS(g) : g in elim | g ne 0];
    for i in [1..#plane] do
        print " TM_GENERATOR",i,"degree",TotalDegree(plane[i]),
              "terms",#Terms(plane[i]);
        if #Terms(plane[i]) le 100 then print plane[i]; end if;
    end for;
    if #plane eq 0 then
        print "TM_ELIMINATION_EMPTY_OUTPUT";
        return;
    end if;
    core := plane[1];
    for i in [2..#plane] do core := GCD(core,plane[i]); end for;
    print "TM_CURVE_GCD degree",TotalDegree(core),"terms",#Terms(core);
    fac := Factorization(core);
    print "TM_CURVE_FACTORS",#fac;

    SW<TW,LW> := PolynomialRing(Q,2);
    squaremap := hom<S -> SW | TW,LW^2>;
    for i in [1..#fac] do
        g := fac[i][1];
        print " TM_FACTOR",i,"multiplicity",fac[i][2],
              "degree",TotalDegree(g),"terms",#Terms(g);
        if #Terms(g) le 100 then print g; end if;
        sq := SW!squaremap(g);
        print "  SQUARE_PULLBACK degree",TotalDegree(sq),
              "factors",#Factorization(sq);
        for se in Factorization(sq) do
            sf := se[1];
            print "   SQ_FACTOR multiplicity",se[2],
                  "degree",TotalDegree(sf),"terms",#Terms(sf);
            if TotalDegree(sf) gt MaxPlaneDegree then continue; end if;
            try
                Caff := Curve(AffineSpace(SW),sf);
                Cp := ProjectiveClosure(Caff);
                print "    projective_degree",Degree(Cp),
                      "genus",Genus(Cp),"nonsingular",IsNonsingular(Cp);
                pts := Points(Cp : Bound := PointBound);
                print "    points_bound",PointBound,"count",#pts;
                if #pts le 40 then print "    points",pts; end if;
            catch err
                print "    plane_geometry_failed",err`Object;
            end try;
        end for;
    end for;
end procedure;

//////////////////////////////////////////////////////////////////////
// Direct rational parametrization of the M=L^2 constant equation.
//////////////////////////////////////////////////////////////////////

procedure ParametrizedSquareCover()
    // On V != 0, put w=(E+L)/V.  Since (E-L)(E+L)=V^3,
    //
    //   L=V*(w^2-V)/(2w), E=V*(w^2+V)/(2w),
    //   M=L^2=V^2*(w^2-V)^2/(4w^2).
    //
    // Impose E_recursive=E_param plus G2=G1=0.  G0 then holds
    // identically.  This is the square cover itself, not only its quotient.
    S<TS,US,VS,WS> := PolynomialRing(Q,4,"grevlex");
    K := FieldOfFractions(S);
    Lpar := K!(VS*(WS^2-VS)/(2*WS));
    Epar := K!(VS*(WS^2+VS)/(2*WS));
    mupar := Lpar^2;
    mp := hom<R -> K | K!TS,mupar,K!US,K!VS>;

    P3 := S!Numerator(mp(E)-Epar);
    P2 := S!Numerator(mp(G2raw));
    P1 := S!Numerator(mp(G1raw));
    assert S!Numerator(Epar^2-VS^3-mupar) eq 0;

    qres := S!Numerator(mp(resqf));
    sboundary := VS*WS*(WS^2-VS)*(US^2-4*VS)*qres*
                 (TS+1)*(TS+3)*(32*TS^3+152*TS^2+173*TS+37);
    Jraw := ideal<S | P3,P2,P1>;
    print "PARAM_SQUARE_COVER raw equation degrees",
          [TotalDegree(g) : g in [P3,P2,P1]],
          "terms",[#Terms(g) : g in [P3,P2,P1]];
    PrintIdeal("PARAM_RAW",Jraw,false);
    time J := Saturation(Jraw,ideal<S | sboundary>);
    PrintIdeal("PARAM_SATURATED",J,true);
end procedure;

if mode eq "summary" then
    Summary();
elif mode eq "vzero" then
    Summary();
    VZeroBranch();
elif mode eq "repeat" then
    Summary();
    RepeatedSupportBranch();
elif mode eq "decompose" then
    Summary();
    DecomposeQuotient();
elif mode eq "eliminate" then
    Summary();
    EliminateToTM();
elif mode eq "param" then
    Summary();
    ParametrizedSquareCover();
end if;

quit;
