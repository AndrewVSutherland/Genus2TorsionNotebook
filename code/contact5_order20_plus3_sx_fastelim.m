//////////////////////////////////////////////////////////////////////
// Fast generic elimination on the normalized cyclic [20]+3 cover.
//
//   g_s=[1+(1+s)x+2s*x^2]^2-4(1-s)x^5,
//   H^2-q^3=M*g_s,  q=x^2+U*x+V.
//
// After matching degrees 5,4,3, the equations G1 and G2 are quadratic
// in V.  Away from their proportional/degree-drop boundary, eliminating
// V^2 gives a linear equation D*V-N=0.  Substitute V=N/D into G1 and
// G0, then eliminate U.  This is much smaller than a four-variable
// primary decomposition and gives the generic (s,M)-projection.  The
// final exact-cyclic cover is M=L^2.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemGB then MemGB := 8; end if;
if Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
if not assigned PointBound then PointBound := 100; end if;
if Type(PointBound) eq MonStgElt then PointBound := StringToInteger(PointBound); end if;
if not assigned MaxPlaneDegree then MaxPlaneDegree := 24; end if;
if Type(MaxPlaneDegree) eq MonStgElt then MaxPlaneDegree := StringToInteger(MaxPlaneDegree); end if;
if not assigned RecoveryAudit then RecoveryAudit := false; end if;
if Type(RecoveryAudit) eq MonStgElt then
    RecoveryAudit := StringToLower(RecoveryAudit) in {"true","1","yes"};
end if;
if not assigned FullFunctionFieldRecovery then FullFunctionFieldRecovery := false; end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();
Z := Integers();
C<s,m,u> := PolynomialRing(Q,3,"grevlex");
PV<v> := PolynomialRing(C);

c0 := C!1;
c1 := 2*(s+1);
c2 := s^2+6*s+1;
c3 := 4*s*(s+1);
c4 := 4*s^2;
c5 := 4*(s-1);

A := (m*c5+3*u)/2;
B := (m*c4+3*(u^2+v)-A^2)/2;
E := (m*c3+u^3+6*u*v-2*A*B)/2;
G2 := B^2+2*A*E-3*(u^2*v+v^2)-m*c2;
G1 := 2*B*E-3*u*v^2-m*c1;
G0 := E^2-v^3-m;

assert Degree(G1) eq 2;
assert Degree(G2) eq 2;
assert Degree(G0) eq 3;

function PrimitiveC(f)
    // Keep the rational scaling.  Coefficients() on a multivariate
    // polynomial is variable-relative in Magma, so integer-content code
    // would needlessly recurse through coefficient polynomials.
    return C!f;
end function;

function ClearEvaluate(f,num,den)
    d := Degree(f);
    return PrimitiveC(&+[Coefficient(f,i)*num^i*den^(d-i)
                         : i in [0..d]]);
end function;

a2 := Coefficient(G1,2); a1 := Coefficient(G1,1); a0 := Coefficient(G1,0);
b2 := Coefficient(G2,2); b1 := Coefficient(G2,1); b0 := Coefficient(G2,0);

// b2*G1-a2*G2 = den*V-num.
denV := PrimitiveC(b2*a1-a2*b1);
numV := PrimitiveC(a2*b0-b2*a0);
assert b2*G1-a2*G2 eq PV!(denV*v-numV);

F1 := ClearEvaluate(G1,numV,denV);
F2 := ClearEvaluate(G2,numV,denV);
F0 := ClearEvaluate(G0,numV,denV);

// Exact quadratic common-root identity.  On denV!=0, V=numV/denV and
// G1=G2=0 are equivalent to Rcore=0.  This removes the possible
// different-U and leading-coefficient artifacts of pairwise resultants.
cCross := a1*b0-a0*b1;
Rcore := PrimitiveC(numV^2+denV*cCross);
assert F1 eq a2*Rcore;
assert F2 eq b2*Rcore;
assert b2*F1 eq a2*F2;

print "CONTACT5 ORDER20+3 NORMALIZED FAST ELIMINATION";
print "V_DEGREES",<Degree(G2),Degree(G1),Degree(G0)>;
print "V_LINEAR_RECOVERY","den_degree",TotalDegree(denV),
      "den_terms",#Terms(denV),"num_degree",TotalDegree(numV),
      "num_terms",#Terms(numV);
print "SUBSTITUTED_EQUATIONS",
      <TotalDegree(Rcore),#Terms(Rcore),Degree(Rcore,u)>,
      <TotalDegree(F1),#Terms(F1),Degree(F1,u)>,
      <TotalDegree(F2),#Terms(F2),Degree(F2,u)>,
      <TotalDegree(F0),#Terms(F0),Degree(F0,u)>;

// Controls: any solution of F1=F0 with the listed open denominators
// recovers a simultaneous solution of G2=G1=G0.
assert denV ne 0 and F1 ne 0 and F2 ne 0 and F0 ne 0;

print "EXACT_RESULTANT_U_BEGIN";
time ResExact := Resultant(Rcore,F0,u);
assert Degree(ResExact,u) eq 0;
print "EXACT_RESULTANT_U_DONE",<TotalDegree(ResExact),#Terms(ResExact)>;

print "PAIRWISE_RESULTANTS_CROSSCHECK_BEGIN";
time Res10 := Resultant(F1,F0,u);
time Res20 := Resultant(F2,F0,u);
// F1 and F2 satisfy the cleared identity b2*F1=a2*F2, so their
// mutual resultant is identically zero.  Intersect each with F0 and
// take the gcd of the two plane resultants.
assert Degree(Res10,u) eq 0 and Degree(Res20,u) eq 0;
print "PAIRWISE_RESULTANTS_CROSSCHECK_DONE",
      <TotalDegree(Res10),#Terms(Res10)>,
      <TotalDegree(Res20),#Terms(Res20)>;

S<ss,mm> := PolynomialRing(Q,2,"grevlex");
toS := hom<C -> S | ss,mm,S!0>;
R10sm := S!toS(Res10);
R20sm := S!toS(Res20);
Rexactsm := S!toS(ResExact);
Rpairsm := GCD(R10sm,R20sm);
assert IsDivisibleBy(Rpairsm,Rexactsm);
Rsm := Rexactsm;
print "PAIRWISE_RESULTANT_GCD","degree",TotalDegree(Rsm),
      "terms",#Terms(Rsm),"pairwise_extra_degree",
      TotalDegree(ExactQuotient(Rpairsm,Rexactsm));

// Strip only parameter boundaries here.  Factors caused by denV=0 or
// leading-coefficient drops are labelled later and retained unless their
// full support is known to be boundary.
Ds := 8*ss^3-59*ss^2-18*ss+197;
known := mm*(ss-1)*(ss-2)*Ds;
procedure StripFactor(~f,h,label)
    n := 0;
    while IsDivisibleBy(f,h) do
        f := ExactQuotient(f,h);
        n +:= 1;
    end while;
    if n gt 0 then print "STRIPPED",label,n; end if;
end procedure;
StripFactor(~Rsm,mm,"M=0");
StripFactor(~Rsm,ss-1,"s=1");
StripFactor(~Rsm,ss-2,"s=2");
StripFactor(~Rsm,Ds,"family cubic");

print "OPEN_RESULTANT","degree",TotalDegree(Rsm),"terms",#Terms(Rsm);
print "DENOMINATOR_BRANCH_RESULTANT_BEGIN";
time DenRes := Resultant(denV,numV,u);
assert Degree(DenRes,u) eq 0;
DenSm := S!toS(DenRes);
denOverlap := GCD(Rsm,DenSm);
print "DENOMINATOR_BRANCH_RESULTANT","degree",TotalDegree(DenSm),
      "terms",#Terms(DenSm),"open_overlap_degree",TotalDegree(denOverlap);
print "OPEN_FACTORIZATION_BEGIN";
time fac := Factorization(Rsm);
print "OPEN_FACTORIZATION_DONE","factors",#fac;

SL<SQ,LQ> := PolynomialRing(Q,2,"grevlex");
squaremap := hom<S -> SL | SQ,LQ^2>;
genericBases := [S | ];
for i in [1..#fac] do
    base := fac[i][1];
    onDenBoundary := IsDivisibleBy(DenSm,base);
    if not onDenBoundary then Append(~genericBases,base); end if;
    print "SM_FACTOR",i,"multiplicity",fac[i][2],
          "degree",TotalDegree(base),"degree_s",Degree(base,ss),
          "degree_m",Degree(base,mm),"terms",#Terms(base),
          "denV_numV_boundary",onDenBoundary;
    if #Terms(base) le 120 then print base; end if;
    pull := SL!squaremap(base);
    pfac := Factorization(pull);
    print " SQUARE_PULLBACK","degree",TotalDegree(pull),"factors",#pfac;
    for pe in pfac do
        pc := pe[1];
        print "  SQ_FACTOR","multiplicity",pe[2],
              "degree",TotalDegree(pc),"degree_s",Degree(pc,SQ),
              "degree_L",Degree(pc,LQ),"terms",#Terms(pc);
        if TotalDegree(pc) gt MaxPlaneDegree then continue; end if;
        try
            Caff := Curve(AffineSpace(SL),pc);
            Cp := ProjectiveClosure(Caff);
            print "   GEOMETRY","projective_degree",Degree(Cp),
                  "genus",Genus(Cp),"nonsingular",IsNonsingular(Cp);
            pts := Points(Cp : Bound := PointBound);
            print "   POINTS_BOUND",PointBound,"count",#pts;
            if #pts le 40 then print pts; end if;
        catch err
            print "   GEOMETRY_FAILED",err`Object;
        end try;
    end for;
end for;

//////////////////////////////////////////////////////////////////////
// Generic recovery audit.
//
// A plane resultant is only a necessary projection until the eliminated
// coordinates are recovered over its function field.  For every
// non-denominator factor P(s,m), view P as an irreducible polynomial in m
// over Q(s), form K=Q(s)[m]/(P), and compute gcd_u(F1,F0).  Degree one
// proves that U is a rational function on the plane curve; V=numV/denV is
// then rational too.  Direct substitution in all G_i closes the audit.
//////////////////////////////////////////////////////////////////////

if RecoveryAudit then
    print "GENERIC_RECOVERY_AUDIT_BEGIN","factors",#genericBases;
    for ib in [1..#genericBases] do
        base := genericBases[ib];
        function FactorValuation(poly,h)
            val := 0;
            while IsDivisibleBy(poly,h) do
                poly := ExactQuotient(poly,h);
                val +:= 1;
            end while;
            return val;
        end function;
        val10 := FactorValuation(R10sm,base);
        val20 := FactorValuation(R20sm,base);
        print " RECOVERY_RESULTANT_MULTIPLICITIES",ib,val10,val20,
              "generic_gcd_degree_one_certificate",Min(val10,val20) eq 1;
        assert Min(val10,val20) eq 1;
        if not FullFunctionFieldRecovery then
            // F1 and F2 are exactly proportional after V=N/D.  A simple
            // resultant divisor means that one of (F1,F0),(F2,F0) has
            // generic gcd degree one over Q(P), hence U is a rational
            // function there.  V=N/D is then rational as well.
            continue;
        end if;
        Ks<st> := FunctionField(Q);
        KM<mt> := PolynomialRing(Ks);
        mapS_KM := hom<S -> KM | KM!st,mt>;
        baseKM := mapS_KM(base);
        print " RECOVERY_FACTOR",ib,"degree_m",Degree(baseKM),
              "irreducible_over_Qs",IsIrreducible(baseKM);
        if not IsIrreducible(baseKM) then
            continue;
        end if;
        Ksm<mbar> := ext<Ks | baseKM>;
        KU<uu> := PolynomialRing(Ksm);
        mapC := hom<C -> KU | KU!(Ksm!st),KU!mbar,uu>;
        f1K := mapC(F1);
        f0K := mapC(F0);
        gu := GCD(f1K,f0K);
        print "  GENERIC_U_GCD degree",Degree(gu),"terms",#Coefficients(gu);
        if Degree(gu) ne 1 then
            continue;
        end if;
        urec := -Coefficient(gu,0)/Coefficient(gu,1);
        nK := Evaluate(mapC(numV),urec);
        dK := Evaluate(mapC(denV),urec);
        print "  GENERIC_DEN_NONZERO",dK ne 0;
        assert dK ne 0;
        vrec := nK/dK;

        function EvalV(F,uv,vv)
            return &+[Evaluate(mapC(C!Coefficient(F,j)),uv)*vv^j
                       : j in [0..Degree(F)]];
        end function;
        assert EvalV(G2,urec,vrec) eq 0;
        assert EvalV(G1,urec,vrec) eq 0;
        assert EvalV(G0,urec,vrec) eq 0;
        print "  GENERIC_UV_RECOVERY_PASS";
    end for;
    print "GENERIC_RECOVERY_AUDIT_DONE";
end if;

//////////////////////////////////////////////////////////////////////
// The exceptional 0/0 branch of V=numV/denV.
//
// Here G1 and G2 are proportional/degree-degenerate as quadratics in V.
// Impose denV=numV=0 and require a common V with G0 by both pairwise
// V-resultants.  This usually cuts the apparent degree-20 plane factor
// down to a zero-dimensional scheme.
//////////////////////////////////////////////////////////////////////

print "DEN_BRANCH_EXACT_BEGIN";
RV10 := C!Resultant(G1,G0);
RV20 := C!Resultant(G2,G0);
JdenRaw := ideal<C | denV,numV,RV10,RV20>;
Cb := m*(s-1)*(s-2)*(8*s^3-59*s^2-18*s+197);
time Jden := Saturation(JdenRaw,ideal<C | Cb>);
bd := Basis(Jden);
print "DEN_BRANCH_EXACT","dimension",Dimension(Jden),
      "basis_size",#bd,"basis_degrees",[TotalDegree(f) : f in bd];
if Dimension(Jden) eq 0 then
    try
        ptsden := Variety(Jden);
        print "DEN_BRANCH_RATIONAL_SMU",#ptsden,ptsden;
        // Recovery of V and the square test for m are done only if a
        // rational quotient point exists.
        for pp in ptsden do
            sv := Q!pp[1]; mv := Q!pp[2]; uv := Q!pp[3];
            print " DEN_POINT","s",sv,"m",mv,"u",uv,
                  "m_square",mv ge 0 and IsSquare(mv);
        end for;
    catch err
        print "DEN_BRANCH_VARIETY_FAILED",err`Object;
    end try;
end if;

print "FAST_ELIMINATION_DONE";
quit;
