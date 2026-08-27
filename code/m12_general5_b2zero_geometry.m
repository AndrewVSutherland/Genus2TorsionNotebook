//////////////////////////////////////////////////////////////////////
// General Mumford order-5 locus on the compact M(12) chart, b2=0.
//
//   L = b+(2b-1)x,
//   H = x+w(1+b*x),
//   F = L*(L*H^2+4*b*(1+x)^2*(w*L-x^2)),
//
// and
//
//   A^2-B^2*F = q^5,
//   q=x^2+U*x+V,
//   B=b0+b1*x,
//   A=x^5+a4*x^4+...+a0.
//
// The coefficients a4,...,a0 are forced, in that order, by the
// coefficients of x^9,...,x^5.  The remaining five coefficients give
// five equations in (b,w,U,V,b0,b1).  This is the degree-one-B chart;
// b1 is inverted, but b0 is deliberately NOT inverted.
//
// Modes:
//   summary    exact equations, factors, boundaries, regression tests;
//   local      exhaustive open point count over GF(p), with a 5-torsion
//              Mumford check at each printed sample;
//   eliminate  exact saturation and projection to the (b,w)-plane;
//   slice      exact saturated sections b0=0 and b1=1.
//
// Examples (from torsion_jac):
//   magma -b mode:="summary" code/m12_general5_b2zero_geometry.m
//   magma -b mode:="local" p:=7 code/m12_general5_b2zero_geometry.m
//   magma -b mode:="eliminate" code/m12_general5_b2zero_geometry.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned mode then mode := "summary"; end if;
if mode notin {"summary","local","eliminate","slice"} then
    error "mode must be summary, local, eliminate, or slice";
end if;
if not assigned p then p := 7;
elif Type(p) eq MonStgElt then p := StringToInteger(p); end if;
if not assigned sample_limit then sample_limit := 8;
elif Type(sample_limit) eq MonStgElt then
    sample_limit := StringToInteger(sample_limit);
end if;
if not assigned MemGB then MemGB := 8;
elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

function Primitive(poly)
    if poly eq 0 then return poly; end if;
    if not IsFinite(BaseRing(Parent(poly))) then
        den := LCM([Denominator(c) : c in Coefficients(poly)]);
        nums := [Integers()!(den*c) : c in Coefficients(poly)];
        cont := GCD(nums); if cont eq 0 then cont := 1; end if;
        return Parent(poly)!((Rationals()!den/Rationals()!Abs(cont))*poly);
    end if;
    return poly;
end function;

function SpecializeUnivariate(G,vals,PP)
    xx:=PP.1;
    ans:=PP!0;
    for i in [0..Degree(G)] do
        cc:=Coefficient(G,i);
        num:=Numerator(cc); den:=Denominator(cc);
        ans +:= (Evaluate(num,vals)/Evaluate(den,vals))*xx^i;
    end for;
    return ans;
end function;

function ForcedA(F,q,B)
    P := Parent(F); x := P.1;
    A := x^5;
    coeffs := [];
    for d in [9,8,7,6,5] do
        E := A^2-B^2*F-q^5;
        aa := -Coefficient(E,d)/2;
        Append(~coeffs,aa);
        A +:= aa*x^(d-5);
    end for;
    E := A^2-B^2*F-q^5;
    assert &and[Coefficient(E,d) eq 0 : d in [5..10]];
    return A,E,coeffs;
end function;

function NormData(k : order := "grevlex")
    R<b,w,U,V,b0,b1> := PolynomialRing(k,6,order);
    K := FieldOfFractions(R);
    P<x> := PolynomialRing(K);
    L := b+(2*b-1)*x;
    H := x+w*(1+b*x);
    F := L*(L*H^2+4*b*(1+x)^2*(w*L-x^2));
    q := x^2+U*x+V;
    B := b0+b1*x;
    A,E,acoeffs := ForcedA(F,q,B);
    residuals := [R!Coefficient(E,i) : i in [0..4]];

    // Exact open chart.  b0 is not a boundary: B=b1*x is a valid
    // degree-one polynomial whenever q(0)=V is nonzero.
    discF := Primitive(R!Discriminant(F));
    discq := Primitive(R!Discriminant(q));
    resBq := Primitive(R!Resultant(B,q));
    resqF := Primitive(R!Resultant(q,F));
    boundary_factors := [b,w,b-1,2*b-1,b1,discF,discq,resBq,resqF];
    boundary_names := ["b","w","b-1","2b-1","b1","discF",
                       "discq","resBq","resqF"];
    boundary := &*boundary_factors;

    // Complete reconstruction test: after the five low residuals there
    // are no other coefficients left in the norm difference.
    assert E eq &+[K!residuals[i+1]*x^i : i in [0..4]];
    return R,P,F,q,B,A,E,acoeffs,residuals,boundary_factors,
           boundary_names,boundary;
end function;

procedure PrintSizes(label, polys)
    print label;
    for i in [1..#polys] do
        print " EQ",i-1,"degree",TotalDegree(polys[i]),
              "terms",#Terms(polys[i]);
    end for;
end procedure;

function IsUnitIdeal(I)
    return &or[g ne 0 and TotalDegree(g) eq 0 : g in Basis(I)];
end function;

procedure IdealReport(label,I)
    print label,"basis_len",#Basis(I),"unit",IsUnitIdeal(I);
    if not IsUnitIdeal(I) then
        try
            d,ds := Dimension(I);
            print label,"dimension",d,"component_degrees",ds;
        catch e
            print label,"dimension_failed",e`Object;
        end try;
        try print label,"degree",Degree(I);
        catch e print label,"degree_unavailable",e`Object; end try;
    end if;
end procedure;

procedure PrintFactorShape(label,g)
    fac := Factorization(g);
    print label,"factor_count",#fac,
          "shapes",[<TotalDegree(fe[1]),#Terms(fe[1]),fe[2]> : fe in fac];
    for fe in fac do
        if #Terms(fe[1]) le 12 then print label,"FACTOR",fe; end if;
    end for;
end procedure;

if mode eq "summary" then
    Q := Rationals();
    R,P,F,q,B,A,E,ac,res,bfs,bns,boundary := NormData(Q);
    x := P.1; b:=R.1; w:=R.2; U:=R.3; V:=R.4;
    b0:=R.5; b1:=R.6;
    print "M12_GENERAL5_B2ZERO_SELF_TEST_PASS";
    print "F_degree",Degree(F),"F_lc",R!Numerator(LeadingCoefficient(F));
    PrintSizes("FORCED_A_COEFFICIENT_SIZES",
               [R!Numerator(aa) : aa in ac]);
    PrintSizes("RESIDUAL_SIZES",res);
    g := res[1]; for i in [2..#res] do g:=GCD(g,res[i]); end for;
    print "RESIDUAL_GCD",g;
    for i in [1..#res] do PrintFactorShape(Sprintf("E%o",i-1),res[i]); end for;
    print "BOUNDARY_FACTORS";
    for i in [1..#bfs] do
        print " BOUNDARY",bns[i],"degree",TotalDegree(bfs[i]),
              "terms",#Terms(bfs[i]);
        if #Terms(bfs[i]) le 12 then print "  POLY",bfs[i]; end if;
    end for;
    PrintFactorShape("DISC_F",bfs[6]);
    print "B0_NOT_INVERTED","resBq_at_b0=0",Evaluate(bfs[8],[b,w,U,V,0,b1]);

    // Regression for the eight finite-sieve false positives.  They all
    // lie on the same multiplicity-six discriminant component.
    candidates := [
      <Q!-2,Q!1/3>, <Q!-1,Q!1/2>, <Q!2,Q!-1>, <Q!3,Q!-1/2>,
      <Q!-1/2,Q!2/3>, <Q!3/2,Q!-2>, <Q!1/3,Q!3/2>, <Q!2/3,Q!3>
    ];
    for bw in candidates do
        lin := bw[1]*bw[2]-bw[2]+1;
        assert lin eq 0;
        PQ<z>:=PolynomialRing(Q);
        fsp := SpecializeUnivariate(F,
                    [bw[1],bw[2],Q!0,Q!0,Q!0,Q!0],PQ);
        assert Discriminant(fsp) eq 0;
        print "SINGULAR_SIEVE_CONTROL",bw,"bw-w+1",lin;
    end for;
    print "M12_GENERAL5_B2ZERO_SUMMARY_DONE";
    quit;
end if;

if mode eq "local" then
    require IsPrime(p) and p ne 2: "p must be an odd prime";
    k := GF(p);
    R,P,F,q,B,A,E,ac,res,bfs,bns,boundary := NormData(k);
    elts := [a:a in k]; raw:=0; open:=0; smooth:=0; samples:=0;
    for bv in elts do for wv in elts do for uv in elts do
    for vv in elts do for b0v in elts do for b1v in elts do
        vals := [bv,wv,uv,vv,b0v,b1v];
        if &or[Evaluate(e,vals) ne 0:e in res] then continue; end if;
        raw +:= 1;
        if &or[Evaluate(h,vals) eq 0:h in bfs] then continue; end if;
        open +:= 1;
        JM := Matrix(k,5,6,
              [Evaluate(Derivative(res[i],j),vals):i in [1..5],j in [1..6]]);
        if Rank(JM) eq 5 then smooth +:= 1; end if;
        if samples lt sample_limit then
            samples +:= 1;
            PP<xx>:=PolynomialRing(k);
            ff:=SpecializeUnivariate(F,vals,PP);
            qq:=SpecializeUnivariate(q,vals,PP);
            BB:=SpecializeUnivariate(B,vals,PP);
            AA:=SpecializeUnivariate(A,vals,PP);
            assert AA^2-BB^2*ff eq qq^5;
            vvpoly := (-AA*InverseMod(BB,qq)) mod qq;
            assert (vvpoly^2-ff) mod qq eq 0;
            C:=HyperellipticCurve(ff); J:=Jacobian(C); D5:=J![qq,vvpoly];
            assert D5 ne J!0 and 5*D5 eq J!0;
            print "LOCAL_SAMPLE",vals,"jac_rank",Rank(JM),
                  "q",qq,"v",vvpoly,"D5_order",Order(D5);
        end if;
    end for; end for; end for; end for; end for; end for;
    print "LOCAL_COUNT","p",p,"raw",raw,"open",open,
          "smooth_open",smooth,"samples",samples;
    print "M12_GENERAL5_B2ZERO_LOCAL_DONE";
    quit;
end if;

Q:=Rationals();
R,P,F,q,B,A,E,ac,res,bfs,bns,boundary := NormData(Q);
Iraw := ideal<R|res>;
IdealReport("RAW",Iraw);

if mode eq "slice" then
    for spec in [<"b0=0",5,0>,<"b1=1",6,1>] do
        gen:=R.(spec[2]);
        Is := Iraw+ideal<R|gen-spec[3]>;
        // Saturating by the full open product is exact even when a factor
        // contains the specialized variable.
        Js := Saturation(Is,ideal<R|boundary>);
        IdealReport(Sprintf("SLICE_%o_SAT",spec[1]),Js);
        try
            comps:=PrimaryDecomposition(Js);
            print "SLICE_COMPONENTS",spec[1],#comps;
            for j in [1..#comps] do
                print " COMPONENT",j,"basis",Basis(comps[j]);
            end for;
        catch e print "SLICE_PRIMARY_FAILED",spec[1],e`Object; end try;
    end for;
    print "M12_GENERAL5_B2ZERO_SLICE_DONE";
    quit;
end if;

// Saturate sequentially: this is usually much cheaper and records exactly
// which open condition removes which dominant component.
J:=Iraw;
for i in [1..#bfs] do
    print "SATURATE_BEGIN",bns[i];
    time J:=Saturation(J,ideal<R|bfs[i]>);
    IdealReport(Sprintf("SAT_%o",bns[i]),J);
end for;

// Eliminate (U,V,b0,b1), retaining the parameter plane (b,w).  A lex
// ring with hidden variables first makes the projection ideal visible as
// those Groebner-basis elements involving only its last two variables.
RL<Ul,Vl,b0l,b1l,bl,wl> := PolynomialRing(Q,6,"lex");
mp:=hom<R->RL|bl,wl,Ul,Vl,b0l,b1l>;
JL:=ideal<RL|[mp(g):g in Basis(J)]>;
print "LEX_GROEBNER_BEGIN";
time Glex:=Basis(JL);
proj:=[g:g in Glex|Degree(g,Ul) eq 0 and Degree(g,Vl) eq 0 and
                         Degree(g,b0l) eq 0 and Degree(g,b1l) eq 0];
print "PROJECTION_BW_COUNT",#proj;
for g in proj do
    gp:=Primitive(g);
    print "PROJECTION_BW","degree",TotalDegree(gp),
          "bidegree",<Degree(gp,bl),Degree(gp,wl)>,"terms",#Terms(gp);
    PrintFactorShape("PROJECTION_BW_FACTOR",gp);
end for;
print "M12_GENERAL5_B2ZERO_ELIMINATE_DONE";
quit;
