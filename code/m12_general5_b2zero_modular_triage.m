//////////////////////////////////////////////////////////////////////
// Memory-bounded modular triage for the b2=0 order-5 norm cover on
// compact M(12).
//
// This file deliberately avoids the six-variable characteristic-zero
// Groebner basis.  It treats two disjoint loci.
//
// generic:
//   b1 != 0.  Use the root-of-B/sign quotient from
//   m12_general5_b2zero_rootquotient.m.  Its generic chart is three
//   equations K2,K3,K4 in (b,w,c,d).  We saturate first by the cheap
//   chart factors and only then by Disc(F) and Res(q,F), reducing each
//   large boundary polynomial modulo the current Groebner basis first.
//
// constant:
//   b1 = 0, b0 != 0.  Put s=b0^2 and form
//       A^2-s*F-q^5=0.
//   This gives five equations in (b,w,U,V,s), expected dimension zero.
//   The quotient points that lift to signed constant B are precisely
//   those with nonzero square s over the finite field.
//
// Typical capped calls (from the repository root):
//   magma -b p:=7 locus:="generic" MemGB:=6 \
//       code/m12_general5_b2zero_modular_triage.m
//   magma -b p:=7 locus:="constant" MemGB:=6 \
//       code/m12_general5_b2zero_modular_triage.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned p then p := 7;
elif Type(p) eq MonStgElt then p := StringToInteger(p); end if;
if not assigned locus then locus := "generic"; end if;
if locus notin {"generic","constant"} then
    error "locus must be generic or constant";
end if;
if not assigned generic_chart then generic_chart := "five"; end if;
if generic_chart notin {"five","four"} then
    error "generic_chart must be five or four";
end if;
if not assigned MemGB then MemGB := 6;
elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
if not assigned do_primary then do_primary := true;
elif Type(do_primary) eq MonStgElt then
    do_primary := do_primary in {"true","True","1"};
end if;
if not assigned full_open then full_open := true;
elif Type(full_open) eq MonStgElt then
    full_open := full_open in {"true","True","1"};
end if;
if not assigned print_low_degree then print_low_degree := 8;
elif Type(print_low_degree) eq MonStgElt then
    print_low_degree := StringToInteger(print_low_degree);
end if;
if not assigned sample_limit then sample_limit := 8;
elif Type(sample_limit) eq MonStgElt then
    sample_limit := StringToInteger(sample_limit);
end if;
if assigned fix_b and Type(fix_b) eq MonStgElt then
    fix_b := StringToInteger(fix_b);
end if;
if assigned fix_w and Type(fix_w) eq MonStgElt then
    fix_w := StringToInteger(fix_w);
end if;
require (assigned fix_b) eq (assigned fix_w):
        "fix_b and fix_w must be supplied together";

require IsPrime(p) and p notin {2,5}: "p must be prime and different from 2,5";
SetMemoryLimit(MemGB*10^9);
k := GF(p);

function IsUnitIdeal(I)
    return &or[g ne 0 and TotalDegree(g) eq 0 : g in Basis(I)];
end function;

procedure IdealReport(label,I,geometry)
    B := Basis(I);
    print "IDEAL",label,"basis_len",#B,"unit",IsUnitIdeal(I),
          "basis_degrees",[TotalDegree(g):g in B];
    if IsUnitIdeal(I) or not geometry then return; end if;
    try
        d,ds := Dimension(I);
        print "IDEAL",label,"dimension",d,"component_degrees",ds;
    catch err
        print "IDEAL",label,"dimension_failed",err`Object;
    end try;
    try
        print "IDEAL",label,"degree",Degree(I);
    catch err
        print "IDEAL",label,"degree_failed",err`Object;
    end try;
end procedure;

function SequentialSaturation(I,names,factors)
    J := I;
    for i in [1..#factors] do
        f := factors[i];
        if f eq 0 then
            print "SAT_SKIP_ZERO",names[i];
            continue;
        end if;
        B := Basis(J);
        nf := NormalForm(f,B);
        if nf eq 0 then
            print "SAT_FACTOR_IN_IDEAL",names[i],
                  "original_degree",TotalDegree(f);
        else
            print "SAT_BEGIN",names[i],"original_degree",TotalDegree(f),
                  "original_terms",#Terms(f),"nf_degree",TotalDegree(nf),
                  "nf_terms",#Terms(nf);
        end if;
        t0 := Cputime();
        J := Saturation(J,ideal<Generic(Parent(f))|nf>);
        print "SAT_END",names[i],"cpu",Cputime(t0),
              "basis_len",#Basis(J),"unit",IsUnitIdeal(J);
        if IsUnitIdeal(J) then break; end if;
    end for;
    return J;
end function;

procedure ComponentReport(label,J)
    if IsUnitIdeal(J) or not do_primary then
        print "PRIMARY_SKIPPED",label;
        return;
    end if;
    print "PRIMARY_BEGIN",label;
    t0 := Cputime();
    try
        comps := PrimaryDecomposition(J);
        print "PRIMARY_END",label,"cpu",Cputime(t0),"count",#comps;
        for i in [1..#comps] do
            C := comps[i]; BC:=Basis(C);
            cdim := -99; cdeg := -1;
            try
                cdim,cds := Dimension(C);
            catch err
                cdim := -99;
            end try;
            try
                cdeg := Degree(C);
            catch err
                cdeg := -1;
            end try;
            print "COMPONENT",label,i,"dimension",cdim,"degree",cdeg,
                  "basis_len",#BC,
                  "basis_degrees",[TotalDegree(g):g in BC],
                  "basis_terms",[#Terms(g):g in BC];
            if cdeg ge 0 and cdeg le print_low_degree and
               #BC le 12 and &and[#Terms(g) le 60:g in BC] then
                print "LOW_DEGREE_BASIS",label,i,BC;
            end if;
        end for;
    catch err
        print "PRIMARY_FAILED",label,"cpu",Cputime(t0),err`Object;
    end try;
end procedure;

function GenericData4(k)
    S<b,w,c,d> := PolynomialRing(k,4,"grevlex");
    PS<Z> := PolynomialRing(S);
    x := d*Z-c;
    L := b+(2*b-1)*x;
    H := x+w*(1+b*x);
    G := PS!(L*(L*H^2+4*b*(1+x)^2*(w*L-x^2)));
    gg := [S!Coefficient(G,i):i in [0..5]];
    G0:=gg[1]; G1:=gg[2]; G2:=gg[3]; G3:=gg[4];
    G4:=gg[5]; G5:=gg[6];
    h := G0+G1-G5;
    L0 := G0-G5;

    K3 := (38*h-30*G0)*L0-40*G0*(h-G4);
    K2 := (128*h^2-330*h*G0+1450*G0^2)*L0
          -1250*G0^3-150*G0*h^2+500*G0^2*G3;
    K4 := (95*h-75*G0)*L0^2
          -10*(25*G0^2+3*h^2)*L0+100*G0*G2*L0
          +(5*G0-h)^3;

    cheap_names := ["b","w","b-1","2b-1","d","G0","L0",
                    "h-5G0","h+5G0"];
    cheap := [b,w,b-1,2*b-1,d,G0,L0,h-5*G0,h+5*G0];

    // Full generic open conditions not already represented above.
    // Disc(q) is a unit multiple of (h-5G0)(h+5G0).
    X := Z;
    LX := b+(2*b-1)*X;
    HX := X+w*(1+b*X);
    Fbase := PS!(LX*(LX*HX^2+4*b*(1+X)^2*(w*LX-X^2)));
    discF := S!Discriminant(Fbase);
    qscaled := (5*G0)*Z^2+2*h*Z+5*G0;
    resqF := S!Resultant(qscaled,G);
    full_names := ["discF","resqF"];
    full := [discF,resqF];
    return S,[K2,K3,K4],cheap_names,cheap,full_names,full;
end function;

function GenericData5(k)
    R<b,w,c,d,e> := PolynomialRing(k,5,"grevlex");
    P<Z> := PolynomialRing(R);
    x := d*Z-c;
    L := b+(2*b-1)*x;
    H := x+w*(1+b*x);
    G := P!(L*(L*H^2+4*b*(1+x)^2*(w*L-x^2)));
    gg := [R!Coefficient(G,i):i in [0..5]];
    g0:=gg[1]; g1:=gg[2]; g2:=gg[3]; g3:=gg[4];
    g4:=gg[5]; g5:=gg[6];
    ell3 := (5*e/2)*g0-g1;
    ell4 := (5/2+15*e^2/8)*g0-g2;
    ell5 := (5/2+15*e^2/8)*g0-g3;
    ell6 := (5*e/2)*g0-g4;
    ell7 := g0-g5;
    delta := (2-e)^3;
    C1 := ell3-ell7;
    C2 := (32*e^2-33*e+58)*ell3-20*ell5;
    C3 := (19*e-6)*ell3-8*ell6;
    C4 := 4*(19*e-6)*ell3^2-32*ell4*ell3+5*delta*g0^2;

    cheap_names := ["b","w","b-1","2b-1","d","e-2","e+2","ell3"];
    cheap := [b,w,b-1,2*b-1,d,e-2,e+2,ell3];

    X := Z;
    LX := b+(2*b-1)*X;
    HX := X+w*(1+b*X);
    Fbase := P!(LX*(LX*HX^2+4*b*(1+X)^2*(w*LX-X^2)));
    discF := R!Discriminant(Fbase);
    qnorm := Z^2+e*Z+1;
    resqF := R!Resultant(qnorm,G);
    full_names := ["discF","resqF"];
    full := [discF,resqF];
    return R,[C1,C2,C3,C4],cheap_names,cheap,full_names,full;
end function;

function ConstantData(k)
    R<b,w,U,V,s> := PolynomialRing(k,5,"grevlex");
    K := FieldOfFractions(R);
    P<x> := PolynomialRing(K);
    L := b+(2*b-1)*x;
    H := x+w*(1+b*x);
    F := L*(L*H^2+4*b*(1+x)^2*(w*L-x^2));
    q := x^2+U*x+V;
    A := x^5;
    for deg in [9,8,7,6,5] do
        E := A^2-s*F-q^5;
        A +:= (-Coefficient(E,deg)/2)*x^(deg-5);
    end for;
    E := A^2-s*F-q^5;
    res := [R!Coefficient(E,i):i in [0..4]];
    assert E eq &+[K!res[i+1]*x^i:i in [0..4]];
    cheap_names := ["b","w","b-1","2b-1","s"];
    cheap := [b,w,b-1,2*b-1,s];
    full_names := ["discF","discq","resqF"];
    full := [R!Discriminant(F),R!Discriminant(q),R!Resultant(q,F)];
    return R,res,cheap_names,cheap,full_names,full;
end function;

procedure CountConstantPoints(R,res,boundaries)
    if p gt 13 then
        print "CONSTANT_BRUTE_SKIPPED","p",p,"reason","p>13";
        return;
    end if;
    elts := [a:a in k]; raw:=0; open:=0; signed:=0; smooth:=0;
    samples:=[];
    for bv in elts do for wv in elts do for uv in elts do
    for vv in elts do for sv in elts do
        vals := [bv,wv,uv,vv,sv];
        if &or[Evaluate(f,vals) ne 0:f in res] then continue; end if;
        raw +:= 1;
        if &or[Evaluate(f,vals) eq 0:f in boundaries] then continue; end if;
        open +:= 1;
        JM := Matrix(k,5,5,
              [Evaluate(Derivative(res[i],j),vals):i in [1..5],j in [1..5]]);
        if Rank(JM) eq 5 then smooth +:= 1; end if;
        if IsSquare(sv) then signed +:= 1; end if;
        if #samples lt sample_limit then
            Append(~samples,<vals,Rank(JM),IsSquare(sv)>);
        end if;
    end for; end for; end for; end for; end for;
    print "CONSTANT_POINT_COUNT","p",p,"raw",raw,"full_open",open,
          "smooth_full_open",smooth,"signed_lifts",signed,"samples",samples;
end procedure;

print "M12_GENERAL5_B2ZERO_MODULAR_TRIAGE_BEGIN","p",p,"locus",locus,
      "generic_chart",generic_chart,"MemGB",MemGB,"full_open",full_open,
      "do_primary",do_primary;

if locus eq "generic" then
    if generic_chart eq "five" then
        R,res,cheap_names,cheap,full_names,full := GenericData5(k);
    else
        R,res,cheap_names,cheap,full_names,full := GenericData4(k);
    end if;
    print "GENERIC_EQUATION_SHAPES",
          [<TotalDegree(f),#Terms(f)>:f in res];
    print "GENERIC_FULL_BOUNDARY_SHAPES",
          [<full_names[i],TotalDegree(full[i]),#Terms(full[i])>:
           i in [1..#full]];
    I := ideal<R|res>;
    if assigned fix_b then
        I +:= ideal<R|R.1-k!fix_b,R.2-k!fix_w>;
        print "GENERIC_FIXED_BASE","b",k!fix_b,"w",k!fix_w;
    end if;
    IdealReport("generic_raw",I,false);
    J := SequentialSaturation(I,cheap_names,cheap);
    IdealReport("generic_cheap_open",J,false);
    if full_open and not IsUnitIdeal(J) then
        J := SequentialSaturation(J,full_names,full);
        IdealReport("generic_full_open",J,true);
    end if;
    ComponentReport("generic_final",J);
else
    R,res,cheap_names,cheap,full_names,full := ConstantData(k);
    print "CONSTANT_EQUATION_SHAPES",
          [<TotalDegree(f),#Terms(f)>:f in res];
    print "CONSTANT_FULL_BOUNDARY_SHAPES",
          [<full_names[i],TotalDegree(full[i]),#Terms(full[i])>:
           i in [1..#full]];
    CountConstantPoints(R,res,cheap cat full);
    I := ideal<R|res>;
    IdealReport("constant_raw",I,false);
    J := SequentialSaturation(I,cheap_names,cheap);
    IdealReport("constant_cheap_open",J,false);
    if full_open and not IsUnitIdeal(J) then
        J := SequentialSaturation(J,full_names,full);
        IdealReport("constant_full_open",J,true);
    end if;
    ComponentReport("constant_final",J);
end if;

print "M12_GENERAL5_B2ZERO_MODULAR_TRIAGE_DONE","p",p,"locus",locus;
quit;
