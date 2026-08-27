//////////////////////////////////////////////////////////////////////
//  Exact A=0 slice of the cubic-contact square cover on
//
//    g_s=(1+(1+s)x+2*s*x^2)^2-4*(1-s)*x^5.
//
//  Since A=(4(s-1)M+3U)/2, this slice has
//
//    U=-(4/3)(s-1)M.
//
//  The full coefficient identity leaves a zero-dimensional ideal in
//  (s,M,V).  We saturate the exact open boundary, eliminate to s, recover
//  all rational quotient points, and require M=L^2.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(6*10^9);

Q := Rationals();
P<s,M,V> := PolynomialRing(Q,3,"grevlex");
PX<x> := PolynomialRing(P);
h := 1+(1+s)*x+2*s*x^2;
g := h^2-4*(1-s)*x^5;
c := [P!Coefficient(g,i) : i in [0..5]];
assert c eq [P!1,2*(s+1),s^2+6*s+1,4*s*(s+1),4*s^2,4*(s-1)];

U := -(Q!4/3)*(s-1)*M;
A := (M*c[6]+3*U)/2;
assert A eq 0;
B := (M*c[5]+3*(U^2+V)-A^2)/2;
E := (M*c[4]+U^3+6*U*V-2*A*B)/2;
q := x^2+U*x+V;
H := x^3+B*x+E;
identity := H^2-q^3-M*g;
G2raw := B^2-3*(U^2*V+V^2)-M*c[3];
G1raw := 2*B*E-3*U*V^2-M*c[2];
G0raw := E^2-V^3-M;
assert &and[Coefficient(identity,i) eq 0 : i in [3..6]];
assert [Coefficient(identity,i) : i in [0..2]] eq [G0raw,G1raw,G2raw];
print "CONTACT_A0_SELF_TEST_PASS";

function PrimitiveP(f)
    if f eq 0 then return f; end if;
    den := LCM([Denominator(a) : a in Coefficients(f)]);
    ints := [Integers()!(den*a) : a in Coefficients(f)];
    cont := GCD(ints); if cont eq 0 then cont := 1; end if;
    return P!((Q!den/Q!Abs(cont))*f);
end function;

G2 := PrimitiveP(G2raw); G1 := PrimitiveP(G1raw); G0 := PrimitiveP(G0raw);
discq := PrimitiveP(U^2-4*V);
resqg := PrimitiveP(P!Resultant(q,g));
Ds := 8*s^3-59*s^2-18*s+197;
boundary := M*(s-1)*(s-2)*Ds*discq*resqg;
Iraw := ideal<P | G2,G1,G0>;
print "A0 equation_degrees",[TotalDegree(f) : f in [G2,G1,G0]],
      "terms",[#Terms(f) : f in [G2,G1,G0]],
      "raw_dimension",Dimension(Iraw);
time Isat := Saturation(Iraw,ideal<P | boundary>);
print "A0_SAT dimension",Dimension(Isat),"basis_size",#Basis(Isat),
      "basis_degrees",[TotalDegree(f) : f in Basis(Isat)],
      "basis_terms",[#Terms(f) : f in Basis(Isat)];

RL<VL,ML,SL> := PolynomialRing(Q,3,"lex");
mp := hom<P -> RL | SL,ML,VL>;
Jraw := ideal<RL | mp(G2),mp(G1),mp(G0)>;
time J := Saturation(Jraw,ideal<RL | mp(boundary)>);
time gb := GroebnerBasis(J);
elim := [f : f in gb | Degree(f,1) eq 0 and Degree(f,2) eq 0];
print "A0_LEX dimension",Dimension(J),"basis_size",#gb,
      "variable_degrees",[[Degree(f,i) : i in [1..3]] : f in gb],
      "elim_generators",#elim;

QS<zS> := PolynomialRing(Q);
toQS := hom<RL -> QS | QS!0,QS!0,zS>;
spols := [QS!toQS(f) : f in elim | f ne 0];
if #spols eq 0 then print "A0_NO_S_ELIMINATION"; quit; end if;
s_core := spols[1];
for i in [2..#spols] do s_core := GCD(s_core,spols[i]); end for;
s_core := s_core/LeadingCoefficient(s_core);
fac := Factorization(s_core);
print "A0_S_ELIM degree",Degree(s_core),"terms",#Coefficients(s_core),
      "factors",#fac;
for fe in fac do
    print " S_FACTOR multiplicity",fe[2],"degree",Degree(fe[1]),
          "terms",#Coefficients(fe[1]),fe[1];
end for;
rat_s := Roots(s_core);
print "A0_RATIONAL_S_ROOTS",rat_s;

QMV<VV,MM> := PolynomialRing(Q,2,"lex");
candidates := [];
for sr in rat_s do
    s0 := sr[1];
    ev := hom<P -> QMV | QMV!s0,MM,VV>;
    K0 := ideal<QMV | ev(G2),ev(G1),ev(G0)>;
    kb := GroebnerBasis(K0);
    monly := [f : f in kb | Degree(f,1) eq 0 and f ne 0];
    if #monly eq 0 then print "A0_FIBER_POSDIM s",s0; continue; end if;
    QM<mvar> := PolynomialRing(Q);
    toM := hom<QMV -> QM | QM!0,mvar>;
    mpoly := QM!toM(monly[1]);
    for i in [2..#monly] do mpoly := GCD(mpoly,QM!toM(monly[i])); end for;
    for mr in Roots(mpoly) do
        m0 := mr[1];
        QV<vvar> := PolynomialRing(Q);
        toV := hom<QMV -> QV | vvar,QV!m0>;
        vpolys := [QV!toV(f) : f in [ev(G2),ev(G1),ev(G0)] | toV(f) ne 0];
        if #vpolys eq 0 then continue; end if;
        vp := vpolys[1];
        for i in [2..#vpolys] do vp := GCD(vp,vpolys[i]); end for;
        for vr in Roots(vp) do
            v0 := vr[1];
            if Evaluate(G2,[s0,m0,v0]) ne 0 or
               Evaluate(G1,[s0,m0,v0]) ne 0 or
               Evaluate(G0,[s0,m0,v0]) ne 0 then continue; end if;
            open := Evaluate(boundary,[s0,m0,v0]) ne 0;
            issq,l0 := IsSquare(m0);
            print "A0_RATIONAL_POINT s",s0,"M",m0,"V",v0,
                  "open",open,"M_square",issq;
            if open and issq then Append(~candidates,<s0,m0,v0,l0>); end if;
        end for;
    end for;
end for;
print "A0_OPEN_SQUARE_CANDIDATES",#candidates,candidates;

QX<X> := PolynomialRing(Q);
for tup in candidates do
    s0,m0,v0,l0 := Explode(tup);
    u0 := -(Q!4/3)*(s0-1)*m0;
    g0 := (1+(1+s0)*X+2*s0*X^2)^2-4*(1-s0)*X^5;
    cc := [Coefficient(g0,i) : i in [0..5]];
    a0 := (m0*cc[6]+3*u0)/2; assert a0 eq 0;
    b0 := (m0*cc[5]+3*(u0^2+v0))/2;
    e0 := (m0*cc[4]+u0^3+6*u0*v0)/2;
    q0 := X^2+u0*X+v0;
    H0 := X^3+b0*X+e0;
    assert H0^2-q0^3 eq m0*g0;
    C := HyperellipticCurve(g0); Jc := Jacobian(C);
    TG,phi := TorsionSubgroup(Jc); inv := Invariants(TG);
    t0 := (1+s0)/(1-s0);
    print "A0_EXACT_CANDIDATE s",s0,"t",t0,"M",m0,"L",l0,
          "U",u0,"V",v0,"g",g0,"q",q0,"H",H0,
          "torsion",inv,"cyclic60",inv eq [60];
end for;

print "A0_DONE";
quit;
