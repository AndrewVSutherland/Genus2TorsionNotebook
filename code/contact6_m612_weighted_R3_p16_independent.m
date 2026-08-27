//////////////////////////////////////////////////////////////////////
// Independent normalization diagnostics for the endpoint R3 P16
// projection component.  This deliberately avoids descent/rank work.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q:=Rationals();
A<e,mu,nu>:=PolynomialRing(Q,3,"grevlex");
PX<X>:=PolynomialRing(A);
D:=3+5*e;
A1:=(-2*e-3*e^2)*X^2+(6*e+10*e^2)*X+(1-3*e^2);
A2:=2*e*X^2-(1+3*e);
R3:=2-3*X^2;
S:=R3*(mu*X+nu)^2-D*A1*A2;
s4:=Coefficient(S,4); s3:=Coefficient(S,3);
s2:=Coefficient(S,2); s1:=Coefficient(S,1); s0:=Coefficient(S,0);
H1:=8*s4^2*s1-s3*(4*s4*s2-s3^2);
H0:=64*s4^3*s0-(4*s4*s2-s3^2)^2;

B<ee,mm>:=PolynomialRing(Q,2,"grevlex");
toB:=hom<A -> B | ee,mm,B!0>;
fac:=Factorization(B!toB(Resultant(H1,H0,nu)));
P16:=[f[1]:f in fac|Degree(f[1],ee) eq 16 and
                         Degree(f[1],mm) eq 8][1];
Caff:=Curve(AffineSpace(B),P16);
Cp:=ProjectiveClosure(Caff);

print "CONTACT6_M612_WEIGHTED_R3_P16_INDEPENDENT";
assert IsIrreducible(P16);
assert Genus(Cp) eq 1;
print "P16_SHAPE",<TotalDegree(P16),Degree(P16,ee),Degree(P16,mm),#Terms(P16)>;
print "P16_IRREDUCIBLE_GENUS",true,Genus(Cp);
print "P16_RATIONAL_SINGULAR_PLACE_DEGREES",
      [<q,[Degree(pl):pl in Places(q)]>:q in SingularPoints(Cp)];

// The rational plane singularity (-6/17:0:1) resolves to one closed
// point of degree 2.  Its divisor should give a degree-2 genus-one
// model; twice it gives a genus-one normal quartic in P^3.
q2:=Cp![-6/17,0,1];
assert #Places(q2) eq 1 and Degree(Places(q2)[1]) eq 2;
pl2:=Places(q2)[1];
D2:=Divisor([<pl2,1>]);
V2,rho2:=RiemannRochSpace(D2);
V4,rho4:=RiemannRochSpace(2*D2);
V8,rho8:=RiemannRochSpace(4*D2);
assert <Degree(D2),Dimension(V2),Dimension(V4),Dimension(V8)> eq <2,2,4,8>;
print "D2_RR_DATA",<Degree(D2),Dimension(V2),Dimension(V4),Dimension(V8)>;
f2:=[rho2(V2.i):i in [1..Dimension(V2)]];
f4:=[rho4(V4.i):i in [1..Dimension(V4)]];

// As a cross-check, 2D2 embeds the normalization as an intersection
// of two quadrics, and that model has the same genus and Jacobian.
mons4:=[<i,j>:i,j in [1..4]|i le j];
prods:=[f4[ij[1]]*f4[ij[2]]:ij in mons4];
prodvecs:=[g @@ rho8:g in prods];
M:=Matrix(Q,#prodvecs,Dimension(V8),&cat[Eltseq(v):v in prodvecs]);
ker:=Nullspace(M);
assert Dimension(ker) eq 2;
P3<x1,x2,x3,x4>:=PolynomialRing(Q,4);
qmons:=[];
for ij in mons4 do Append(~qmons,P3.ij[1]*P3.ij[2]); end for;
qs:=[&+[v[i]*qmons[i]:i in [1..#qmons]]:v in Basis(ker)];
C4:=Curve(ProjectiveSpace(Q,3),qs);
assert Dimension(C4) eq 1 and Degree(C4) eq 4 and Genus(C4) eq 1;
G4:=GenusOneModel(C4);
E4:=MinimalModel(Jacobian(G4));
assert aInvariants(E4) eq [0,0,0,3,0];
print "2D2_QUARTIC_CHECK",<Dimension(ker),Degree(C4),Genus(C4)>;
print "2D2_JACOBIAN",E4;

    // Recover the actual degree-2 model attached to D2.  If x,z is a
    // basis of L(D2), then L(2D2) is spanned by x^2,xz,z^2 and one
    // additional function y.  The single relation in L(4D2) is a
    // generalized binary quartic.
xx:=f2[1]; zz:=f2[2];
sym2:=[xx^2,xx*zz,zz^2];
symvecs:=[g @@ rho4:g in sym2];
assert Rank(Matrix(Q,3,4,&cat[Eltseq(v):v in symvecs])) eq 3;
yi:=0;
for i in [1..4] do
    testvecs:=symvecs cat [V4.i];
    if Rank(Matrix(Q,4,4,&cat[Eltseq(v):v in testvecs])) eq 4 then
        yi:=i; break;
    end if;
end for;
assert yi ne 0;
yy:=f4[yi];
deg2prods:=[yy^2,yy*xx^2,yy*xx*zz,yy*zz^2,
           xx^4,xx^3*zz,xx^2*zz^2,xx*zz^3,zz^4];
deg2vecs:=[g @@ rho8:g in deg2prods];
M2:=Matrix(Q,#deg2vecs,Dimension(V8),&cat[Eltseq(v):v in deg2vecs]);
ker2:=Nullspace(M2);
assert Dimension(ker2) eq 1;
rel:=Eltseq(Basis(ker2)[1]);
assert rel[1] ne 0;
rel:=[c/rel[1]:c in rel];
coeffs2:=[rel[2],rel[3],rel[4],-rel[5],-rel[6],
         -rel[7],-rel[8],-rel[9]];
G2:=GenusOneModel(coeffs2);
assert IsIsomorphic(MinimalModel(Jacobian(G2)),E4);
G2min,tmin,bad:=Minimise(G2);
G2red,tred:=Reduce(G2min);
assert Eltseq(G2red) eq [2,0,0,0,-6];
assert aInvariants(MinimalModel(Jacobian(G2red))) eq [0,0,0,3,0];
print "D2_BINARY_QUARTIC_RELATION_DIMENSION",Dimension(ker2);
print "D2_MINIMISE_BAD_PRIMES",bad;
print "D2_MINIMISE_TRANSFORMATION",tmin;
print "D2_REDUCE_TRANSFORMATION",tred;
print "D2_MINIMISED_REDUCED",G2red;
print "D2_REDUCED_JACOBIAN",MinimalModel(Jacobian(G2red));

local_results:=[];
for p in [2,3,5,7,11,13,17] do
    ok,pt:=IsLocallySoluble(G2red,p);
    Append(~local_results,<p,ok>);
end for;
assert local_results[1] eq <2,false>;
assert local_results[2] eq <3,false>;
assert local_results[3] eq <5,true>;
print "D2_LOCAL_RESULTS",local_results;

// Transparent finite congruence certificates for the two obstructions.
N2:=#[<x,z,y>:x,z,y in [0..15] |
      (IsOdd(x) or IsOdd(z)) and
      (y^2-2*x^4+6*z^4) mod 16 eq 0];
N3:=#[<x,z,y>:x,z,y in [0..8] |
      (x mod 3 ne 0 or z mod 3 ne 0) and
      (y^2-2*x^4+6*z^4) mod 9 eq 0];
assert N2 eq 0 and N3 eq 0;
print "PRIMITIVE_CONGRUENCE_COUNTS_MOD16_MOD9",N2,N3;
print "CONCLUSION",
      "normalization(P16)(Q_2)=normalization(P16)(Q_3)=empty, hence no Q-point";

print "CONTACT6_M612_WEIGHTED_R3_P16_INDEPENDENT_DONE";
quit;
