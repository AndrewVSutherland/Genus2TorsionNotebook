//////////////////////////////////////////////////////////////////////
// Local deformation system through the explicit split HLP [6,12] seed.
//
// The order-4 layer is imposed by
//   ell^2-f = lambda*G^2*Q,  ell=Q*(M*x+N),
// so H=[G,ell mod G] doubles to the rational 2-class [Q,0].
// Two independent order-3 directions are imposed by
//   h_i^2-f = k_i*q_i^3.
// This script verifies the seed and computes the Zariski tangent dimension,
// both before and after a four-condition coordinate/scaling gauge.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
QQ:=Rationals();
R<q1,q0,g1,g0,M,N,lam,
  u1,v1,a13,a12,a11,a10,k1,
  u2,v2,b13,b12,b11,b10,k2>:=PolynomialRing(QQ,21,"grevlex");
PX<x>:=PolynomialRing(R);

Q2:=x^2+q1*x+q0;
G2:=x^2+g1*x+g0;
ell:=Q2*(M*x+N);
f:=ell^2-lam*G2^2*Q2;
qa:=x^2+u1*x+v1;
ha:=a13*x^3+a12*x^2+a11*x+a10;
qb:=x^2+u2*x+v2;
hb:=b13*x^3+b12*x^2+b11*x+b10;
Ea:=ha^2-f-k1*qa^3;
Eb:=hb^2-f-k2*qb^3;
eqs:=[Coefficient(Ea,i):i in [0..6]] cat
     [Coefficient(Eb,i):i in [0..6]];

seed:=[QQ!0,QQ!1,QQ!0,QQ!-32/29,QQ!183,QQ!0,QQ!-153903,
       QQ!-61/8,QQ!1,QQ!3904/9,QQ!-61/3,QQ!61/3,QQ!-3904/9,
       QQ!62464/81,
       QQ!0,QQ!-13/48,QQ!0,QQ!-183,QQ!0,QQ!2623/6,QQ!-187392];
badseed:=[<i,Evaluate(eqs[i],seed)>:i in [1..#eqs]
          | Evaluate(eqs[i],seed) ne 0];
print "SEED_RESIDUALS",badseed;
assert #badseed eq 0;

J:=Matrix(QQ,#eqs,21,
          [Evaluate(Derivative(e,j),seed):e in eqs,j in [1..21]]);

// Fix Q=x^2+1, the x-coefficient of G, and M=183.  These are four
// independent coordinate/scaling normalizations at the seed.
fixed:={1,2,3,5};
freecols:=[j:j in [1..21]|j notin fixed];
Jg:=Submatrix(J,1,freecols[1],Nrows(J),1);
// The columns are not contiguous, so build the restricted matrix directly.
Jg:=Matrix(QQ,#eqs,#freecols,[J[i,j]:i in [1..#eqs],j in freecols]);

print "M612_HLP_DEFORMATION_TANGENT";
print "ORDER4_IDENTITY","ell^2-f=lambda*G^2*Q";
PQQ<X>:=PolynomialRing(QQ);
function EvalPolyAtSeed(poly,pt)
    return PQQ!(&+[Evaluate(Coefficient(poly,i),pt)*X^i
                   :i in [0..Degree(poly)]]);
end function;
print "F_SEED",EvalPolyAtSeed(f,seed);
print "CONTACT_A",EvalPolyAtSeed(ha,seed),EvalPolyAtSeed(qa,seed),
      Evaluate(k1,seed);
print "CONTACT_B",EvalPolyAtSeed(hb,seed),EvalPolyAtSeed(qb,seed),
      Evaluate(k2,seed);
print "JACOBIAN_RANK",Rank(J),"AMBIENT",21,"TANGENT_DIM",21-Rank(J);
print "GAUGE_FIXED_VARIABLE_INDICES",Sort(Setseq(fixed));
print "GAUGE_FREE_VARIABLE_INDICES",freecols;
print "GAUGE_JACOBIAN_RANK",Rank(Jg),"AMBIENT",#freecols,
      "TANGENT_DIM",#freecols-Rank(Jg);
rightker:=Nullspace(Transpose(Jg));
print "GAUGE_NULLSPACE_DIM",Dimension(rightker);
for w in Basis(rightker) do
    print "TANGENT_BASIS",w;
end for;
quit;
