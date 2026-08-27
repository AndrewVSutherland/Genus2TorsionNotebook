//////////////////////////////////////////////////////////////////////
// Exact R3-halving cover on the weighted endpoint b=0, a=1/e.
//
// For the dual curve y^2=Delta*R1*R2*R3, halving the class R3 is
// equivalent on the square-quartic chart to
//
//   R3*(M*x+N)^2-Delta*R1*R2 = scalar*(quadratic)^2.
//
// Near e=0 the balanced scaling is M=mu/e^2, N=nu/e^2.  Multiplying
// the quartic by e^4 gives the integral quartic S below.  The two standard
// square-quartic covariants are therefore an exact cover over Q(e).
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Z:=Integers();
Q:=Rationals();
p:=5;
k:=GF(p);

A<e,mu,nu>:=PolynomialRing(Q,3,"grevlex");
PX<X>:=PolynomialRing(A);

// Cleared endpoint factors: D=e*Delta, A1=e^2*R1, A2=e*R2.
D:=3+5*e;
A1:=(-2*e-3*e^2)*X^2+(6*e+10*e^2)*X+(1-3*e^2);
A2:=2*e*X^2-(1+3*e);
R3:=2-3*X^2;

S:=R3*(mu*X+nu)^2-D*A1*A2;
s4:=Coefficient(S,4);
s3:=Coefficient(S,3);
s2:=Coefficient(S,2);
s1:=Coefficient(S,1);
s0:=Coefficient(S,0);
H1:=8*s4^2*s1-s3*(4*s4*s2-s3^2);
H0:=64*s4^3*s0-(4*s4*s2-s3^2)^2;

// Verify S=e^4*(R3*(M*x+N)^2-Delta*R1*R2) exactly.
K:=FieldOfFractions(A);
KX<x>:=PolynomialRing(K);
aa:=1/(K!e);
Delta:=3*aa+5;
R1:=(-2*aa-3)*x^2+(6*aa+10)*x+(aa^2-3);
R2:=2*x^2-(aa+3);
R3K:=2-3*x^2;
M:=K!mu/(K!e)^2;
N:=K!nu/(K!e)^2;
Sorig:=R3K*(M*x+N)^2-Delta*R1*R2;
assert KX!S eq (K!e)^4*Sorig;

print "CONTACT6_M612_WEIGHTED_E9_R3";
print "ENDPOINT", "a=1/e b=0";
print "HALVING_SCALING", "M=mu/e^2 N=nu/e^2";
print "INTEGRAL_QUARTIC_S",S;
print "SQUARE_QUARTIC_H1",H1;
print "SQUARE_QUARTIC_H0",H0;
print "H1_SHAPE",<TotalDegree(H1),Degree(H1,e),Degree(H1,mu),
                       Degree(H1,nu),#Terms(H1)>;
print "H0_SHAPE",<TotalDegree(H0),Degree(H0,e),Degree(H0,mu),
                       Degree(H0,nu),#Terms(H0)>;
print "OPEN_CONDITION", "s4 != 0",s4;

// The boundary fiber is small enough to identify exactly.
B<m,n>:=PolynomialRing(Q,2,"grevlex");
h:=hom<A -> B | B!0,m,n>;
H1e0:=h(H1);
H0e0:=h(H0);
s4e0:=h(s4);
assert H1e0 eq 144*m^5*n;
assert H0e0 eq -576*m^6*(m^2+6*n^2+9);
assert s4e0 eq -3*m^2;
print "E_ZERO_H1_FACTORED",Factorization(H1e0);
print "E_ZERO_H0_FACTORED",Factorization(H0e0);

function EvalZ(g,pt)
    z:=Evaluate(g,[Q!a:a in pt]);
    assert Denominator(z) eq 1;
    return Z!z;
end function;

// On the quartic-open chart mu != 0, there are two smooth F_5 points.
r3_points:=[];
r3_ranks:=[];
for m0 in [1..p-1] do for n0 in [0..p-1] do
    pt:=[0,m0,n0];
    if (EvalZ(s4,pt) mod p) eq 0 then continue; end if;
    if (EvalZ(H1,pt) mod p) ne 0 or (EvalZ(H0,pt) mod p) ne 0 then
        continue;
    end if;
    J:=Matrix(k,2,2,[
        k!EvalZ(Derivative(H1,mu),pt),
        k!EvalZ(Derivative(H1,nu),pt),
        k!EvalZ(Derivative(H0,mu),pt),
        k!EvalZ(Derivative(H0,nu),pt)
    ]);
    Append(~r3_points,[m0,n0]);
    Append(~r3_ranks,Rank(J));
    print "R3_POINT_MOD5",[m0,n0],"MU_NU_JAC_RANK",Rank(J);
end for; end for;

assert r3_points eq [[1,0],[4,0]];
assert r3_ranks eq [2,2];

// Since the E9 core has e=25E, e is already zero modulo 25.  The two
// smooth R3 branches lift uniquely to the following residues modulo 25.
r3_mod25:=[];
for m0 in [0..p^2-1] do for n0 in [0..p^2-1] do
    if m0 mod p eq 0 then continue; end if;
    pt:=[0,m0,n0];
    if (EvalZ(H1,pt) mod p^2) eq 0 and
       (EvalZ(H0,pt) mod p^2) eq 0 and
       (EvalZ(s4,pt) mod p) ne 0 then
        Append(~r3_mod25,[m0,n0]);
    end if;
end for; end for;

assert r3_mod25 eq [[4,0],[21,0]];
print "R3_POINTS_MOD25_mu_nu",r3_mod25;
print "R3_LOCAL_VERDICT",
      "two smooth R3-halving branches exist for every e in 5*Z_5";
print "E9_COMBINED_VERDICT",
      "each of the four next-layer E9 core points is locally compatible with R3";
print "GLOBAL_WARNING",
      "rational [6,12] still requires rational solutions of H1=H0=0, not only Q_5 points";
print "CONTACT6_M612_WEIGHTED_E9_R3_DONE";
quit;
