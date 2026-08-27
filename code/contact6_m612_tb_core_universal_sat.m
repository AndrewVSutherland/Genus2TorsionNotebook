//////////////////////////////////////////////////////////////////////
// Universal saturation test for the T_B-halving + [3,3] fiber product.
// Variables are (s,v,r,m,L,U), with K=m^2 and M=L^2.  Saturating before
// finite specialization is the correct way to retain vertical boundary
// limits while removing components which are generically extraneous.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemGB then MemGB:=8;
elif Type(MemGB) eq MonStgElt then MemGB:=StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q:=Rationals();
R<s,v,r,m,L,U>:=PolynomialRing(Q,6,"grevlex");
K:=m^2; M:=L^2;
b:=2*s^2-3;

// Keep a as a univariate elimination variable over Frac(Q[s,v,r,m,L,U]).
FF:=FieldOfFractions(R);
PA<a>:=PolynomialRing(FF);
sF:=FF!s; vF:=FF!v; rF:=FF!r; KF:=FF!K; MF:=FF!M; UF:=FF!U;
bF:=2*sF^2-3;
A3:=a-3+4*sF^2*rF-2*KF;
H1:=sF*((a-3)*rF^2+4*rF-KF*(a+3))-rF*A3;
H2:=8*sF^2*(2*sF^2*rF^2+2*(a-3)*rF+2-KF*(2*sF^2-6))
    -A3^2-32*sF^3*rF;
D:=Coefficient(H1,1);
N:=-Coefficient(H1,0);

function SubNum(poly)
    d:=Degree(poly);
    val:=FF!0;
    for i in [0..d] do val+:=Coefficient(poly,i)*N^i*D^(d-i); end for;
    return R!Numerator(val);
end function;

G0:=SubNum(H2);
while IsDivisibleBy(G0,m) do G0:=ExactQuotient(G0,m); end while;

c1:=2*a+6;
c2:=a^2+2*bF-15;
c3:=2*a*bF+22;
c4:=2*a+bF^2-15;
c5:=2*bF+6;
B3:=c5*MF+3*UF;
Delta3:=4*c4*MF+12*(UF^2+vF^2)-B3^2;
F3:=B3*Delta3+16*vF^3-8*c3*MF-8*UF^3-48*UF*vF^2;
F2:=Delta3^2+64*B3*vF^3-64*c2*MF-192*(UF^2*vF^2+vF^4);
F1:=Delta3*vF^3-4*c1*MF-12*UF*vF^4;
G1:=SubNum(F1); G2:=SubNum(F2); G3:=SubNum(F3);

DR:=R!Numerator(D); NR:=R!Numerator(N);
boundary:=s*v*r*m*L*DR*(U^2-4*v^2)*(NR+(b+2)*DR);
I:=ideal<R|G0,G1,G2,G3>;

print "CONTACT6_M612_TB_CORE_UNIVERSAL_SAT";
print "raw_shapes",[
    <TotalDegree(g),#Terms(g)>:g in [G0,G1,G2,G3]
];
print "raw_dimension"; time print Dimension(I);
print "saturation_begin";
time Isat:=Saturation(I,ideal<R|boundary>);
print "saturation_done";
print "sat_dimension"; time print Dimension(Isat);
print "basis_size",#Basis(Isat);
print "basis_shapes",[<TotalDegree(g),#Terms(g)>:g in Basis(Isat)];
quit;
