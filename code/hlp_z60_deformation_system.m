//////////////////////////////////////////////////////////////////////
// Simultaneous marked order-5/order-3/order-4 deformation system at
// the exact split HLP cyclic-[60] seed.
//
// The order-5 norm block uses the local gauge B5(0)=1.  This leaves
// 32 variables and 25 coefficient equations in total.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemGB then MemGB:=2; end if;
if Type(MemGB) eq MonStgElt then MemGB:=StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q:=Rationals();
R<f0,f1,f2,f3,f4,f5,f6,
  q50,q51,a0,a1,a2,a3,a4,a5,b1,b2,k5,
  q30,q31,h0,h1,h2,h3,k3,
  q20,q21,u40,u41,l0,l1,k4>:=PolynomialRing(Q,32);
P<x>:=PolynomialRing(R);

F:=&+[[f0,f1,f2,f3,f4,f5,f6][i+1]*x^i:i in [0..6]];

Q5:=x^2+q51*x+q50;
A5:=a5*x^5+a4*x^4+a3*x^3+a2*x^2+a1*x+a0;
B5:=1+b1*x+b2*x^2;
E5:=[Coefficient(A5^2-B5^2*F-k5*Q5^5,i):i in [0..10]];

Q3:=x^2+q31*x+q30;
H3:=h3*x^3+h2*x^2+h1*x+h0;
E3:=[Coefficient(H3^2-F-k3*Q3^3,i):i in [0..6]];

Q2:=x^2+q21*x+q20;
U4:=x^2+u41*x+u40;
L4:=l1*x+l0;
Ell:=Q2*L4;
E4:=[Coefficient(Ell^2-F-k4*U4^2*Q2,i):i in [0..6]];
eqs:=E5 cat E3 cat E4;

seed:=[
  94277468160,0,-22332312000,0,1761500625,0,-46250000,
  -4608/395,0,
  0,66600,0,-1558625/144,0,48752125/110592,
  0,-125/1728,1778923230671875/4076863488,
  -316/25,0,-29600,0,2775,0,46250000,
  -1728/125,0,-1506/125,0,2775,0,46250000
];
assert #seed eq 32 and #eqs eq 25;
assert &and[Evaluate(e,seed) eq 0:e in eqs];

Jac:=Matrix(Q,[[Evaluate(Derivative(e,j),seed):j in [1..32]]:e in eqs]);
print "HLP_Z60_SIMULTANEOUS_DEFORMATION_SYSTEM";
print "variables",32,"equations",#eqs,
      "Jacobian_rank",Rank(Jac),"tangent_dimension",32-Rank(Jac);
print "block_ranks",
      Rank(Submatrix(Jac,1,8,11,11)),
      Rank(Submatrix(Jac,12,19,7,7)),
      Rank(Submatrix(Jac,19,26,7,7));

// Exact curve slice F_t=F_seed+t*(1+x), with t=f1.
slice_eqs:=[
  f0-f1-94277468160,
  f2+22332312000,
  f3,
  f4-1761500625,
  f5,
  f6+46250000
];
all_eqs:=eqs cat slice_eqs;
assert &and[Evaluate(e,seed) eq 0:e in all_eqs];
SliceJac:=Matrix(Q,
  [[Evaluate(Derivative(e,j),seed):j in [1..32]]:e in all_eqs]);
print "transverse_slice equations",#all_eqs,
      "Jacobian_rank",Rank(SliceJac),
      "local_dimension",32-Rank(SliceJac);
print "slice F_t=F_seed+t*(1+x), t=f1";
print "HLP_Z60_SIMULTANEOUS_DEFORMATION_SYSTEM_DONE";
quit;
