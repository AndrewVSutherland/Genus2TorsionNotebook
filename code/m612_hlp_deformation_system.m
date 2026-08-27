//////////////////////////////////////////////////////////////////////
// Simultaneous two-contact plus 2-halving deformation system at the
// explicit split HLP [6,12] seed.
//
// This is a derivation/geometry file, not a search job.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
R<f0,f1,f2,f3,f4,f5,f6,
  a0,a1,ha0,ha1,ha2,ha3,ka,
  b0,b1,hb0,hb1,hb2,hb3,kb,
  c0,c1,u0,u1,l0,l1,kg> := PolynomialRing(Q,28);
P<x> := PolynomialRing(R);

F  := &+[ [f0,f1,f2,f3,f4,f5,f6][i+1]*x^i : i in [0..6] ];
qA := x^2+a1*x+a0;
hA := ha3*x^3+ha2*x^2+ha1*x+ha0;
qB := x^2+b1*x+b0;
hB := hb3*x^3+hb2*x^2+hb1*x+hb0;

// The rational 2-class is [q0,0].  Writing ell=q0*(l1*x+l0) forces
// q0|F.  The identity ell^2-F=kg*uG^2*q0 makes [uG,ell mod uG]
// a half of [q0,0].
q0 := x^2+c1*x+c0;
uG := x^2+u1*x+u0;
ell := q0*(l1*x+l0);

EA := [Coefficient(hA^2-F-ka*qA^3,i) : i in [0..6]];
EB := [Coefficient(hB^2-F-kb*qB^3,i) : i in [0..6]];
EG := [Coefficient(ell^2-F-kg*uG^2*q0,i) : i in [0..6]];
eqs := EA cat EB cat EG;

seed := [
    187392,0,-118767,0,-118767,0,187392,
    1,-61/8,-3904/9,61/3,-61/3,3904/9,62464/81,
    -13/48,0,2623/6,0,-183,0,-187392,
    1,0,-32/29,0,0,183,-153903
];
assert &and[Evaluate(e,seed) eq 0 : e in eqs];

Jac := Matrix(Q,[[Evaluate(Derivative(e,j),seed) : j in [1..28]]
                 : e in eqs]);
print "M612_HLP_DEFORMATION_SYSTEM";
print "variables",28,"equations",#eqs,
      "Jacobian_rank",Rank(Jac),"tangent_dimension",28-Rank(Jac);

// Exact algebraic curve slice with tangent dF=1+x:
// F_t=F_seed+t*(1+x).  These six linear equations leave f0/f1 linked.
slice_eqs := [
    f0-f1-187392,
    f2+118767,
    f3,
    f4+118767,
    f5,
    f6-187392
];
slice_all := eqs cat slice_eqs;
assert &and[Evaluate(e,seed) eq 0 : e in slice_all];
SliceJac := Matrix(Q,
    [[Evaluate(Derivative(e,j),seed) : j in [1..28]] : e in slice_all]);
print "transverse_slice equations",#slice_all,
      "Jacobian_rank",Rank(SliceJac),
      "local_dimension",28-Rank(SliceJac);
print "slice curve polynomial F_t = F_seed+t*(1+x), t=f1";
print "The rational seed is t=0; finding another rational point is open.";
quit;
