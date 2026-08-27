//////////////////////////////////////////////////////////////////////
// Exact model and verification of the rational points found on the
// genus-6 trigonal quotient of the contact-30 C3-root cover.
//
// Run:
//   magma -b code/contact30_g6_quotient_points_verify.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
L<z> := FunctionField(Q);
P<T> := PolynomialRing(L);
M<R> := ext<L | T^2-z*T+(5*z-7)/3>;

// Reconstruct the eps=-1 quotient directly from the order-30 family.
t := (5*R^2-20*R+19)/(R^2-5);
Y := -2*(5*R^2-22*R+25)/(R^2-5);
u := t^3;
s := t^5+t^4+(M!5/2)*t^3+(M!1/2)*t
   - t*(t-M!1/2)*(t+1)*Y;
C := (u^2+1)/(2*u);
denq := u^6+6*u^4*s-2*u^4+15*u^3*s-u*s^3+u^2;
numq := 15*u^5+90*u^4+20*u^3*s-6*u^2*s^2+231*u^3
      + 2*u^2*s-15*u*s^2+90*u^2-20*u*s+15*u-2*s;
q := numq/denq;
A := (s+q)/2;
B := (15-s*q)/2;
assert Eltseq(A)[2] eq 0 and Eltseq(B)[2] eq 0 and Eltseq(C)[2] eq 0;
A := L!Eltseq(A)[1];
B := L!Eltseq(B)[1];
C := L!Eltseq(C)[1];

// Compact integral model.  Put d=5z^2-20z-16 and
// A=PA/d^5, B=PB/d^4, C=PC/d^3.
ZP<zz> := PolynomialRing(Q);
d := 5*zz^2-20*zz-16;
PA := 8461353*zz^10-322506900*zz^9+5471086032*zz^8
    -54340244352*zz^7+349528189824*zz^6-1519391001600*zz^5
    +4514219274240*zz^4-9038973468672*zz^3+11658061873152*zz^2
    -8735943819264*zz+2886221168640;
PB := -445299*zz^8+13749864*zz^7-183283224*zz^6
    +1375090368*zz^5-6337639680*zz^4+18329174016*zz^3
    -32393066496*zz^2+31885688832*zz-13343981568;
PC := 7813*zz^6-187500*zz^5+1849992*zz^4-9583360*zz^3
    +27411072*zz^2-40903680*zz+24776704;
assert A eq L!Evaluate(PA,z)/(L!Evaluate(d,z))^5;
assert B eq L!Evaluate(PB,z)/(L!Evaluate(d,z))^4;
assert C eq L!Evaluate(PC,z)/(L!Evaluate(d,z))^3;

Prho<rho> := PolynomialRing(L);
f := 2*rho^3+(A-3)*rho^2+(B+3)*rho+(C-1);
assert IsIrreducible(f);
F<a> := FunctionField(f);
assert Genus(F) eq 6;
print "GENUS",Genus(F);
print "MODEL d=",d;
print "MODEL PA=",PA;
print "MODEL PB=",PB;
print "MODEL PC=",PC;
print "EQUATION 2*d^5*rho^3+(PA-3*d^5)*rho^2+d*(PB+3*d^4)*rho+d^2*(PC-d^3)=0";

// The projective searches leave these four z-values.  Decompose their
// fibers in the normalization, not merely in the singular plane order.
known_z := [Q!2,Q!5,Q!14/3,Q!32/7];
PQ<X> := PolynomialRing(Q);
for z0 in known_z do
    f0 := PQ![ Evaluate(Coefficient(f,i),z0) : i in [0..3] ];
    dz := Divisor(F!(z-z0));
    dec := [ <Degree(pl),Valuation(dz,pl)> : pl in Support(dz)
             | Valuation(dz,pl) gt 0 ];
    liftpoly := 3*X^2-3*z0*X+5*z0-7;
    delta := Discriminant(liftpoly);
    print "POINT_FIBER", "z",z0,
          "plane_factorization",Factorization(f0),
          "normalization_degree_ramification",dec,
          "lift_discriminant",delta,
          "rational_R_lifts",Roots(liftpoly);
end for;

// No rational point in either projective infinity fiber.
finf := 6250*X^3+8451978*X^2-2217120*X+192200;
assert IsIrreducible(finf);
assert not IsSquare(720); // discriminant of d=5z^2-20z-16
print "Z_INFINITY_FACTORIZATION",Factorization(finf);
print "RHO_INFINITY_Z_EQUATION",d,"discriminant",Discriminant(d);

// Independently verify the rho-projection survivors and all their rational
// z-roots.  Numerators suffice because d has no rational zero.
for rho0 in [Q!-1,Q!0,Q!1] do
    hz := Numerator(Evaluate(f,rho0));
    print "RHO_FIBER",rho0,"factorization",Factorization(hz),
          "rational_z_roots",Roots(hz);
end for;

// The extra quotient point (z,rho)=(32/7,1) has no rational R-lift.  Checking
// its two lifts over Q(sqrt(-3)) shows that u=-1, hence c=0, and both the
// numerator and denominator used to recover q vanish.  It is another
// boundary point, rather than an open source-family point.
K<w> := QuadraticField(-3);
RR := (16+w)/7;
tt := (5*RR^2-20*RR+19)/(RR^2-5);
YY := -2*(5*RR^2-22*RR+25)/(RR^2-5);
uu := tt^3;
ss := tt^5+tt^4+(K!5/2)*tt^3+(K!1/2)*tt
    - tt*(tt-K!1/2)*(tt+1)*YY;
CC := (uu^2+1)/(2*uu);
cc := (uu^2-1)/(2*uu);
dq := uu^6+6*uu^4*ss-2*uu^4+15*uu^3*ss-uu*ss^3+uu^2;
nq := 15*uu^5+90*uu^4+20*uu^3*ss-6*uu^2*ss^2+231*uu^3
    +2*uu^2*ss-15*uu*ss^2+90*uu^2-20*uu*ss+15*uu-2*ss;
print "NEW_LIFT_RAW", "t",tt,"u",uu,"c",cc,"qden",dq,"qnum",nq;
assert uu eq -1 and cc eq 0 and nq eq 0 and dq eq 0;
print "NEW_QUOTIENT_POINT", "z",Q!32/7,"rho",1,
      "R_lifts",[(16+w)/7,(16-w)/7],
      "u",uu,"c_zero",cc eq 0,
      "q_numerator_zero",nq eq 0,"q_denominator_zero",dq eq 0;

print "EXACT_GENUS6_QUOTIENT_POINT_VERIFICATION_OK";
quit;
