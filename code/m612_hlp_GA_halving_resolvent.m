//////////////////////////////////////////////////////////////////////
// Exact 2-halving support resolvent on the optimized HLP slice
//
//   F_t = F_HLP + t*(2+x-x^2+x^3+x^4+x^5+x^6).
//
// Write q0=x^2+c1*x+c0, F=q0*R, and quotient L -> -L by
//
//   L^2 = w*(x+z)^2.
//
// The identity q0*L^2-R=k*(x^2+u1*x+u0)^2 says that the displayed
// quartic is a scalar square.  We retain the four low-degree square
// coefficient equations rather than clear recursive denominators.
// The inverse variable saturates k!=0.  The marked HLP point is
// (c0,c1,z,w)=(1,0,0,183^2).
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);

Q := Rationals();
K<t> := FunctionField(Q);

// Compute in grevlex, then use the zero-dimensional order change to lex;
// z remains the last (trial-resolvent) variable.
P<inv,w,c0,c1,u0,u1,z> := PolynomialRing(K,7,"grevlex");

f0 := K!187392 + 2*t;
f1 := t;
f2 := K!-118767 - t;
f3 := t;
f4 := K!-118767 + t;
f5 := t;
f6 := K!187392 + t;

// Quotient coefficients R=r0+...+r4*x^4, recovered from the top down.
r4 := f6;
r3 := f5-c1*r4;
r2 := f4-c1*r3-c0*r4;
r1 := f3-c1*r2-c0*r3;
r0 := f2-c1*r1-c0*r2;

// The last two coefficients impose q0 | F.
D1 := f1-c1*r0-c0*r1;
D0 := f0-c0*r0;

// S=q0*w*(x+z)^2-R.
s0 := w*c0*z^2-r0;
s1 := w*(2*c0*z+c1*z^2)-r1;
s2 := w*(c0+2*c1*z+z^2)-r2;
s3 := w*(c1+2*z)-r3;
s4 := w-r4;

k  := s4;
C3 := s3-2*k*u1;
C2 := s2-k*(u1^2+2*u0);
C1 := s1-2*k*u1*u0;
C0 := s0-k*u0^2;
eqs := [D1,D0,C3,C2,C1,C0,inv*k-1];

I := ideal<P | eqs>;

seed := [ K!(1/(-153903)), K!(183^2), K!1, K!0,
          K!(-32/29), K!0, K!0 ];
seedvals := [Evaluate(e,seed) : e in eqs];
assert &and[Evaluate(Numerator(e),0) eq 0 : e in seedvals];

print "M612_HLP_GA_HALVING_RESOLVENT";
print "reduced_equation_shapes",
      [<TotalDegree(e),#Terms(e)> : e in eqs];

time GB0 := GroebnerBasis(I);
print "grevlex_basis_size",#GB0;
dim,parameters := Dimension(I);
print "ideal_dimension",dim,"parameters",parameters;
print "ideal_degree",Dimension(quo<P|I>);


time IL := ChangeOrder(I,"lex");
time GB := GroebnerBasis(IL);
print "lex_basis_size",#GB;

univ := [g : g in GB |
    &and[Degree(g,j) eq 0 : j in [1..6]] and Degree(g,7) gt 0];
if #univ eq 0 then
    print "z_is_not_a_finite_resolvent";
    quit;
end if;
RZ := univ[#univ] / LeadingCoefficient(univ[#univ]);
print "resolvent_degree_z",Degree(RZ,7);
print "resolvent_term_count",#Terms(RZ);
print "resolvent_t_denominator_degree",
      Max([Degree(Denominator(Coefficient(RZ,7,i)))
           : i in [0..Degree(RZ,7)]]);
print "resolvent_t_numerator_degree",
      Max([Degree(Numerator(Coefficient(RZ,7,i)))
           : i in [0..Degree(RZ,7)]]);

time fac := Factorization(UnivariatePolynomial(RZ));
print "factor_degrees_Qt",[<Degree(a[1]),a[2]> : a in fac];

den := LCM([Denominator(Coefficient(RZ,7,i))
            : i in [0..Degree(RZ,7)]]);
RZint := den*RZ;
Pt<zz> := PolynomialRing(Q);
coeff0 := [];
for i in [0..Degree(RZint,7)] do
    ci := Coefficient(RZint,7,i);
    assert Denominator(ci) eq 1;
    Append(~coeff0,Evaluate(Numerator(ci),0));
end for;
RZ0 := Pt!coeff0;
RZ0sf := RZ0 div GCD(RZ0,Derivative(RZ0));
print "special_fiber_degree_with_multiplicity",Degree(RZ0);
print "squarefree_special_fiber_degree",Degree(RZ0sf);
print "special_fiber_factor_degrees",
      [<Degree(a[1]),a[2]> : a in Factorization(RZ0)];
print "marked_G_z_is_root",Evaluate(RZ0sf,0) eq 0;
print "squarefree_special_fiber_resolvent",RZ0sf;
quit;
