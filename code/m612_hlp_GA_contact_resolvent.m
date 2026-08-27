//////////////////////////////////////////////////////////////////////
// Exact cubic-contact support resolvent on the optimized HLP slice
//
//   F_t = F_HLP + t*(2+x-x^2+x^3+x^4+x^5+x^6).
//
// We quotient H -> -H by normalizing its nonzero constant term:
//
//   H^2 = W*(1+r1*x+r2*x^2+r3*x^3)^2.
//
// Both marked HLP contacts lie in this chart.  We retain the seven
// low-degree coefficient equations: in practice their grevlex basis is
// much cheaper than equations obtained by clearing recursive denominators.
// inv*w-1 saturates away from the discarded H(0)=0 chart.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);

Q := Rationals();
K<t> := FunctionField(Q);

// Compute the generic fiber first in grevlex, then use the zero-dimensional
// order change to lex; U remains the last variable.
P<inv,r1,r2,r3,w,k,V,U> := PolynomialRing(K,8,"grevlex");

f0 := K!187392 + 2*t;
f1 := t;
f2 := K!-118767 - t;
f3 := t;
f4 := K!-118767 + t;
f5 := t;
f6 := K!187392 + t;

C0 := w-f0-k*V^3;
C1 := 2*w*r1-f1-3*k*U*V^2;
C2 := w*(r1^2+2*r2)-f2-3*k*(U^2*V+V^2);
C3 := 2*w*(r3+r1*r2)-f3-k*(U^3+6*U*V);
C4 := w*(r2^2+2*r1*r3)-f4-3*k*(U^2+V);
C5 := 2*w*r2*r3-f5-3*k*U;
C6 := w*r3^2-f6-k;

eqs := [C0,C1,C2,C3,C4,C5,C6,inv*w-1];
I := ideal<P | eqs>;

print "M612_HLP_GA_CONTACT_RESOLVENT";
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
    &and[Degree(g,j) eq 0 : j in [1..7]] and Degree(g,8) gt 0];
assert #univ gt 0;
RU := univ[#univ] / LeadingCoefficient(univ[#univ]);
print "resolvent_degree_U",Degree(RU,8);
print "resolvent_term_count",#Terms(RU);
print "resolvent_t_denominator_degree",
      Max([Degree(Denominator(Coefficient(RU,8,i)))
           : i in [0..Degree(RU,8)]]);
print "resolvent_t_numerator_degree",
      Max([Degree(Numerator(Coefficient(RU,8,i)))
           : i in [0..Degree(RU,8)]]);

time fac := Factorization(UnivariatePolynomial(RU));
print "factor_degrees_Qt",[<Degree(a[1]),a[2]> : a in fac];

// Convert the Q(t) polynomial to a primitive polynomial in Q[t,U], then
// specialize exactly at the split HLP fiber.
den := LCM([Denominator(Coefficient(RU,8,i))
            : i in [0..Degree(RU,8)]]);
RUint := den*RU;
Pt<z> := PolynomialRing(Q);
coeff0 := [];
for i in [0..Degree(RUint,8)] do
    ci := Coefficient(RUint,8,i);
    assert Denominator(ci) eq 1;
    Append(~coeff0,Evaluate(Numerator(ci),0));
end for;
RU0 := Pt!coeff0;
RU0sf := RU0 div GCD(RU0,Derivative(RU0));
print "special_fiber_degree_with_multiplicity",Degree(RU0);
print "squarefree_special_fiber_degree",Degree(RU0sf);
print "special_fiber_factor_degrees",
      [<Degree(a[1]),a[2]> : a in Factorization(RU0)];
print "marked_A_U_is_root",Evaluate(RU0sf,-61/8) eq 0;
print "marked_B_U_is_root",Evaluate(RU0sf,0) eq 0;
print "squarefree_special_fiber_resolvent",RU0sf;
quit;
