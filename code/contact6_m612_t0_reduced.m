//////////////////////////////////////////////////////////////////////
// Reduced exact halving cover for T0=[x,0] in the contact-6 family.
//
// Starting from
//   B*C - x*(m*x+n)^2 = c4*(x^2+u*x+w)^2,
// the constant coefficient gives a+3=(b+3)w^2.  Put s=b+3 and
// W=m^2.  Eliminating u and n gives one quartic F(s,w,W)=0, followed
// by the square lift W=m^2.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
R<s,w,W> := PolynomialRing(Q,3);

b := s-3;
a := s*w^2-3;
c4 := 2*s;
c3 := b^2+2*a-15;
c2 := 2*a*b+22;
c1 := a^2+2*b-15;

// u=(c3-W)/(4s), n^2=N, and 2mn=R8/(8s).
N := c1-w*c3+w*W;
R8 := 8*s*(c2-4*s*w)-(c3-W)^2;
F := R8^2-256*s^2*W*N;

print "CONTACT6_M612_T0_REDUCED";
print "a",a,"b",b;
print "c3",c3;
print "c2",c2;
print "c1",c1;
print "N",N;
print "R8",R8;
print "F_total_degree",TotalDegree(F),"F_degree_W",Degree(F,3),
      "F_terms",#Terms(F);
print "F_factorization",Factorization(F);
print "F",F;

Q2<ss,WW> := PolynomialRing(Q,2);
for wfix in [-2,-1,1,2] do
    spec := Evaluate(F,[Q2.1,Q!wfix,Q2.2]);
    print "SPECIAL_W",wfix,"degree",TotalDegree(spec),
          "terms",#Terms(spec),"factorization",Factorization(spec);
end for;

// Symbolic reconstruction check after adjoining m with m^2=W.
S<s0,w0,m> := PolynomialRing(Q,3);
hom := hom<R->S | s0,w0,m^2>;
c3s := hom(c3);
c2s := hom(c2);
c1s := hom(c1);
Ns := hom(N);
R8s := hom(R8);
us := (c3s-m^2)/(4*s0);
ns := R8s/(16*s0*m);
print "RECONSTRUCTION";
print "F_AFTER_W_M2_degree_m",Degree(hom(F),3),
      "factorization",Factorization(hom(F));
print "u",us;
print "n",ns;
print "n2_minus_N_identity",
      Numerator(ns^2-Ns-hom(F)/(256*s0^2*m^2)) eq 0;

quit;
