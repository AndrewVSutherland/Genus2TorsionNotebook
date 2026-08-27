//////////////////////////////////////////////////////////////////////
// Exact good-characteristic geometry certificates for the optimized
// HLP G_A slice.  Usage examples:
//
//   magma -b block:=contact p:=101 code/m612_hlp_GA_modular_geometry.m
//   magma -b block:=halving p:=101 code/m612_hlp_GA_modular_geometry.m
//
// An irreducible full-degree reduction certifies generic irreducibility
// in characteristic zero (after the displayed degree check).  Its
// function-field genus is also a lower-bound obstruction to a rational
// or elliptic characteristic-zero component through a good reduction.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);

if not assigned block then block := "contact"; end if;
if not assigned p then p := 101;
elif Type(p) eq MonStgElt then p := StringToInteger(p);
end if;
if not assigned do_genus then do_genus := false;
elif Type(do_genus) eq MonStgElt then
    do_genus := do_genus eq "true" or do_genus eq "1";
end if;

Fp := GF(p);
K<t> := FunctionField(Fp);

Fpt<tt> := PolynomialRing(Fp);
Fpx<xx> := PolynomialRing(Fpt);
Fcheck := 187392*xx^6-118767*xx^4-118767*xx^2+187392
          +tt*(2+xx-xx^2+xx^3+xx^4+xx^5+xx^6);
Dcheck := Discriminant(Fcheck);
Gcheck := 2+xx-xx^2+xx^3+xx^4+xx^5+xx^6;
print "good_pencil_reduction",
      Degree(Dcheck) eq 10 and GCD(Dcheck,Derivative(Dcheck)) eq 1
      and Discriminant(Gcheck) ne 0;

f0 := K!187392 + 2*t;
f1 := t;
f2 := K!-118767 - t;
f3 := t;
f4 := K!-118767 + t;
f5 := t;
f6 := K!187392 + t;

if block eq "contact" then
    P<inv,r1,r2,r3,w,k,V,U> := PolynomialRing(K,8,"grevlex");
    C0 := w-f0-k*V^3;
    C1 := 2*w*r1-f1-3*k*U*V^2;
    C2 := w*(r1^2+2*r2)-f2-3*k*(U^2*V+V^2);
    C3 := 2*w*(r3+r1*r2)-f3-k*(U^3+6*U*V);
    C4 := w*(r2^2+2*r1*r3)-f4-3*k*(U^2+V);
    C5 := 2*w*r2*r3-f5-3*k*U;
    C6 := w*r3^2-f6-k;
    I := ideal<P | C0,C1,C2,C3,C4,C5,C6,inv*w-1>;
    expected := 40;
    last := 8;
else
    assert block eq "halving";
    P<inv,w,c0,c1,u0,u1,z> := PolynomialRing(K,7,"grevlex");
    r4 := f6;
    r3 := f5-c1*r4;
    r2 := f4-c1*r3-c0*r4;
    r1 := f3-c1*r2-c0*r3;
    r0 := f2-c1*r1-c0*r2;
    D1 := f1-c1*r0-c0*r1;
    D0 := f0-c0*r0;
    s0 := w*c0*z^2-r0;
    s1 := w*(2*c0*z+c1*z^2)-r1;
    s2 := w*(c0+2*c1*z+z^2)-r2;
    s3 := w*(c1+2*z)-r3;
    kk := w-r4;
    C3 := s3-2*kk*u1;
    C2 := s2-kk*(u1^2+2*u0);
    C1 := s1-2*kk*u1*u0;
    C0 := s0-kk*u0^2;
    I := ideal<P | D1,D0,C3,C2,C1,C0,inv*kk-1>;
    expected := 120;
    last := 7;
end if;

print "M612_HLP_GA_MODULAR_GEOMETRY",block,"p",p;
time G0 := GroebnerBasis(I);
print "grevlex_basis_size",#G0;
dim,parameters := Dimension(I);
deg := Dimension(quo<P|I>);
print "dimension",dim,"parameters",parameters,
      "degree",deg,"expected",expected;
time IL := ChangeOrder(I,"lex");
time G := GroebnerBasis(IL);
print "lex_basis_size",#G;

univ := [g : g in G |
    &and[Degree(g,j) eq 0 : j in [1..last-1]] and Degree(g,last) gt 0];
assert #univ gt 0;
R := UnivariatePolynomial(univ[#univ]);
R /:= LeadingCoefficient(R);
fac := Factorization(R);
print "resolvent_degree",Degree(R);
print "factor_degrees",[<Degree(a[1]),a[2]> : a in fac];
print "full_degree_coordinate",Degree(R) eq expected;
if do_genus and #fac eq 1 and fac[1][2] eq 1 then
    print "function_field_genus_begin";
    time L := FunctionField(R : Check := true);
    time g := Genus(L);
    print "function_field_genus",g;
end if;
quit;
