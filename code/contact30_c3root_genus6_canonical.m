//////////////////////////////////////////////////////////////////////
// Canonical-model probe for the genus-6 trigonal quotient of the
// contact-30 C3-root cover.
//
//   magma -b code/contact30_c3root_genus6_canonical.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(8*10^9);

if not assigned PrintEquation then
    PrintEquation := false;
elif Type(PrintEquation) eq MonStgElt then
    PrintEquation := StringToLower(PrintEquation) in
        { "true", "t", "yes", "y", "1" };
end if;
if not assigned ComputePrincipalTest then
    ComputePrincipalTest := false;
elif Type(ComputePrincipalTest) eq MonStgElt then
    ComputePrincipalTest := StringToLower(ComputePrincipalTest) in
        { "true", "t", "yes", "y", "1" };
end if;

Q := Rationals();
K<z> := FunctionField(Q);
PT<T> := PolynomialRing(K);
M<r> := ext<K | T^2-z*T+(5*z-7)/3>;

RR := r;
t := (5*RR^2-20*RR+19)/(RR^2-5);
Y := -2*(5*RR^2-22*RR+25)/(RR^2-5);
u := t^3;
s := t^5+t^4+(M!5/2)*t^3+(M!1/2)*t
   - t*(t-M!1/2)*(t+1)*Y;
C0 := (u^2+1)/(2*u);
den := u^6+6*u^4*s-2*u^4+15*u^3*s-u*s^3+u^2;
num := 15*u^5+90*u^4+20*u^3*s-6*u^2*s^2+231*u^3
   +2*u^2*s-15*u*s^2+90*u^2-20*u*s+15*u-2*s;
q := num/den;
A := (s+q)/2;
B := (15-s*q)/2;
assert Eltseq(A)[2] eq 0 and Eltseq(B)[2] eq 0
    and Eltseq(C0)[2] eq 0;
aa := K!Eltseq(A)[1];
bb := K!Eltseq(B)[1];
cc := K!Eltseq(C0)[1];
PK<x> := PolynomialRing(K);
f := 2*x^3+(aa-3)*x^2+(bb+3)*x+(cc-1);
assert IsIrreducible(f);
Ffun<a> := FunctionField(f);
assert Genus(Ffun) eq 6;

print "NORMALIZED BOUNDARY PLACE DECOMPOSITION";
boundary_places := [];
zero_fibers := [];
for zv in [Q!2,Q!5,Q!14/3,Q!32/7] do
    Dz := Divisor(Ffun!(z-zv));
    Zz := Numerator(Dz);
    Append(~zero_fibers,Zz);
    pls := Support(Zz);
    print "z",zv,"places",
        [ <Degree(pl),Valuation(Zz,pl),Valuation(Ffun!a,pl),
           Valuation(Ffun!(a-1),pl),Valuation(Ffun!(a+1),pl)>
          : pl in pls ];
    boundary_places cat:= [ pl : pl in pls | Degree(pl) eq 1 ];
end for;
z_poles := Support(Denominator(Divisor(Ffun!z)));
print "z_infinity_places",
    [ <Degree(pl),-Valuation(Divisor(Ffun!z),pl)> : pl in z_poles ];
boundary_places cat:= [ pl : pl in z_poles | Degree(pl) eq 1 ];
print "known_rational_normalized_places",#boundary_places;
assert #boundary_places eq 6;
PA0,PB0,PC0,PD0,PE0,PF0 := Explode(boundary_places);
assert zero_fibers[1] eq Divisor(PA0)+2*Divisor(PB0);
assert zero_fibers[2] eq Divisor(PC0)+2*Divisor(PD0);
assert zero_fibers[3] eq 3*Divisor(PE0);
assert zero_fibers[1]-zero_fibers[2] eq
    Divisor(Ffun!((z-2)/(z-5)));
assert zero_fibers[1]-zero_fibers[3] eq
    Divisor(Ffun!((z-2)/(z-Q!14/3)));
print "exact_rational_point_fiber_relations",
    "A+2B=C+2D=3E";

// At p=13 the shortest possible relation among the three visible
// directions B-E, D-E, F-E is (78,-13,42).  Test it exactly over Q.
G1 := Divisor(PB0)-Divisor(PE0);
G2 := Divisor(PD0)-Divisor(PE0);
G3 := Divisor(PF0)-Divisor(PE0);
if ComputePrincipalTest then
    // This general-function-field principal-divisor test is expensive.
    time short_is_principal := IsPrincipal(78*G1-13*G2+42*G3);
    print "shortest_p13_modular_relation_78_m13_42_principal",
        short_is_principal;
end if;

common_den := LCM([ Denominator(Coefficient(f,i)) : i in [0..3] ]);
R2<Z,X> := PolynomialRing(Q,2);
function RatFunNumeratorToR2(h)
    n := Numerator(h);
    return &+[ Q!Coefficient(n,i)*Z^i : i in [0..Degree(n)] ];
end function;
Faff := &+[ RatFunNumeratorToR2(Coefficient(f,i)*common_den)*X^i
             : i in [0..3] ];
assert IsIrreducible(Faff);

print "CONTACT30 GENUS6 CANONICAL PROBE";
print "affine_bidegree_z_x",Degree(Faff,Z),Degree(Faff,X);
print "affine_total_degree",TotalDegree(Faff),"terms",#Terms(Faff);
if PrintEquation then print "affine_equation",Faff; end if;
print "boundary_fibers";
QX<xx> := PolynomialRing(Q);
for zv in [Q!2,Q!5,Q!14/3] do
    ev := hom<R2 -> QX | zv,xx>;
    print "z",zv,"factorization",
        Factorization(ev(Faff));
end for;

A2 := AffineSpace(Q,2);
Caff := Curve(A2,Faff);
Cp := ProjectiveClosure(Caff);
print "plane_degree",Degree(Cp);
time g,ishyp,canmap := GenusAndCanonicalMap(Cp);
assert g eq 6 and not ishyp;
print "normalized_genus",g,"hyperelliptic",ishyp;
// GenusAndCanonicalMap already returns a map whose codomain is the
// canonical image (unlike CanonicalMap on a smooth input model).
Ccan := Codomain(canmap);
print "canonical_ambient_dimension",Dimension(Ambient(Ccan));
print "canonical_degree",Degree(Ccan);
eqs := DefiningEquations(Ccan);
print "canonical_equation_degrees",Sort([ TotalDegree(h) : h in eqs ]);
print "canonical_equation_count",#eqs;
quadrics := [ h : h in eqs | TotalDegree(h) eq 2 ];
Xscroll := Scheme(Ambient(Ccan),quadrics);
print "canonical_scroll_dimension",Dimension(Xscroll);
print "canonical_scroll_degree",Degree(Xscroll);
print "canonical_scroll_nonsingular",IsNonsingular(Xscroll);

// Recover the unique trigonal map from the canonical model as an
// independent check that the canonical ideal has no extra low-degree
// structure.
time trigmap := CliffordIndexOne(Ccan);
print "clifford_index_one_codomain_dimension",
    Dimension(Ambient(Codomain(trigmap)));

quit;
