//////////////////////////////////////////////////////////////////////
// Exact Frobenius polynomial of the genus-6 trigonal quotient of the
// contact-30 C3-root cover.
//
// Examples:
//   magma -b Prime:=13 ComputeGalois:=true \
//       code/contact30_c3root_genus6_frobenius.m
//   magma -b Prime:=17 code/contact30_c3root_genus6_frobenius.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(3*10^9);

if not assigned Prime then
    Prime := 13;
elif Type(Prime) eq MonStgElt then
    Prime := StringToInteger(Prime);
end if;
if not assigned ComputeGalois then
    ComputeGalois := false;
elif Type(ComputeGalois) eq MonStgElt then
    ComputeGalois := StringToLower(ComputeGalois) in
        { "true", "t", "yes", "y", "1" };
end if;

Q := Rationals();
L<z> := FunctionField(Q);
PT<T> := PolynomialRing(L);
M<r> := ext<L | T^2-z*T+(5*z-7)/3>;

// eps=-1 quotient model.  The coefficients below are checked to lie in
// the invariant field Q(z), exactly as in the geometry script.
RR := r;
t := (5*RR^2-20*RR+19)/(RR^2-5);
Y := -2*(5*RR^2-22*RR+25)/(RR^2-5);
u := t^3;
s := t^5+t^4+(M!5/2)*t^3+(M!1/2)*t
   - t*(t-M!1/2)*(t+1)*Y;
C := (u^2+1)/(2*u);
den := u^6+6*u^4*s-2*u^4+15*u^3*s-u*s^3+u^2;
num := 15*u^5+90*u^4+20*u^3*s-6*u^2*s^2+231*u^3
   +2*u^2*s-15*u*s^2+90*u^2-20*u*s+15*u-2*s;
q := num/den;
A := (s+q)/2;
B := (15-s*q)/2;
assert Eltseq(A)[2] eq 0 and Eltseq(B)[2] eq 0 and Eltseq(C)[2] eq 0;
aa := L!Eltseq(A)[1];
bb := L!Eltseq(B)[1];
cc := L!Eltseq(C)[1];
PL<x> := PolynomialRing(L);
f := 2*x^3+(aa-3)*x^2+(bb+3)*x+(cc-1);
assert IsIrreducible(f);

p := Prime;
assert IsPrime(p) and p notin {2,3,5};
k := GF(p);
K<zz> := FunctionField(k);
PK<xx> := PolynomialRing(K);

function ReduceFunction(c)
    return Evaluate(ChangeRing(Numerator(c),k),zz)
        / Evaluate(ChangeRing(Denominator(c),k),zz);
end function;

fp := &+[ ReduceFunction(Coefficient(f,i))*xx^i : i in [0..3] ];
if not IsIrreducible(fp) then
    print "Prime",p,"reduced_cubic_irreducible",false;
    quit;
end if;
Fp<a> := FunctionField(fp);
gp := Genus(Fp);
print "Prime",p,"reduced_genus",gp;
if gp ne 6 then
    print "not_good_reduction_in_this_model";
    quit;
end if;

time zeta := ZetaFunction(Fp);
Lpoly := Numerator(zeta);
assert Degree(Lpoly) eq 12;
QX<X> := PolynomialRing(Q);
LQ := QX!Lpoly;
charpoly := ReciprocalPolynomial(LQ);
fac := Factorization(charpoly);
print "L_polynomial",Lpoly;
print "Frobenius_characteristic_polynomial",charpoly;
print "factor_degrees_Q",[ <Degree(fe[1]),fe[2]> : fe in fac ];
print "irreducible_Q",#fac eq 1 and fac[1][2] eq 1;
print "ordinary",Coefficient(Lpoly,6) mod p ne 0;
print "curve_points_Fp",p+1+Integers()!Coefficient(Lpoly,1);
print "jacobian_points_Fp",Integers()!Evaluate(Lpoly,1);

if ComputeGalois then
    time G,roots := GaloisGroup(charpoly);
    print "Galois_group_order",Order(G);
    try
        print "Galois_transitive_id",TransitiveGroupIdentification(G);
    catch e
        print "Galois_transitive_id_unavailable";
    end try;
    print "expected_full_hyperoctahedral_order",2^6*Factorial(6);
    for n in [2..12] do
        Pn := WeilPolynomialOverFieldExtension(charpoly,n);
        fn := Factorization(Pn);
        print "extension_degree",n,
              "factor_degrees",[ <Degree(fe[1]),fe[2]> : fe in fn ];
    end for;
end if;

quit;
