//////////////////////////////////////////////////////////////////////
// Exact, non-eliminated degree-16 halving cover for the marked
// order-22 class.  This intentionally performs no Groebner basis or
// resultant computation.
//
// Run one of the eight charts, for example:
//   magma -b family:="Flynn" source_sign:=1 half_sign:=1 \
//       code/order44_from_order22_halving_cover.m
//
// family is "Flynn" or "DaowsudSchmidt"; source_sign chooses the
// order-22 parameter branch; half_sign chooses ell(0)=+/-1.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();

if not assigned family then family := "Flynn"; end if;
if not assigned source_sign then source_sign := 1; end if;
if not assigned half_sign then half_sign := 1; end if;
if Type(source_sign) eq MonStgElt then source_sign := StringToInteger(source_sign); end if;
if Type(half_sign) eq MonStgElt then half_sign := StringToInteger(half_sign); end if;
require family in {"Flynn","DaowsudSchmidt"}: "unknown family";
require source_sign in {-1,1}: "source_sign must be +/-1";
require half_sign in {-1,1}: "half_sign must be +/-1";

Ks<s> := FunctionField(Q);
PX<X> := PolynomialRing(Ks);

num := -s^2*(s^2+1)*(s^4-s^2+1) + 2*source_sign*s^5;
den := (s^2-1)^2;

if family eq "Flynn" then
    par := num/den;
    r := s^2;
    f := X^6 + 2*X^5 + (2*par+3)*X^4 + 2*X^3
         + (par^2+1)*X^2 + 2*par*(1-par)*X + par^2;
else
    par := num/(4*den);
    r := 1+s^2;
    f := X^6 - 4*X^5 + 8*(1+par)*X^4 - (10+32*par)*X^3
         + 8*(1+6*par+2*par^2)*X^2
         - 4*(1+6*par+16*par^2)*X + 64*par^2+1;
end if;

// g(U)=U^6*f(r+1/U), retained in the same polynomial variable X.
g := PX!(&+[Coefficient(f,i)*(r*X+1)^i*X^(6-i) : i in [0..6]]);
assert Degree(g) eq 5;
assert Coefficient(g,0) eq 1;
assert Evaluate(f,r) eq 0;

Rmn<m,n> := PolynomialRing(Ks,2);
gi := [Rmn!Coefficient(g,i) : i in [0..5]];
k := Rmn!half_sign;

// If u=X^2+aX+b and ell=mX^2+nX+k, then
//       g-ell^2 = g5*X*u^2.
// The X^5 and X^4 coefficients solve for a,b; K2=K1=0 are
// precisely the remaining X^2 and X coefficients.
a := (gi[5]-m^2)/(2*gi[6]);
b := (gi[4]-2*m*n-gi[6]*a^2)/(2*gi[6]);
K2 := gi[3]-n^2-2*m*k-2*gi[6]*a*b;
K1 := gi[2]-2*n*k-gi[6]*b^2;

function ClearFunctionDenominators(h)
    coeffs := Coefficients(h);
    denh := LCM([Denominator(c) : c in coeffs]);
    return Rmn!(denh*h),denh;
end function;

K2num,K2den := ClearFunctionDenominators(K2);
K1num,K1den := ClearFunctionDenominators(K1);

print "ORDER44_FROM_ORDER22_HALVING_COVER";
print "family",family,"source_sign",source_sign,"half_sign",half_sign;
print "source_parameter",par;
print "marked_root",r;
print "odd_quintic_g",g;
print "g_coefficients_0_to_5",[Coefficient(g,i) : i in [0..5]];
print "a",a;
print "b",b;
print "K2_denominator",K2den;
print "K2_numerator",K2num;
print "K1_denominator",K1den;
print "K1_numerator",K1num;
print "bidegrees_m_n",<
    <Degree(K2num,1),Degree(K2num,2)>,
    <Degree(K1num,1),Degree(K1num,2)>>;
print "cover_open_conditions s*(s^2-1)*Disc(f)*g5 != 0";
print "No saturation or elimination was run.";

quit;
