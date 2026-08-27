//////////////////////////////////////////////////////////////////////
//  Fixed-Weierstrass normalization of the exact [20]+3 square cover.
//
//  In the low-leading-coefficient chart
//
//    F_s(x)=((1-s)^2+(1-s^2)x+2s*x^2)^2-4*x^5
//           =(1-s)^4 f_t(x),   s=(t-1)/(t+1),
//
//  x=1 is a rational simple branch point and
//
//    F_s'(1)=-4(s-2)^2.
//
//  Shift z=x-1.  Since H(0)^2=q(0)^3, write
//
//    q=z^2+a*z+r^2,       H=z^3+A*z^2+B*z+r^3.
//
//  (The opposite sign of H(0) is absorbed by r -> -r.)  The z^1 equation
//  and the square condition M=L^2 are rationally parameterized by j:
//
//    L=r*j/2,
//    A=(3a-r^2*j^2)/2,
//    B=(3a*r^2-(s-2)^2*j^2)/(2r).
//
//  Thus the square condition is built in without the constant-conic
//  parameter w.  Only the z^4,z^3,z^2 equations remain in (s,r,a,j).
//  This is a birational discovery chart for the generic exact cover.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned mode then mode := "summary"; end if;
if mode notin {"summary","saturate"} then
    error "mode must be summary or saturate";
end if;
if not assigned MemGB then MemGB := 8;
elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
if not assigned PrintEquations then PrintEquations := false; end if;
if Type(PrintEquations) eq MonStgElt then
    PrintEquations := StringToLower(PrintEquations) in {"true","1","yes"};
end if;

Q := Rationals();
R<s,r,a,j> := PolynomialRing(Q,4,"grevlex");
K := FieldOfFractions(R);
KZ<z> := PolynomialRing(K);
x := z+1;

Fs := ((1-s)^2+(1-s^2)*x+2*s*x^2)^2-4*x^5;
fc := [K!Coefficient(Fs,i) : i in [0..5]];
assert fc[1] eq 0;
assert fc[2] eq -4*(s-2)^2;
assert fc[4] eq 4*(s-2)*(-s^2+2*s+5);
assert fc[5] eq 4*(s^2-5);
assert fc[6] eq -4;

M := K!(r^2*j^2/4);
A := K!((3*a-r^2*j^2)/2);
B := K!((3*a*r^2-(s-2)^2*j^2)/(2*r));
q := z^2+a*z+r^2;
H := z^3+A*z^2+B*z+r^3;
identity := H^2-q^3-M*Fs;

assert Coefficient(identity,6) eq 0;
assert Coefficient(identity,5) eq 0;
assert Coefficient(identity,1) eq 0;
assert Coefficient(identity,0) eq 0;

function PrimitiveR(f)
    if f eq 0 then return f; end if;
    den := LCM([Denominator(c) : c in Coefficients(f)]);
    ints := [Integers()!(den*c) : c in Coefficients(f)];
    cont := GCD(ints); if cont eq 0 then cont := 1; end if;
    return R!((Q!den/Q!Abs(cont))*f);
end function;

F4 := PrimitiveR(R!Numerator(Coefficient(identity,4)));
F3 := PrimitiveR(R!Numerator(Coefficient(identity,3)));
F2 := PrimitiveR(R!Numerator(Coefficient(identity,2)));

// Remove only individual powers of r introduced by denominator clearing;
// r=0 means q meets the fixed branch point and is impossible on the open
// smooth cover.
rpows := [];
for idx in [1..3] do
    ff := [F4,F3,F2][idx]; n := 0;
    while IsDivisibleBy(ff,r) do ff := ExactQuotient(ff,r); n +:= 1; end while;
    if idx eq 1 then F4 := ff; elif idx eq 2 then F3 := ff; else F2 := ff; end if;
    Append(~rpows,n);
end for;

// Direct identity is the positive control; every remaining coefficient is
// exactly one of these three equations up to the removed open unit r.
print "FIXED_WEIERSTRASS_SELF_TEST_PASS";
print "mode",mode,"equation_degrees",[TotalDegree(f) : f in [F4,F3,F2]],
      "terms",[#Terms(f) : f in [F4,F3,F2]],"removed_r_powers",rpows,
      "degrees_in_a",[Degree(f,a) : f in [F4,F3,F2]];
if PrintEquations then
    print "F4",F4;
    print "F3",F3;
    print "F2",F2;
end if;

discq := a^2-4*r^2;
resqF := PrimitiveR(R!Resultant(q,Fs));
Ds := 8*s^3-59*s^2-18*s+197;
boundary := r*j*(s-1)*(s-2)*Ds*discq*resqF;
Iraw := ideal<R | F4,F3,F2>;
print "raw_dimension",Dimension(Iraw),"basis_size",#Basis(Iraw),
      "basis_degrees",[TotalDegree(f) : f in Basis(Iraw)];

if mode eq "summary" then quit; end if;
print "SATURATING_FIXED_WEIERSTRASS";
time Isat := Saturation(Iraw,ideal<R | boundary>);
print "sat_dimension",Dimension(Isat),"basis_size",#Basis(Isat),
      "basis_degrees",[TotalDegree(f) : f in Basis(Isat)],
      "basis_terms",[#Terms(f) : f in Basis(Isat)];
print "FIXED_WEIERSTRASS_DONE";
quit;
