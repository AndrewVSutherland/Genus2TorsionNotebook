//////////////////////////////////////////////////////////////////////
// Optimized exact [6,12] search on the T0-halved contact-6 [3,3] core.
//
// Fix (omega,nu).  The cubic-contact equations form a zero-dimensional
// system in (s,M,U), where
//
//   b=s-3, a=s*omega^2-3, M=L^2.
//
// For each rational core point, solve the independent quartic
// Fhalf(s,omega,W)=0 and require W=m^2.  This block decomposition is
// equivalent to solving all four equations in (s,W,M,U), but is much
// faster and exposes the two square lifts separately.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned height then height := 2;
elif Type(height) eq MonStgElt then height := StringToInteger(height); end if;
if not assigned max_slices then max_slices := 0;
elif Type(max_slices) eq MonStgElt then max_slices := StringToInteger(max_slices); end if;
if not assigned max_hits then max_hits := 10;
elif Type(max_hits) eq MonStgElt then max_hits := StringToInteger(max_hits); end if;
if not assigned progress_interval then progress_interval := 1;
elif Type(progress_interval) eq MonStgElt then progress_interval := StringToInteger(progress_interval); end if;

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

function RationalParametersOfHeight(B)
    vals := [];
    seen := {};
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num,den) ne 1 then continue; end if;
            z := Q!num/den;
            if Sprint(z) notin seen then
                Include(~seen,Sprint(z));
                Append(~vals,z);
            end if;
        end for;
    end for;
    return vals;
end function;

function ContactPolynomial(a,b)
    h := 1+a*x+b*x^2+x^3;
    return h^2-(x-1)^6,h;
end function;

function GoodPolynomial(f)
    return Degree(f) eq 5 and Discriminant(f) ne 0;
end function;

function FactorDegrees(f)
    return Sort([Degree(ff[1]) : ff in Factorization(f)]);
end function;

function IntegralModel(f)
    d := LCM([Denominator(Coefficient(f,i)) : i in [0..Degree(f)]]);
    return P!(d^2*f),d;
end function;

function OnRequiredBoundaryAtPrime(a,b,p)
    if Numerator(Denominator(a)) mod p eq 0 or
       Numerator(Denominator(b)) mod p eq 0 then
        return true;
    end if;
    K := GF(p);
    aa := K!Numerator(a)/K!Denominator(a);
    bb := K!Numerator(b)/K!Denominator(b);
    DB := (aa-3)^2-8*(bb+3);
    DC := (bb-3)^2-8*(aa+3);
    RR := ((bb+3)*(aa+3)-4)^2
          -(bb^2-2*aa-3)*(aa^2-2*bb-3);
    return (bb+3)*(aa+3)*DB*DC*RR eq 0;
end function;

function PassBoundary57(a,b)
    return OnRequiredBoundaryAtPrime(a,b,5) and
           OnRequiredBoundaryAtPrime(a,b,7);
end function;

function SimpleCertificate(f)
    C := HyperellipticCurve(f);
    for p in [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97] do
        try
            fp := ChangeRing(f,GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then continue; end if;
            Lp := LPolynomial(ChangeRing(C,GF(p)));
            fac := Factorization(Lp);
            if #fac eq 1 and fac[1][2] eq 1 and Degree(fac[1][1]) eq 4 then
                return true,p,Lp;
            end if;
        catch e
            continue;
        end try;
    end for;
    return false,0,P!0;
end function;

function PrimitivePolynomial(f)
    if f eq 0 then return f; end if;
    den := LCM([Denominator(c) : c in Coefficients(f)]);
    g := Parent(f)!(den*f);
    nums := [Z!c : c in Coefficients(g)];
    cont := GCD([Abs(n) : n in nums | n ne 0]);
    return cont gt 1 select Parent(f)!(g/cont) else g;
end function;

function CoreSolutions(omega0,nu0)
    R<s,M,U> := PolynomialRing(Q,3);
    a := s*omega0^2-3;
    b := s-3;
    c5 := 2*s;
    c4 := b^2+2*a-15;
    c3 := 2*a*b+22;
    c2 := a^2+2*b-15;
    c1 := 2*(a+3);
    B3 := c5*M+3*U;
    Delta3 := 4*c4*M+12*(U^2+nu0^2)-B3^2;
    F3 := PrimitivePolynomial(B3*Delta3+16*nu0^3-8*c3*M
                              -8*U^3-48*U*nu0^2);
    F2 := PrimitivePolynomial(Delta3^2+64*B3*nu0^3-64*c2*M
                              -192*(U^2*nu0^2+nu0^4));
    F1 := PrimitivePolynomial(Delta3*nu0^3-4*c1*M-12*U*nu0^4);
    I := ideal<R | F1,F2,F3>;
    boundary := s*M*(U^2-4*nu0^2)*(s*(omega0^2+1)-4);
    sat_ok := true;
    try
        I := Saturation(I,ideal<R|boundary>);
    catch e
        sat_ok := false;
    end try;
    dim := Dimension(I);
    sols := [];
    if dim eq 0 then
        try
            sols := Variety(I);
        catch e
            return [],dim,sat_ok,false;
        end try;
    end if;
    return sols,dim,sat_ok,true;
end function;

function HalfQuarticRoots(s0,omega0)
    PW<W> := PolynomialRing(Q);
    a := s0*omega0^2-3;
    b := s0-3;
    c4 := b^2+2*a-15;
    c3 := 2*a*b+22;
    c2 := a^2+2*b-15;
    N := c2-omega0*c4+omega0*W;
    R8 := 8*s0*(c3-4*s0*omega0)-(c4-W)^2;
    F := R8^2-256*s0^2*W*N;
    return Roots(F);
end function;

function VerifyCandidate(s0,omega0,W0,M0,U0,nu0)
    if s0 eq 0 or W0 eq 0 or M0 eq 0 then return false,_,_,_,_,_,_,_; end if;
    sqW,m := IsSquare(W0);
    sqM,L := IsSquare(M0);
    if not sqW or not sqM then return false,_,_,_,_,_,_,_; end if;
    a := s0*omega0^2-3;
    b := s0-3;
    f,h := ContactPolynomial(a,b);
    if not GoodPolynomial(f) or FactorDegrees(f) ne [1,2,2] then
        return false,_,_,_,_,_,_,_;
    end if;
    if not PassBoundary57(a,b) then return false,_,_,_,_,_,_,_; end if;

    c4 := b^2+2*a-15;
    c3 := 2*a*b+22;
    c2 := a^2+2*b-15;
    u0 := (c4-W0)/(4*s0);
    R8 := 8*s0*(c3-4*s0*omega0)-(c4-W0)^2;
    n0 := R8/(16*s0*m);
    qh := x^2+u0*x+omega0;
    vh := (m*u0-n0)*x+m*omega0;

    // Verify the tangent identity before entering the Jacobian.
    g := ExactQuotient(f,x);
    if g-x*(m*x+n0)^2 ne 2*s0*qh^2 then
        return false,_,_,_,_,_,_,_;
    end if;

    // Verify the cubic contact and the marked classes on the rational model.
    B3 := 2*s0*M0+3*U0;
    Delta3 := 4*c4*M0+12*(U0^2+nu0^2)-B3^2;
    q3 := x^2+U0*x+nu0^2;
    h3 := (1/L)*x^3+(B3/(2*L))*x^2+(Delta3/(8*L))*x+nu0^3/L;
    J := Jacobian(HyperellipticCurve(f));
    D := J![x-1,Evaluate(h,1)];
    E := J![q3,h3 mod q3];
    ordD := Order(D);
    ordE := Order(E);
    if ordD ne 6 or ordE ne 3 then return false,_,_,_,_,_,_,_; end if;

    fI,d := IntegralModel(f);
    JI := Jacobian(HyperellipticCurve(fI));
    T0 := JI![x,0];
    H4 := JI![qh,d*vh];
    half_ok := 2*H4 eq T0 and Order(H4) eq 4;
    if not half_ok then return false,_,_,_,_,_,_,_; end if;

    simple,pcert,Lp := SimpleCertificate(fI);
    if not simple then
        return true,[],false,0,ordD,ordE,fI,H4;
    end if;
    G,mp := TorsionSubgroup(JI);
    inv := Invariants(G);
    return true,inv,true,pcert,ordD,ordE,fI,H4;
end function;

params := RationalParametersOfHeight(height);
checked := 0;
dim0 := 0;
positive_dim := 0;
core_points := 0;
M_squares := 0;
half_roots := 0;
W_squares := 0;
boundary57 := 0;
smooth_boundary := 0;
verified := 0;
simple_count := 0;
exact_tests := 0;
hits := [];
total_start := Cputime();

print "CONTACT6_M612_T0_CORE_SLICE_SEARCH";
print "height",height,"parameter_count",#params,"max_slices",max_slices;

for omega0 in params do
    if omega0 in {-1,0,1} then continue; end if;
    for nu0 in params do
        if nu0 eq 0 or nu0^3 eq 1 then continue; end if;
        if max_slices gt 0 and checked ge max_slices then break omega0; end if;
        if #hits ge max_hits then break omega0; end if;
        checked +:= 1;
        t0 := Cputime();
        sols,dim,sat_ok,variety_ok := CoreSolutions(omega0,nu0);
        elapsed := Cputime(t0);
        if dim eq 0 then dim0 +:= 1; else positive_dim +:= 1; end if;
        if progress_interval gt 0 and checked mod progress_interval eq 0 then
            print "SLICE",checked,"omega",omega0,"nu",nu0,"dim",dim,
                  "sat",sat_ok,"variety",variety_ok,"solutions",#sols,
                  "seconds",elapsed;
        end if;
        for pt in sols do
            core_points +:= 1;
            s0 := Q!pt[1]; M0 := Q!pt[2]; U0 := Q!pt[3];
            okM,L0 := IsSquare(M0);
            if not okM or M0 eq 0 then continue; end if;
            M_squares +:= 1;
            a0 := s0*omega0^2-3;
            b0 := s0-3;
            if not PassBoundary57(a0,b0) then continue; end if;
            boundary57 +:= 1;
            f0,h0 := ContactPolynomial(a0,b0);
            if not GoodPolynomial(f0) or FactorDegrees(f0) ne [1,2,2] then
                print "BOUNDARY_SINGULAR_OR_WRONG_FACTOR","omega",omega0,
                      "nu",nu0,"s",s0,"M",M0,"U",U0,
                      "factor_degrees",FactorDegrees(f0);
                continue;
            end if;
            smooth_boundary +:= 1;
            rootsW := HalfQuarticRoots(s0,omega0);
            half_roots +:= #rootsW;
            print "BOUNDARY_CORE","omega",omega0,"nu",nu0,"s",s0,
                  "M",M0,"L",L0,"U",U0,"half_roots",rootsW;
            for rt in rootsW do
                W0 := Q!rt[1];
                okW,m0 := IsSquare(W0);
                if okW then
                    print "HALF_ROOT","W",W0,"multiplicity",rt[2],
                          "is_square",true,"square_root",m0;
                else
                    print "HALF_ROOT","W",W0,"multiplicity",rt[2],
                          "is_square",false;
                end if;
                if not okW or W0 eq 0 then continue; end if;
                W_squares +:= 1;
                ok,inv,simple,pcert,ordD,ordE,fI,H4 :=
                    VerifyCandidate(s0,omega0,W0,M0,U0,nu0);
                if not ok then continue; end if;
                verified +:= 1;
                if simple then
                    simple_count +:= 1;
                    exact_tests +:= 1;
                end if;
                print "CANDIDATE","omega",omega0,"nu",nu0,"s",s0,
                      "W",W0,"M",M0,"U",U0,"simple",simple,
                      "pcert",pcert,"torsion",inv,"ordD",ordD,"ordE",ordE;
                print " curve",fI;
                if simple and inv eq [6,12] then
                    Append(~hits,<omega0,nu0,s0,W0,M0,U0,inv,pcert,fI>);
                    print "HIT_6_12",hits[#hits];
                end if;
            end for;
        end for;
    end for;
end for;

print "DONE";
print "checked_slices",checked,"dim0",dim0,"positive_dim",positive_dim,
      "core_points",core_points,"M_squares",M_squares,
      "boundary57",boundary57,"smooth_boundary",smooth_boundary,
      "half_roots",half_roots,
      "W_squares",W_squares,"verified",verified,
      "simple",simple_count,"exact_tests",exact_tests,"hits",#hits,
      "seconds",Cputime(total_start);
for H in hits do print "H",H; end for;

quit;
