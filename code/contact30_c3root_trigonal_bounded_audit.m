//////////////////////////////////////////////////////////////////////
//  Bounded geometry audit for the degree-3 C3-root cover of the
//  simultaneous contact-(5,6) order-30 family.
//
//  This deliberately works over the original genus-zero parameter R.
//  It does not form or normalize the total-degree-32 plane projection.
//
//  Typical bounded runs:
//    magma -b mode:=summary eps:=-1 \
//      code/contact30_c3root_trigonal_bounded_audit.m
//    magma -b mode:=modgenus eps:=-1 p:=101 \
//      code/contact30_c3root_trigonal_bounded_audit.m
//    magma -b mode:=exactgenus eps:=-1 \
//      code/contact30_c3root_trigonal_bounded_audit.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);
SetMemoryLimit(4*10^9);

if not assigned mode then mode := "summary"; end if;
if not assigned eps then eps := -1;
elif Type(eps) eq MonStgElt then eps := StringToInteger(eps);
end if;
assert eps in {-1,1};
if not assigned p then p := 101;
elif Type(p) eq MonStgElt then p := StringToInteger(p);
end if;

//////////////////////////////////////////////////////////////////////
// Build the primitive cubic F(R,rho).  The return value lies in
// k[R][rho], where k is Q or F_p.  No elimination is used.
//////////////////////////////////////////////////////////////////////

function PrimitiveCubic(k, branch)
    kR<R> := PolynomialRing(k);
    K := FieldOfFractions(kR);
    KRho<rho> := PolynomialRing(K);

    t := (5*R^2-20*R+19)/(R^2-5);
    Y := -2*(5*R^2-22*R+25)/(R^2-5);
    u := t^3;
    s := t^5+t^4+(k!5/2)*t^3+(k!1/2)*t
       + branch*t*(t-k!1/2)*(t+1)*Y;
    C := (u^2+1)/(2*u);
    D := u^6+6*u^4*s-2*u^4+15*u^3*s-u*s^3+u^2;
    N := 15*u^5+90*u^4+20*u^3*s-6*u^2*s^2+231*u^3
       + 2*u^2*s-15*u*s^2+90*u^2-20*u*s+15*u-2*s;
    q := N/D;
    A := (s+q)/2;
    B := (15-s*q)/2;
    f := 2*rho^3+(A-3)*rho^2+(B+3)*rho+(C-1);

    den := &*[ Denominator(Coefficient(f,i)) : i in [0..3] ];
    cs := [ kR!(den*Coefficient(f,i)) : i in [0..3] ];
    cont := cs[1];
    for i in [2..4] do cont := GCD(cont,cs[i]); end for;
    cs := [ ExactQuotient(c,cont) : c in cs ];
    // Divide the integer/scalar content as well in characteristic zero.
    if k cmpeq Rationals() then
        zcont := GCD([ GCD([ Abs(Numerator(cc)) : cc in Coefficients(c) ])
                       : c in cs | c ne 0 ]);
        if zcont gt 1 then cs := [ c/zcont : c in cs ]; end if;
    end if;
    kRRho<z> := PolynomialRing(kR);
    F := &+[ kRRho!cs[i+1]*z^i : i in [0..3] ];
    return F,R,t,Y,u,s,D,N;
end function;

function FactorProfile(f)
    return [ <Degree(fe[1]),fe[2]> : fe in Factorization(f) ];
end function;

if mode eq "summary" then
    Q := Rationals();
    F,R,t,Y,u,s,D,N := PrimitiveCubic(Q,eps);
    QR := Parent(R);
    disc := Discriminant(F);
    lc := LeadingCoefficient(F);
    c0 := Coefficient(F,0);
    print "CONTACT30_C3ROOT_TRIGONAL_SUMMARY", "eps", eps;
    print "BIDEGREE", <Max([Degree(Coefficient(F,i)):i in [0..3]]),Degree(F)>,
          "TERMS_BY_RHO", [#Terms(Coefficient(F,i)):i in [0..3]];
    print "IRREDUCIBLE_OVER_QR", IsIrreducible(F);
    print "LEADING_COEFFICIENT_DEGREE",Degree(lc),
          "PROFILE",FactorProfile(lc),"FACTORIZATION",Factorization(lc);
    print "CONSTANT_COEFFICIENT_DEGREE",Degree(c0),
          "PROFILE",FactorProfile(c0);
    print "DISCRIMINANT_DEGREE",Degree(disc),
          "TERMS",#Terms(disc),
          "PROFILE",FactorProfile(disc);
    print "DISCRIMINANT_FACTORIZATION",Factorization(disc);
    print "DISCRIMINANT_SQUAREFREE_DEGREE",Degree(SquarefreePart(disc));
    print "GCD_DISC_DERIVATIVE_DEGREE",Degree(GCD(disc,Derivative(disc)));

    // The two signs are the same cover: the deck involution of the
    // degree-two conic parametrization is R |-> (11R-25)/(5R-11).
    K := FieldOfFractions(QR);
    Ri := (11*R-25)/(5*R-11);
    print "CONIC_INVOLUTION_SQUARE", Evaluate(
        Numerator(Evaluate(Ri,Ri)-R),R) eq 0;
    print "CONIC_INVOLUTION_T_FIXED", t eq Evaluate(t,Ri);
    print "CONIC_INVOLUTION_Y_NEGATED", Y eq -Evaluate(Y,Ri);

    // Base-locus profiles used to distinguish genuine branch points from
    // poles or degenerate fibers in this affine cubic model.
    bden := R^2-5;
    tnum := 5*R^2-20*R+19;
    tminus1num := Numerator(t-1);
    tplus1num := Numerator(t+1);
    Dnum := Numerator(D);
    print "BOUNDARY_R_DEN",Factorization(bden);
    print "BOUNDARY_U_ZERO",Factorization(tnum);
    print "BOUNDARY_U2_MINUS1",Factorization(tminus1num*tplus1num);
    print "BOUNDARY_RECOVERY_D_DEGREE",Degree(Dnum),
          "PROFILE",FactorProfile(Dnum),"FACTORIZATION",Factorization(Dnum);
    print "DISC_GCD_R_DEN",FactorProfile(GCD(disc,bden));
    print "DISC_GCD_U_ZERO",FactorProfile(GCD(disc,tnum));
    print "DISC_GCD_U2_MINUS1",FactorProfile(GCD(disc,tminus1num*tplus1num));
    print "DISC_GCD_RECOVERY_D",FactorProfile(GCD(disc,Dnum));
    print "CONTACT30_C3ROOT_TRIGONAL_SUMMARY_DONE";
    quit;
end if;

if mode eq "modgenus" then
    assert IsPrime(p) and p gt 5;
    k := GF(p);
    F,R,t,Y,u,s,D,N := PrimitiveCubic(k,eps);
    kR := Parent(R);
    K<T> := FunctionField(k);
    KT<z> := PolynomialRing(K);
    FK := KT![ K!Coefficient(F,i) : i in [0..3] ];
    print "CONTACT30_C3ROOT_TRIGONAL_MODGENUS", "eps",eps,"p",p;
    print "BIDEGREE",<Max([Degree(Coefficient(F,i)):i in [0..3]]),Degree(F)>;
    print "FACTOR_DEGREES_OVER_FP_R",FactorProfile(F);
    disc := Discriminant(F);
    print "DISCRIMINANT_DEGREE",Degree(disc),
          "PROFILE",FactorProfile(disc),
          "SQUAREFREE_DEGREE",Degree(SquarefreePart(disc));
    assert IsIrreducible(FK);
    time L<w> := FunctionField(FK : Check := true);
    time g := Genus(L);
    time Diff := DifferentDivisor(L);
    print "GENUS",g,"DIFFERENT_DEGREE",Degree(Diff);
    print "RH_GENUS",1-3+Degree(Diff) div 2;
    print "CONTACT30_C3ROOT_TRIGONAL_MODGENUS_DONE";
    quit;
end if;

assert mode eq "exactgenus";
Q := Rationals();
F,R,t,Y,u,s,D,N := PrimitiveCubic(Q,eps);
K<T> := FunctionField(Q);
KT<z> := PolynomialRing(K);
FK := KT![ K!Coefficient(F,i) : i in [0..3] ];
print "CONTACT30_C3ROOT_TRIGONAL_EXACTGENUS", "eps",eps;
assert IsIrreducible(FK);
time L<w> := FunctionField(FK : Check := true);
time g := Genus(L);
time Diff := DifferentDivisor(L);
print "GENUS",g,"DIFFERENT_DEGREE",Degree(Diff);
print "RH_GENUS",1-3+Degree(Diff) div 2;
Ps,ns := Support(Diff);
print "DIFFERENT_SUPPORT_PROFILE",
      [<Degree(Ps[i]),RamificationIndex(Ps[i]),ns[i]>:i in [1..#Ps]];
print "CONTACT30_C3ROOT_TRIGONAL_EXACTGENUS_DONE";
quit;
