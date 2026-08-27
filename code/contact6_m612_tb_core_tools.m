//////////////////////////////////////////////////////////////////////
// Shared fixed-curve cubic-contact core tools for the T_B-halved
// contact-6 [6,12] search.  This file intentionally has no driver/quit.
//////////////////////////////////////////////////////////////////////

function M612_PrimitivePolynomial(f)
    if f eq 0 then return f; end if;
    ZZ := Integers();
    den := LCM([Denominator(c) : c in Coefficients(f)]);
    g := Parent(f)!(den*f);
    nums := [ZZ!c : c in Coefficients(g)];
    cont := GCD([Abs(n) : n in nums | n ne 0]);
    return cont gt 1 select Parent(f)!(g/cont) else g;
end function;

function M612_FixedCoreSolutions(a0,b0)
    QQ := Rationals();
    R<M,U,v> := PolynomialRing(QQ,3);
    c1 := 2*a0+6;
    c2 := a0^2+2*b0-15;
    c3 := 2*a0*b0+22;
    c4 := 2*a0+b0^2-15;
    c5 := 2*b0+6;
    B3 := c5*M+3*U;
    Delta3 := 4*c4*M+12*(U^2+v^2)-B3^2;
    F3 := M612_PrimitivePolynomial(B3*Delta3+16*v^3-8*c3*M
                                   -8*U^3-48*U*v^2);
    F2 := M612_PrimitivePolynomial(Delta3^2+64*B3*v^3-64*c2*M
                                   -192*(U^2*v^2+v^4));
    F1 := M612_PrimitivePolynomial(Delta3*v^3-4*c1*M-12*U*v^4);
    I := ideal<R | F1,F2,F3>;
    sat_ok := true;
    try
        I := Saturation(I,ideal<R | M*v*(U^2-4*v^2)>);
    catch e
        sat_ok := false;
    end try;
    dim := Dimension(I);
    if dim ne 0 then return [],dim,sat_ok,false; end if;
    try
        return Variety(I),dim,sat_ok,true;
    catch e
        return [],dim,sat_ok,false;
    end try;
end function;

function M612_VerifyCorePoint(a0,b0,M0,U0,v0)
    QQ := Rationals();
    PP<xx> := PolynomialRing(QQ);
    if M0 eq 0 or v0 eq 0 or U0^2 eq 4*v0^2 then
        return false,_,_,_,_;
    end if;
    sq,L := IsSquare(M0);
    if not sq or L eq 0 then return false,_,_,_,_; end if;
    h6 := 1+a0*xx+b0*xx^2+xx^3;
    f := h6^2-(xx-1)^6;
    if Degree(f) ne 5 or Discriminant(f) eq 0 or Evaluate(h6,1) eq 0 then
        return false,_,_,_,_;
    end if;
    q3 := xx^2+U0*xx+v0^2;
    if Discriminant(q3) eq 0 or Degree(GCD(q3,f)) gt 0 then
        return false,_,_,_,_;
    end if;
    c4 := 2*a0+b0^2-15;
    c5 := 2*b0+6;
    B3 := c5*M0+3*U0;
    Delta3 := 4*c4*M0+12*(U0^2+v0^2)-B3^2;
    h3 := (1/L)*xx^3+(B3/(2*L))*xx^2+(Delta3/(8*L))*xx+v0^3/L;
    if h3^2-f ne (1/M0)*q3^3 then return false,_,_,_,_; end if;
    J := Jacobian(HyperellipticCurve(f));
    D := J![xx-1,Evaluate(h6,1)];
    E := J![q3,h3 mod q3];
    ordD := Order(D);
    ordE := Order(E);
    independent := ordD eq 6 and ordE eq 3 and E notin {J!0,2*D,4*D};
    return independent,L,ordD,ordE,E;
end function;
