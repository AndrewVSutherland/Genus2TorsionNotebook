//////////////////////////////////////////////////////////////////////
// Fixed-curve exact verifier for intersecting a T_B-halved contact-6
// curve with the independent cubic-contact [3,3] core.
//
// For fixed (a,b), solve the three cubic-contact equations directly in
// (M,U,v), require M=L^2, reconstruct the order-3 Mumford class, and
// verify its order and independence on the Jacobian.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

function PrimitivePolynomial(f)
    if f eq 0 then return f; end if;
    den := LCM([Denominator(c) : c in Coefficients(f)]);
    g := Parent(f)!(den*f);
    nums := [Z!c : c in Coefficients(g)];
    cont := GCD([Abs(n) : n in nums | n ne 0]);
    return cont gt 1 select Parent(f)!(g/cont) else g;
end function;

function ContactPolynomial(a,b)
    h := 1+a*x+b*x^2+x^3;
    return h^2-(x-1)^6,h;
end function;

function FixedCoreSolutions(a0,b0)
    R<M,U,v> := PolynomialRing(Q,3);
    c1 := 2*a0+6;
    c2 := a0^2+2*b0-15;
    c3 := 2*a0*b0+22;
    c4 := 2*a0+b0^2-15;
    c5 := 2*b0+6;
    B3 := c5*M+3*U;
    Delta3 := 4*c4*M+12*(U^2+v^2)-B3^2;
    F3 := PrimitivePolynomial(B3*Delta3+16*v^3-8*c3*M
                              -8*U^3-48*U*v^2);
    F2 := PrimitivePolynomial(Delta3^2+64*B3*v^3-64*c2*M
                              -192*(U^2*v^2+v^4));
    F1 := PrimitivePolynomial(Delta3*v^3-4*c1*M-12*U*v^4);
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

function VerifyCorePoint(a0,b0,M0,U0,v0)
    if M0 eq 0 or v0 eq 0 or U0^2 eq 4*v0^2 then
        return false,_,_,_,_;
    end if;
    sq,L := IsSquare(M0);
    if not sq or L eq 0 then return false,_,_,_,_; end if;
    f,h6 := ContactPolynomial(a0,b0);
    if Degree(f) ne 5 or Discriminant(f) eq 0 or Evaluate(h6,1) eq 0 then
        return false,_,_,_,_;
    end if;
    q3 := x^2+U0*x+v0^2;
    if Discriminant(q3) eq 0 or Degree(GCD(q3,f)) gt 0 then
        return false,_,_,_,_;
    end if;
    c4 := 2*a0+b0^2-15;
    c5 := 2*b0+6;
    B3 := c5*M0+3*U0;
    Delta3 := 4*c4*M0+12*(U0^2+v0^2)-B3^2;
    h3 := (1/L)*x^3+(B3/(2*L))*x^2+(Delta3/(8*L))*x+v0^3/L;
    if h3^2-f ne (1/M0)*q3^3 then return false,_,_,_,_; end if;
    J := Jacobian(HyperellipticCurve(f));
    D := J![x-1,Evaluate(h6,1)];
    E := J![q3,h3 mod q3];
    ordD := Order(D);
    ordE := Order(E);
    independent := ordD eq 6 and ordE eq 3 and E notin {J!0,2*D,4*D};
    return independent,L,ordD,ordE,E;
end function;

tests := [
    <"positive_simple_66",Q!133/39,Q!-7/13>,
    <"tb_h10_1",Q!49/32,Q!-1>,
    <"tb_h10_2",Q!-229/96,Q!5/9>,
    <"tb_h10_3",Q!-37/8,Q!25/8>,
    <"tb_h10_4",Q!3/8,Q!-15/8>,
    <"tb_h10_5",Q!519/25,Q!-57/25>,
    <"tb_h10_6",Q!1,Q!-29/18>,
    <"tb_h10_7",Q!-11/9,Q!-5/18>
];

print "CONTACT6_M612_TB_CORE_VERIFIER";
for rec in tests do
    label,a0,b0 := Explode(rec);
    t0 := Cputime();
    sols,dim,sat_ok,variety_ok := FixedCoreSolutions(a0,b0);
    print "CURVE",label,"a",a0,"b",b0,"dim",dim,"sat",sat_ok,
          "variety",variety_ok,"rational_points",#sols,
          "seconds",Cputime(t0);
    verified := 0;
    for pt in sols do
        M0 := Q!pt[1]; U0 := Q!pt[2]; v0 := Q!pt[3];
        ok,L,ordD,ordE,E := VerifyCorePoint(a0,b0,M0,U0,v0);
        if ok then
            verified +:= 1;
            print " CORE_HIT","M",M0,"L",L,"U",U0,"v",v0,
                  "ordD",ordD,"ordE",ordE;
        end if;
    end for;
    print "CURVE_DONE",label,"verified_core_hits",verified;
end for;

quit;
