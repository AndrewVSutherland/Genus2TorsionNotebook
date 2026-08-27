//////////////////////////////////////////////////////////////////////
// The rational two-Weierstrass-root surface in the contact-7 family.
//
// Put
//     h = 1 - (7/2)x + a*x^2 + b*x^3,
//     f = (h^2 + (x-1)^7)/x^2.
// For r=1-s^2 and w=1-t^2, impose
//     h(r)=eps*s^7,  h(w)=delta*t^7.
// The signs are redundant: replace s by eps*s and t by delta*t.
// If
//     A(u)=(2u^5+4u^4+6u^3+8u^2+10u+5)/(2(u+1)^2),
// then the cancelled parametrization is
//     b=(A(t)-A(s))/(s^2-t^2),  a=A(s)-b(1-s^2).
// Its only parameter poles are s=-1, t=-1, and s=-t.  The
// resulting genus-2 curve has the contact-7 class and two independent
// rational Weierstrass classes, hence visible Z/2 x Z/14.
//
// Modes:
//   magma -b mode:=symbolic code/contact7_two_root_surface.m
//   magma -b mode:=sample code/contact7_two_root_surface.m
// Optional sample parameters: s, t, eps, delta.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned mode then mode := "sample"; end if;
if not assigned s then s := 2; elif Type(s) eq MonStgElt then s := StringToInteger(s); end if;
if not assigned t then t := 3; elif Type(t) eq MonStgElt then t := StringToInteger(t); end if;
if not assigned eps then eps := 1; elif Type(eps) eq MonStgElt then eps := StringToInteger(eps); end if;
if not assigned delta then delta := 1; elif Type(delta) eq MonStgElt then delta := StringToInteger(delta); end if;

Q := Rationals();
P<x> := PolynomialRing(Q);

function IntegralModel(f)
    den := 1;
    for i in [0..Degree(f)] do
        den := LCM(den, Denominator(Coefficient(f,i)));
    end for;
    // Smallest positive d such that d^2*f is integral.
    d := 1;
    for q in Factorization(den) do
        d *:= q[1]^Ceiling(Q!q[2]/2);
    end for;
    return P!(d^2*f), d;
end function;

function FrobeniusPolynomial(C,p)
    ef := EulerFactor(C,p);
    d := Degree(ef);
    return &+[ Q!Coefficient(ef,i)*x^(d-i) : i in [0..d] ];
end function;

function AbsoluteSimplicityCertificate(C)
    plist := [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97];
    fC,hC := HyperellipticPolynomials(C);
    // A D4 Frobenius polynomial is the standard genus-2 geometric-
    // simplicity certificate used elsewhere in this project.
    for p in plist do
        try
            fp := ChangeRing(fC,GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then continue; end if;
            Phi := FrobeniusPolynomial(C,p);
            fac := Factorization(Phi);
            if #fac ne 1 or fac[1][2] ne 1 or Degree(fac[1][1]) ne 4 then
                continue;
            end if;
            G := GaloisGroup(Phi);
            desc := "";
            try desc := TransitiveGroupDescription(G); catch e2 desc := "unknown"; end try;
            if Order(G) eq 8 and desc eq "D(4)" then
                K<pi> := NumberField(Phi);
                power_degrees := [ Degree(MinimalPolynomial(pi^n)) : n in [2..12] ];
                if (Integers()!Coefficient(Phi,2)) mod p ne 0 and
                   &and[ d eq 4 : d in power_degrees ] then
                    return true,"ordinary_D4+root_power_2_12",p,Phi;
                end if;
            end if;
        catch e
            continue;
        end try;
    end for;
    // Fallback: no Frobenius-root power through the abelian-surface
    // bound 12 drops degree.
    for p in plist do
        try
            fp := ChangeRing(fC,GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then continue; end if;
            Phi := FrobeniusPolynomial(C,p);
            fac := Factorization(Phi);
            if #fac ne 1 or fac[1][2] ne 1 or Degree(fac[1][1]) ne 4 then
                continue;
            end if;
            K<pi> := NumberField(Phi);
            if &and[ Degree(MinimalPolynomial(pi^n)) eq 4 : n in [2..12] ] then
                return true,"root_power",p,Phi;
            end if;
        catch e
            continue;
        end try;
    end for;
    return false,"none",0,P!0;
end function;

function TwoRootModel(ss,tt,ee,dd)
    ss := Q!ss; tt := Q!tt; ee := Q!ee; dd := Q!dd;
    if ee notin {Q!-1,Q!1} or dd notin {Q!-1,Q!1} then
        return false,"sign",P!0,P!0,Q!0,Q!0,Q!0,Q!0;
    end if;
    // Absorb the signs into signed square roots.
    S := ee*ss;
    T := dd*tt;
    r := 1-S^2;
    w := 1-T^2;
    if S eq -1 or T eq -1 or S eq -T then
        return false,"parameter pole",P!0,P!0,Q!0,Q!0,r,w;
    end if;
    if S eq T then
        return false,"coincident roots",P!0,P!0,Q!0,Q!0,r,w;
    end if;
    AS := (2*S^5+4*S^4+6*S^3+8*S^2+10*S+5)/(2*(S+1)^2);
    AT := (2*T^5+4*T^4+6*T^3+8*T^2+10*T+5)/(2*(T+1)^2);
    b := (AT-AS)/(S^2-T^2);
    a := AS-b*r;
    h := 1-(Q!7/2)*x+a*x^2+b*x^3;
    num := h^2+(x-1)^7;
    assert Coefficient(num,0) eq 0 and Coefficient(num,1) eq 0;
    f := ExactQuotient(num,x^2);
    assert h^2-x^2*f eq -(x-1)^7;
    assert Evaluate(h,r) eq S^7 and Evaluate(h,w) eq T^7;
    assert Evaluate(f,r) eq 0 and Evaluate(f,w) eq 0;
    return true,"",f,h,a,b,r,w;
end function;

procedure PrintModel(ss,tt,ee,dd,DoExact)
    ok,msg,f,h,a,b,r,w := TwoRootModel(ss,tt,ee,dd);
    print "parameters",ss,tt,ee,dd;
    print "model_ok",ok,"boundary",msg;
    if not ok then return; end if;
    print "r",r,"w",w;
    print "a",a;
    print "b",b;
    print "h",h;
    print "f",f;
    print "factorization",Factorization(f);
    print "contact_identity",h^2-x^2*f eq -(x-1)^7;
    print "root_identities",Evaluate(f,r) eq 0,Evaluate(f,w) eq 0;
    print "marked_y",Evaluate(h,1);
    print "degree",Degree(f),"discriminant",Discriminant(f);
    if Degree(f) ne 5 or Discriminant(f) eq 0 or Evaluate(h,1) eq 0 then
        print "smooth_contact_fiber",false;
        return;
    end if;
    print "smooth_contact_fiber",true;
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    D7 := J![x-1,Evaluate(h,1)];
    Dr := J![x-r,0];
    Dw := J![x-w,0];
    print "visible_orders",Order(D7),Order(Dr),Order(Dw),Order(Dr+Dw);
    assert Order(D7) eq 7;
    assert Order(Dr) eq 2 and Order(Dw) eq 2 and Order(Dr+Dw) eq 2;
    assert Dr ne Dw;
    if DoExact then
        fI,L := IntegralModel(f);
        CI := HyperellipticCurve(fI);
        JI := Jacobian(CI);
        G,phi := TorsionSubgroup(JI);
        print "integral_scale",L;
        print "integral_model",fI;
        print "torsion_invariants",Invariants(G),"torsion_order",#G;
        try
            Cmin,tomin := MinimalWeierstrassModel(CI);
            Cred,tored := ReducedModel(Cmin : Smallest:=true, Height:=true);
            fred,hred := HyperellipticPolynomials(Cred);
            print "reduced_minimal_model_f",fred;
            print "reduced_minimal_model_h",hred;
        catch e
            print "reduced_minimal_model_unavailable",e`Object;
        end try;
        simple,kind,p,Phi := AbsoluteSimplicityCertificate(CI);
        print "geometrically_simple",simple,"certificate",kind,"prime",p;
        if simple then
            print "Frobenius",Phi;
            Kpi<pi> := NumberField(Phi);
            print "ordinary",(Integers()!Coefficient(Phi,2)) mod p ne 0;
            print "Frobenius_Galois_group",TransitiveGroupDescription(GaloisGroup(Phi));
            print "root_power_degrees_2_to_12",
                  [ Degree(MinimalPolynomial(pi^n)) : n in [2..12] ];
        end if;
    end if;
end procedure;

if mode eq "symbolic" then
    F<S,T> := RationalFunctionField(Q,2);
    rr := 1-S^2; ww := 1-T^2;
    AS := (2*S^5+4*S^4+6*S^3+8*S^2+10*S+5)/(2*(S+1)^2);
    AT := (2*T^5+4*T^4+6*T^3+8*T^2+10*T+5)/(2*(T+1)^2);
    bb := (AT-AS)/(S^2-T^2);
    aa := AS-bb*rr;
    print "r =",rr;
    print "w =",ww;
    print "A(S) =",AS;
    print "A(T) =",AT;
    print "a =",aa;
    print "b =",bb;
    assert aa*rr^2+bb*rr^3 eq S^7-1+(F!7/2)*rr;
    assert aa*ww^2+bb*ww^3 eq T^7-1+(F!7/2)*ww;
    print "sign_absorption: (s,eps) is replaced by eps*s";
    print "parameter_poles: S=-1, T=-1, S=-T";
    print "singular_boundary includes S=0, T=0, S=T and residual discriminant factors";
else
    PrintModel(Q!s,Q!t,Q!eps,Q!delta,true);
end if;

quit;
