//////////////////////////////////////////////////////////////////////
// Diagnostic audit for contact6_m612_dual_split_core_search.m.
//
// This repeats only the rational slice reconstruction and classifies every
// square-M point before the dual discriminant-cover filter.  In particular,
// it distinguishes a genuine dependent order-3 direction from a spurious
// solution introduced while clearing the recovery denominator.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned height then height:=6;
elif Type(height) eq MonStgElt then height:=StringToInteger(height); end if;
if not assigned example_limit then example_limit:=20;
elif Type(example_limit) eq MonStgElt then
    example_limit:=StringToInteger(example_limit);
end if;

Q:=Rationals();
P<x>:=PolynomialRing(Q);
load "code/contact6_m612_tb_core_tools.m";

function RationalParametersOfHeight(B)
    vals:=[]; seen:={};
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num,den) ne 1 then continue; end if;
            q:=Q!num/den;
            if Sprint(q) notin seen then
                Include(~seen,Sprint(q)); Append(~vals,q);
            end if;
        end for;
    end for;
    return vals;
end function;

function GenericSliceSolutions(b0,v0)
    R<M,U>:=PolynomialRing(Q,2,"grevlex");
    K:=FieldOfFractions(R); PA<a>:=PolynomialRing(K);
    b:=K!b0; v:=K!v0;
    c1:=2*a+6;
    c2:=a^2+2*b-15;
    c3:=2*a*b+22;
    c4:=2*a+b^2-15;
    c5:=2*b+6;
    B3:=c5*M+3*U;
    Delta3:=4*c4*M+12*(U^2+v^2)-B3^2;
    F3:=B3*Delta3+16*v^3-8*c3*M-8*U^3-48*U*v^2;
    F2:=Delta3^2+64*B3*v^3-64*c2*M-192*(U^2*v^2+v^4);
    F1:=Delta3*v^3-4*c1*M-12*U*v^4;
    D:=Coefficient(F1,1); N:=-Coefficient(F1,0);
    if D eq 0 then return [],false,-2; end if;

    function SubNum(poly)
        d:=Degree(poly); value:=K!0;
        for i in [0..d] do
            value+:=Coefficient(poly,i)*N^i*D^(d-i);
        end for;
        return M612_PrimitivePolynomial(R!Numerator(value));
    end function;

    G2:=SubNum(F2); G3:=SubNum(F3);
    boundary:=M*(U^2-4*v0^2)*(R!(N+(b0+2)*D));
    I:=ideal<R|G2,G3>; sat_ok:=true;
    try I:=Saturation(I,ideal<R|boundary>);
    catch e sat_ok:=false; end try;
    dim:=Dimension(I);
    if not sat_ok or dim ne 0 then return [],sat_ok,dim; end if;
    pts:=[];
    try pts:=Variety(I);
    catch e return [],sat_ok,-3; end try;
    out:=[];
    for pt in pts do
        M0:=Q!pt[1]; U0:=Q!pt[2];
        if M0 eq 0 then continue; end if;
        sq,L0:=IsSquare(M0);
        if not sq or L0 eq 0 then continue; end if;
        Nnum:=R!Numerator(N); Nden:=R!Denominator(N);
        Dnum:=R!Numerator(D); Dden:=R!Denominator(D);
        denN:=Evaluate(Nden,<M0,U0>); denD:=Evaluate(Dden,<M0,U0>);
        if denN eq 0 or denD eq 0 or Evaluate(Dnum,<M0,U0>) eq 0 then
            continue;
        end if;
        a0:=(Evaluate(Nnum,<M0,U0>)/denN) /
            (Evaluate(Dnum,<M0,U0>)/denD);
        Append(~out,<a0,b0,M0,U0,v0>);
    end for;
    return out,sat_ok,dim;
end function;

function SpecialSliceSolutions(b0)
    R<L,U>:=PolynomialRing(Q,2,"grevlex");
    K:=FieldOfFractions(R); PA<a>:=PolynomialRing(K);
    b:=K!b0; v:=K!1; M:=L^2;
    c1:=2*a+6;
    c2:=a^2+2*b-15;
    c3:=2*a*b+22;
    c4:=2*a+b^2-15;
    c5:=2*b+6;
    B3:=c5*M+3*U;
    Delta3:=4*c4*M+12*(U^2+v^2)-B3^2;
    F3:=B3*Delta3+16*v^3-8*c3*M-8*U^3-48*U*v^2;
    F2:=Delta3^2+64*B3*v^3-64*c2*M-192*(U^2*v^2+v^4);
    F1:=Delta3*v^3-4*c1*M-12*U*v^4;
    if Degree(F1) ne 0 then return [],false,-2; end if;
    H:=M612_PrimitivePolynomial(R!Coefficient(F1,0));
    D:=Coefficient(F3,1); N:=-Coefficient(F3,0);
    if D eq 0 then return [],false,-2; end if;
    d:=Degree(F2); value:=K!0;
    for i in [0..d] do value+:=Coefficient(F2,i)*N^i*D^(d-i); end for;
    J:=M612_PrimitivePolynomial(R!Numerator(value));
    DR:=R!Numerator(D); NR:=R!Numerator(N);
    I:=ideal<R|H,J>; sat_ok:=true;
    boundary:=L*(U^2-4)*(b0+3)*DR*(NR+(b0+2)*DR);
    try I:=Saturation(I,ideal<R|boundary>);
    catch e sat_ok:=false; end try;
    dim:=Dimension(I);
    if not sat_ok or dim ne 0 then return [],sat_ok,dim; end if;
    pts:=[];
    try pts:=Variety(I);
    catch e return [],sat_ok,-3; end try;
    out:=[];
    for pt in pts do
        L0:=Q!pt[1]; U0:=Q!pt[2];
        if L0 eq 0 then continue; end if;
        d0:=Evaluate(DR,<L0,U0>);
        if d0 eq 0 then continue; end if;
        a0:=Evaluate(NR,<L0,U0>)/d0;
        Append(~out,<a0,b0,L0^2,U0,Q!1>);
    end for;
    return out,sat_ok,dim;
end function;

function Diagnosis(a,b,M,U,v)
    if M eq 0 or v eq 0 then return "parameter_boundary"; end if;
    sq,L:=IsSquare(M);
    if not sq or L eq 0 then return "nonsquare_M"; end if;
    h:=1+a*x+b*x^2+x^3;
    f:=h^2-(x-1)^6;
    if Degree(f) ne 5 then return "degree_drop"; end if;
    if Discriminant(f) eq 0 then return "singular_curve"; end if;
    if Evaluate(h,1) eq 0 then return "marked_boundary"; end if;
    q:=x^2+U*x+v^2;
    if Discriminant(q) eq 0 then return "double_q"; end if;
    if Degree(GCD(q,f)) gt 0 then return "q_meets_branch"; end if;
    c4:=2*a+b^2-15;
    c5:=2*b+6;
    B3:=c5*M+3*U;
    Delta3:=4*c4*M+12*(U^2+v^2)-B3^2;
    h3:=(1/L)*x^3+(B3/(2*L))*x^2+(Delta3/(8*L))*x+v^3/L;
    if h3^2-f ne (1/M)*q^3 then return "contact_identity_failure"; end if;
    J:=Jacobian(HyperellipticCurve(f));
    D:=J![x-1,Evaluate(h,1)];
    E:=J![q,h3 mod q];
    ordD:=Order(D); ordE:=Order(E);
    if ordD ne 6 then return "marked_order_not_6"; end if;
    if ordE ne 3 then return "contact_order_not_3"; end if;
    if E eq J!0 then return "zero_class"; end if;
    if E eq 2*D then return "dependent_2D"; end if;
    if E eq 4*D then return "dependent_4D"; end if;
    return "independent_core";
end function;

// The shared verifier must recognize the three presentations on the known
// simple [6,6] source before the slice audit is meaningful.
controls:=[
    <Q!133/39,Q!-7/13,Q!9/256,Q!5/12,Q!5/6>,
    <Q!133/39,Q!-7/13,Q!841/256,Q!-9/4,Q!5/2>,
    <Q!133/39,Q!-7/13,Q!169/16,Q!-17/2,Q!-5/4>
];
for rec in controls do
    a,b,M,U,v:=Explode(rec);
    ok,L,ordD,ordE,E:=M612_VerifyCorePoint(a,b,M,U,v);
    print "POSITIVE_CONTROL","a",a,"b",b,"M",M,"U",U,"v",v,
          "ok",ok,"diagnosis",Diagnosis(a,b,M,U,v),
          "ordD",ordD,"ordE",ordE;
end for;

params:=RationalParametersOfHeight(height);
checked:=0; raw:=0; cover:=0; printed:=0;
diagnoses:=AssociativeArray();
for b0 in params do
    for v0 in params do
        if v0 eq 0 then continue; end if;
        checked+:=1;
        if v0 eq 1 then
            recs,sat_ok,dim:=SpecialSliceSolutions(b0);
        else
            recs,sat_ok,dim:=GenericSliceSolutions(b0,v0);
        end if;
        for rec in recs do
            raw+:=1;
            a,b,M,U,v:=Explode(rec);
            DB:=(a-3)^2-8*(b+3);
            DC:=(b-3)^2-8*(a+3);
            oncover:=(DB ne 0 and IsSquare(DB)) or
                     (DC ne 0 and IsSquare(DC)) or
                     (DB*DC ne 0 and IsSquare(DB*DC));
            if oncover then cover+:=1; end if;
            why:=Diagnosis(a,b,M,U,v);
            if IsDefined(diagnoses,why) then
                diagnoses[why]+:=1;
            else
                diagnoses[why]:=1;
            end if;
            if printed lt example_limit then
                print "AUDIT_POINT","a",a,"b",b,"M",M,"U",U,"v",v,
                      "cover",oncover,"diagnosis",why;
                printed+:=1;
            end if;
            if why eq "independent_core" then
                print "INDEPENDENT_CORE","a",a,"b",b,"M",M,"U",U,
                      "v",v,"cover",oncover,"DB",DB,"DC",DC;
            end if;
            if v eq 1 then
                print "SPECIAL_V1_POINT","a",a,"b",b,"M",M,"U",U,
                      "cover",oncover,"diagnosis",why;
            end if;
        end for;
    end for;
end for;

print "AUDIT_DONE","height",height,"checked",checked,"raw",raw,
      "cover",cover,"diagnoses",
      Sort([<key,diagnoses[key]>:key in Keys(diagnoses)]);
quit;
