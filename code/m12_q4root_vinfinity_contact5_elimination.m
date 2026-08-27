//////////////////////////////////////////////////////////////////////
// Exact contact-5 elimination on the v=infinity boundary of the
// rational-Q4-root cover of M(12).
//
// The Q4 leading coefficient is a*(4*a+1).  The a=0 component is
// degenerate, while a=-1/4 gives the valid one-parameter odd model
//
//   f = 16*((x-r)^2*(T+1)^2 - x^2*T*(T+1)),
//   T = -x^2/4-x+r.
//
// Its marked point (0,-4*r*(r+1)) relative to infinity has order 12
// on the smooth M(12) open.  A point-contact-5 class at u is detected
// exactly by the two standard Taylor covariants E3(u)=E4(u)=0 and the
// square condition f(u)=c^2 != 0.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(8*10^9);

Q := Rationals();
Z := Integers();
QR<r> := PolynomialRing(Q);
S<x> := PolynomialRing(QR);
PQ<X> := PolynomialRing(Q);
PF<Tfrob> := PolynomialRing(Q);

T := -x^2/4-x+r;
h := (x-r)*(T+1);
W16 := 16*(h^2-x^2*T*(T+1));

f := -2*r*x^5
     +(r^2-16*r-4)*x^4
     +(24*r^2-16*r-16)*x^3
     +(-8*r^3+72*r^2+80*r+16)*x^2
     -32*r*(2*r^2+3*r+1)*x
     +16*r^2*(r+1)^2;

assert f eq W16;
print "VINF_IDENTITY_SELF_TEST_PASS", "degree", Degree(f);

function ContactCovariants(g)
    A0 := g;
    A1 := Derivative(g);
    // Scale by rational constants.  Dividing by an element of QR would
    // otherwise promote S=QR[x] to its fraction field.
    A2 := (Q!1/2)*Derivative(g,2);
    A3 := (Q!1/6)*Derivative(g,3);
    A4 := (Q!1/24)*Derivative(g,4);
    Q2 := 4*A0*A2-A1^2;
    E3 := 8*A0^2*A3-A1*Q2;
    E4 := 64*A0^3*A4-Q2^2;
    return E3,E4;
end function;

procedure ContactSelfTest()
    h0 := 2+3*X+5*X^2;
    g0 := h0^2+7*X^5;
    E3,E4 := ContactCovariants(g0);
    assert Evaluate(E3,0) eq 0 and Evaluate(E4,0) eq 0;
    assert Evaluate(g0,0) eq 4;
    print "CONTACT5_PLANTED_SELF_TEST_PASS";
end procedure;

function Specialize(poly,rv)
    return &+[Q!Evaluate(Coefficient(poly,i),rv)*X^i
              : i in [0..Degree(poly)]];
end function;

procedure SmoothFamilySelfTest()
    rv := Q!1;
    fq := Specialize(f,rv);
    assert Degree(fq) eq 5 and Discriminant(fq) ne 0;
    assert Evaluate(fq,0) eq (4*rv*(rv+1))^2;
    C := HyperellipticCurve(fq);
    J := Jacobian(C);
    D12 := J![X,-4*rv*(rv+1)];
    assert Order(D12) eq 12;
    print "VINF_SMOOTH_D12_SELF_TEST_PASS", "r",rv,
          "disc",Discriminant(fq),"D12",D12;
end procedure;

function StripSupport(g,b)
    if IsZero(g) then return g; end if;
    bb := GCD(g,b);
    while Degree(bb) gt 0 do
        g := ExactQuotient(g,bb);
        bb := GCD(g,b);
    end while;
    return g;
end function;

function ContactQuadratic(g,u,c)
    A1 := Evaluate(Derivative(g),u);
    A2 := Evaluate(Derivative(g,2),u)/2;
    dd := A1/(2*c);
    ee := (A2-dd^2)/(2*c);
    hh := c+dd*(X-u)+ee*(X-u)^2;
    return hh;
end function;

function IntegralModel(g)
    L := 1;
    for i in [0..Degree(g)] do
        L := LCM(L,Denominator(Coefficient(g,i)));
    end for;
    return PQ!(L^2*g),L;
end function;

function FrobeniusPolynomial(C,p)
    ef := EulerFactor(C,p);
    dd := Degree(ef);
    return &+[Q!Coefficient(ef,i)*Tfrob^(dd-i) : i in [0..dd]];
end function;

function FullSimplicityCertificate(C)
    plist := [p : p in PrimesUpTo(97) | p ge 3];
    for p in plist do
        try
            Phi := FrobeniusPolynomial(C,p);
            fac := Factorization(Phi);
            if #fac ne 1 or fac[1][2] ne 1 or Degree(fac[1][1]) ne 4 then
                continue;
            end if;
            Gal := GaloisGroup(Phi);
            desc := "";
            try
                desc := TransitiveGroupDescription(Gal);
            catch e2
                desc := "unknown";
            end try;
            if Order(Gal) eq 8 and desc eq "D(4)" then
                return true,"D4",p,Phi;
            end if;
        catch e
            continue;
        end try;
    end for;
    for p in plist do
        try
            Phi := FrobeniusPolynomial(C,p);
            fac := Factorization(Phi);
            if #fac ne 1 or fac[1][2] ne 1 or Degree(fac[1][1]) ne 4 then
                continue;
            end if;
            K<pi> := NumberField(Phi);
            if &and[Degree(MinimalPolynomial(pi^n)) eq 4 : n in [2..12]] then
                return true,"root_power",p,Phi;
            end if;
        catch e
            continue;
        end try;
    end for;
    return false,"none",0,PF!0;
end function;

ContactSelfTest();
SmoothFamilySelfTest();

disc := Discriminant(f);
disc_expected := r^6*(r+1)^4*(r+2)*(r^3-14*r^2-12*r-4);
disc_quot,disc_rem := Quotrem(disc,disc_expected);
assert IsZero(disc_rem) and Degree(disc_quot) eq 0;
print "VINF_DISCRIMINANT", Factorization(disc);

E3,E4 := ContactCovariants(f);
print "CONTACT_COVARIANT_DEGREES",
      "E3_u",Degree(E3),"E3_r",Max([Degree(Coefficient(E3,i)) : i in [0..Degree(E3)]]),
      "E4_u",Degree(E4),"E4_r",Max([Degree(Coefficient(E4,i)) : i in [0..Degree(E4)]]);
lc3 := Coefficient(E3,Degree(E3));
lc4 := Coefficient(E4,Degree(E4));
infinity_support := GCD(lc3,lc4);
print "CONTACT_INFINITY_SUPPORT", Factorization(infinity_support);

print "CONTACT_RESULTANT_BEGIN";
time contact_res := Resultant(E3,E4);
print "CONTACT_RESULTANT_DONE", "degree",Degree(contact_res),
      "terms",#Terms(contact_res);

boundary_support := disc*r*(r+1)*(r+2)*infinity_support;
open_res := StripSupport(contact_res,boundary_support);
print "CONTACT_RESULTANT_OPEN", "degree",Degree(open_res),
      "terms",#Terms(open_res);
print "CONTACT_RESULTANT_FACTORIZATION_BEGIN";
time open_fac := Factorization(open_res);
print "CONTACT_RESULTANT_FACTORIZATION_DONE", "count",#open_fac;
for i in [1..#open_fac] do
    print "OPEN_FACTOR",i,"degree",Degree(open_fac[i][1]),
          "exponent",open_fac[i][2],"poly",open_fac[i][1];
end for;

rational_r := [];
for fe in open_fac do
    if Degree(fe[1]) eq 1 then
        Append(~rational_r,-Coefficient(fe[1],0)/Coefficient(fe[1],1));
    end if;
end for;
rational_r := Setseq({rr : rr in rational_r});
print "RATIONAL_OPEN_R",rational_r;

hits := 0;
for rv in rational_r do
    fq := Specialize(f,rv);
    if Degree(fq) ne 5 or Discriminant(fq) eq 0 then
        print "R_CANDIDATE_BAD",rv;
        continue;
    end if;
    e3q,e4q := ContactCovariants(fq);
    gg := GCD(e3q,e4q);
    print "R_CANDIDATE",rv,"contact_gcd_degree",Degree(gg),
          "contact_gcd",gg;
    for rt in Roots(gg) do
        u := rt[1];
        f0 := Evaluate(fq,u);
        if f0 eq 0 or not IsSquare(f0) then
            print "  U_REJECT",u,"f0",f0,"square",IsSquare(f0);
            continue;
        end if;
        for c in [SquareRoot(f0),-SquareRoot(f0)] do
            hh := ContactQuadratic(fq,u,c);
            if fq-hh^2 ne LeadingCoefficient(fq)*(X-u)^5 then
                continue;
            end if;
            fI,L := IntegralModel(fq);
            C := HyperellipticCurve(fI);
            J := Jacobian(C);
            D12 := J![X,-L*4*rv*(rv+1)];
            D5 := J![X-u,L*c];
            o12 := Order(D12);
            o5 := Order(D5);
            P60 := D12+D5;
            o60 := Order(P60);
            Gtors,phi := TorsionSubgroup(J);
            invs := Invariants(Gtors);
            simple,stype,sp,Phi := FullSimplicityCertificate(C);
            hits +:= 1;
            print "VINF_CONTACT5_HIT", "r",rv,"u",u,"c",c,
                  "orders",<o12,o5,o60>,"torsion",invs,
                  "simple",simple,"simple_type",stype,"simple_prime",sp;
            print "  f_integral",fI;
            print "  D12",D12;
            print "  D5",D5;
            print "  P60",P60;
            if simple then print "  Frobenius",Phi; end if;
        end for;
    end for;
end for;

print "VINF_ELIMINATION_DONE", "rational_r",#rational_r,"hits",hits;
quit;
