//////////////////////////////////////////////////////////////////////
// Exact verification of deep finite-mask survivors in the full cyclic
// contact-5/order-20 plus 3-torsion lane.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(8*10^9);
Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);
PF<Tfrob> := PolynomialRing(Q);

candidates := [
    Q!4356/4199, -Q!4817/8425,
    -Q!22485/10693, Q!75226/34441,
    Q!37925/71877, -Q!69764/79661
];

function FamilyPolynomial(t)
    b := (t^2-1)/2;
    h := 1+t*x+b*x^2;
    return h^2-((t+1)^4/4)*x^5,h,b;
end function;

function IntegralModel(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L,Denominator(Coefficient(f,i)));
    end for;
    return P!(L^2*f),L;
end function;

function FrobeniusPolynomial(C,p)
    ef := EulerFactor(C,p);
    d := Degree(ef);
    return &+[Q!Coefficient(ef,i)*Tfrob^(d-i) : i in [0..d]];
end function;

function FullSimplicityCertificate(C)
    plist := [p : p in PrimesUpTo(151) | p ge 3];
    for p in plist do
        try
            Phi := FrobeniusPolynomial(C,p);
            fac := Factorization(Phi);
            if #fac ne 1 or fac[1][2] ne 1 or Degree(fac[1][1]) ne 4 then
                continue;
            end if;
            Gal := GaloisGroup(Phi);
            desc := "unknown";
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

print "CONTACT5 ORDER20+3 DEEP CANDIDATE VERIFICATION";
for t in candidates do
    f,h,b := FamilyPolynomial(t);
    assert Degree(f) eq 5 and Discriminant(f) ne 0;
    fI,L := IntegralModel(f);
    C := HyperellipticCurve(fI);
    J := Jacobian(C);

    D5 := J![x,Q!L];
    H4 := J![x^2+(2/(t+1))*x,L*((t+2)*x+1)];
    assert Order(D5) eq 5 and Order(H4) eq 4;
    D20 := D5+H4;
    assert Order(D20) eq 20;

    first_kill := 0;
    gcd_orders := 0;
    used := [];
    for p in [ell : ell in PrimesUpTo(251) | ell notin {2,3,5}] do
        try
            fp := ChangeRing(fI,GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then continue; end if;
            Np := Z!#Jacobian(HyperellipticCurve(fp));
            gcd_orders := gcd_orders eq 0 select Np else GCD(gcd_orders,Np);
            Append(~used,<p,Np,Np mod 3,gcd_orders>);
            if first_kill eq 0 and Np mod 3 ne 0 then first_kill := p; end if;
        catch e
            continue;
        end try;
    end for;
    print "CANDIDATE","t",t,"b",b,
          "good_primes",#used,"first_3_kill",first_kill,
          "point_count_gcd",gcd_orders;

    G,phi := TorsionSubgroup(J);
    invs := Invariants(G);
    print " EXACT_TORSION",invs,"D20_order",Order(D20);
    if &or[n mod 3 eq 0 : n in invs] then
        gens := [phi(G.i) : i in [1..Ngens(G)]];
        print " TORSION_GENERATORS",gens,
              "generator_orders",[Order(g) : g in gens];
    end if;
    simple,stype,sp,Phi := FullSimplicityCertificate(C);
    print " SIMPLICITY",simple,stype,sp;
    if simple then print " FROBENIUS",Phi; end if;
    if invs eq [60] and simple then
        print "CYCLIC_Z60_CERTIFIED","t",t;
        print " f_integral",fI;
        print " D20",D20;
    end if;
end for;
print "CANDIDATE_VERIFICATION_DONE";
quit;
