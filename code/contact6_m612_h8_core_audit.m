//////////////////////////////////////////////////////////////////////
// Exact audit of the two curves underlying the four p=5,p=7-boundary
// core points found by the height-8 optimized T0 slice search.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
P<x> := PolynomialRing(Q);

function IntegralModel(f)
    d := LCM([Denominator(Coefficient(f,i)) : i in [0..Degree(f)]]);
    return P!(d^2*f),d;
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

recs := [
    <"omega=+-4",Q!5,Q!(-5/2),Q!10,Q!(-14),Q!(-2)>,
    <"omega=+-1/4",Q!(-5/2),Q!5,Q!(5/4),Q!(-7/2),Q!(-1/2)>
];

for rec in recs do
    label,a,b,L,U,nu := Explode(rec);
    h := 1+a*x+b*x^2+x^3;
    f := h^2-(x-1)^6;
    M := L^2;
    c5 := 2*b+6;
    c4 := 2*a+b^2-15;
    B3 := c5*M+3*U;
    Delta3 := 4*c4*M+12*(U^2+nu^2)-B3^2;
    q3 := x^2+U*x+nu^2;
    h3 := (1/L)*x^3+(B3/(2*L))*x^2+(Delta3/(8*L))*x+nu^3/L;
    print "CORE",label,"a",a,"b",b,"L",L,"U",U,"nu",nu;
    print "factorization",Factorization(f),"discriminant",Discriminant(f);
    print "contact_identity",h3^2-f eq (1/M)*q3^3;
    if Discriminant(f) eq 0 then
        print "SINGULAR_SKIP";
        continue;
    end if;
    J := Jacobian(HyperellipticCurve(f));
    D := J![x-1,Evaluate(h,1)];
    E := J![q3,h3 mod q3];
    print "orders D,E",Order(D),Order(E);
    fI,d := IntegralModel(f);
    G,mp := TorsionSubgroup(Jacobian(HyperellipticCurve(fI)));
    simple,p,Lp := SimpleCertificate(fI);
    print "torsion",Invariants(G),"simple",simple,"p",p,"Lp",Lp;
    print "integral_curve",fI;
end for;

quit;
