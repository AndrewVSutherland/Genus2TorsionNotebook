//////////////////////////////////////////////////////////////////////
//  Probe whether the contact-6 extra-root family is generically simple.
//  This does not impose [2,6,6]; it only tests if forcing an extra rational
//  root structurally forces decomposability.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Q := Rationals();
P<x> := PolynomialRing(Q);

function ExtraRootA(eps, r, b)
    return (eps*(r-1)^3 - 1 - r^3 - b*r^2)/r;
end function;

function Contact6Polynomial(a, b)
    h := 1 + a*x + b*x^2 + x^3;
    return h^2 - (x-1)^6, h;
end function;

function SimpleCertificate(f)
    C := HyperellipticCurve(f);
    for p in [5,7,11,13,17,19,23,29,31,37,41,43,47] do
        try
            fp := ChangeRing(f, GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then
                continue;
            end if;
            Lp := LPolynomial(ChangeRing(C, GF(p)));
            fac := Factorization(Lp);
            if #fac eq 1 and fac[1][2] eq 1 and Degree(fac[1][1]) eq 4 then
                return true, p, Lp;
            end if;
        catch e
            continue;
        end try;
    end for;
    return false, 0, P!0;
end function;

samples := [
    <1, Q!2, Q!0>,
    <1, Q!2, Q!1>,
    <1, Q!(3/2), Q!(-1/2)>,
    <-1, Q!2, Q!1>,
    <1, Q!(5/3), Q!(2/5)>,
    <-1, Q!(-2), Q!(1/3)>
];

print "Generic simplicity probe for contact-6 extra-root family";
for S in samples do
    eps := S[1]; r := S[2]; b := S[3];
    a := ExtraRootA(eps,r,b);
    f,h := Contact6Polynomial(a,b);
    if Degree(f) ne 5 or Discriminant(f) eq 0 or Evaluate(h,Q!1) eq 0 then
        print "sample", S, "bad";
        continue;
    end if;
    simple, pcert, Lp := SimpleCertificate(f);
    print "sample", S, "a", a, "factor_degrees",
          Sort([Degree(fe[1]) : fe in Factorization(f)]),
          "simple_cert", simple, "p", pcert, "Lp", Lp;
end for;

quit;
