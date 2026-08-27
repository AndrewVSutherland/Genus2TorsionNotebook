//////////////////////////////////////////////////////////////////////
//  Finite-field smoke test for the contact-5/contact-5 obstruction.
//
//  Over Q the normalized route is blocked because
//
//      Phi = x^5 - (x-1)^5
//
//  is irreducible.  Over finite fields where Phi has a quadratic
//  divisor, the same sum/difference construction should produce curves
//  with two independent rational 5-contact classes.  This script checks
//  that sanity test without writing any output files.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned prime_bound then
    prime_bound := 50;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;

Z := Integers();

function RelationPairs(J, D0, D1)
    rels := [];
    for a in [0..4] do
        for b in [0..4] do
            if a eq 0 and b eq 0 then
                continue;
            end if;
            if a*D0 + b*D1 eq J!0 then
                Append(~rels, <a,b>);
            end if;
        end for;
    end for;
    return rels;
end function;

function DegreeTwoDivisors(Phi, fac)
    divisors := [];
    n := #fac;
    for mask in [0..2^n-1] do
        U := Parent(Phi)!1;
        for i in [1..n] do
            if ((mask div 2^(i-1)) mod 2) eq 1 then
                U *:= fac[i][1];
            end if;
        end for;
        if Degree(U) eq 2 then
            Append(~divisors, U);
        end if;
    end for;
    return divisors;
end function;

print "# finite-field contact-5/contact-5 smoke test";
print "# Looking for p where Phi has a degree-2 divisor.";

for p in PrimesUpTo(prime_bound) do
    if p in {2,5} then
        continue;
    end if;

    F := GF(p);
    P<x> := PolynomialRing(F);
    Phi := x^5 - (x-1)^5;
    fac := Factorization(Phi);
    degs := Sort([Degree(pair[1]) : pair in fac]);
    divisors := DegreeTwoDivisors(Phi, fac);

    printf "p=%o factor_degrees=%o degree2_divisors=%o", p, degs, #divisors;
    if #divisors eq 0 then
        print "";
        continue;
    end if;

    smooth_seen := false;
    independent_seen := false;
    for U in divisors do
        V := ExactQuotient(Phi, U);
        if Degree(V) gt 2 then
            continue;
        end if;

        h0 := (U + V)/2;
        h1 := (V - U)/2;
        f := h0^2 - x^5;
        if Degree(f) ne 5 or Discriminant(f) eq 0 then
            continue;
        end if;

        smooth_seen := true;
        C := HyperellipticCurve(f);
        J := Jacobian(C);
        D0 := J![x, Evaluate(h0, 0)];
        D1 := J![x-1, Evaluate(h1, 1)];
        rels := RelationPairs(J, D0, D1);
        ord0 := Order(D0);
        ord1 := Order(D1);
        ordsum := Order(D0 + D1);

        printf "\n  U=%o\n  V=%o\n", U, V;
        printf "  h0=%o\n  h1=%o\n  f=%o\n", h0, h1, f;
        printf "  #J=%o orders=(%o,%o,%o) relations=%o\n",
                #J, ord0, ord1, ordsum, rels;

        if ord0 eq 5 and ord1 eq 5 and #rels eq 0 then
            independent_seen := true;
            print "  INDEPENDENT_F5_PAIR";
            break;
        end if;
    end for;

    if not smooth_seen then
        print " no_smooth_model_from_degree2_divisor";
    elif not independent_seen then
        print " no_independent_pair_seen";
    end if;
end for;

quit;
