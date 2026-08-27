//////////////////////////////////////////////////////////////////////
//  Low-memory marked 7-division sieve on the contact-7 family.
//
//      h = 1 - (7/2)*x + a*x^2 + b*x^3,
//      f = (h^2 + (x - 1)^7)/x^2,
//      D7 = [(x-1), h(1)].
//
//  For an exact cyclic [49] specialization, the nonzero marked class
//  D7 must lie in 7*J(Q).  At every good prime p != 7 this implies
//  D7_bar lies in 7*J(F_p).  The finite test below uses full abelian
//  coordinates of D7_bar.  This is stronger than either 49 | #J(F_p)
//  or 49 | exponent(J(F_p)).
//
//  Rational parameter pairs are streamed.  Only O(#parameters*#primes)
//  integer residues and O(sum p^2) Boolean masks are stored; no sequence
//  of all O(#parameters^2) pairs is built.
//
//  Typical run:
//      magma -b height:=50 prime_bound:=43 \
//          code/z49_local_marked_h50_sieve.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned height then
    height := 50;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned prime_bound then
    prime_bound := 43;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;
if not assigned max_record then
    max_record := 12;
elif Type(max_record) eq MonStgElt then
    max_record := StringToInteger(max_record);
end if;
if not assigned progress_interval then
    progress_interval := 250000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;

Z := Integers();
Q := Rationals();
P<x> := PolynomialRing(Q);
filter_primes := [ p : p in PrimesUpTo(prime_bound) | p ge 3 and p ne 7 ];

function HasExponent49(invs)
    return exists{ n : n in invs | (Z!n) mod 49 eq 0 };
end function;

function MarkedDivisibleBy7(D, G, phi)
    g := D @@ phi;
    coords := Eltseq(g);
    invs := Invariants(G);
    for i in [1..#coords] do
        if (Z!coords[i]) mod GCD(7, Z!invs[i]) ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;

// Sanity checks for the coordinate criterion n*q=g in Z/mZ.
function CoordinateDivisible(c, m, n)
    return (Z!c) mod GCD(Z!m, Z!n) eq 0;
end function;
assert not CoordinateDivisible(1, 7, 7);
assert CoordinateDivisible(7, 49, 7);
assert CoordinateDivisible(14, 98, 7);

function Contact7PolynomialFinite(F, a, b)
    PF<X> := PolynomialRing(F);
    h := 1 - (F!7/F!2)*X + a*X^2 + b*X^3;
    num := h^2 + (X-1)^7;
    if Coefficient(num,0) ne 0 or Coefficient(num,1) ne 0 then
        return false, PF!0, PF!0;
    end if;
    f := ExactQuotient(num, X^2);
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        return false, f, h;
    end if;
    if Evaluate(h,F!1) eq 0 then
        return false, f, h;
    end if;
    return true, f, h;
end function;

function Contact7Polynomial(a, b)
    h := 1 - (Q!7/2)*x + a*x^2 + b*x^3;
    num := h^2 + (x-1)^7;
    if Coefficient(num,0) ne 0 or Coefficient(num,1) ne 0 then
        return false, P!0, P!0;
    end if;
    return true, ExactQuotient(num,x^2), h;
end function;

function RationalParametersOfHeight(B)
    vals := [];
    seen := {};
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num,den) ne 1 then
                continue;
            end if;
            q := Q!num/Q!den;
            key := Sprint(q);
            if key notin seen then
                Include(~seen,key);
                Append(~vals,q);
            end if;
        end for;
    end for;
    return vals;
end function;

// The second irreducible discriminant factor in the normalization used
// here.  The full identity is
//   256*Disc(f) = (2*a+2*b-5)^7 * Q5(a,b).
function Q5Component(a,b)
    return 432*a^4 - 64*a^3*b^3 + 1008*a^3*b + 3024*a^3
        - 448*a^2*b^3 + 224*a^2*b^2 + 21168*a^2*b - 32536*a^2
        - 2016*a*b^4 + 4480*a*b^3 + 38416*a*b^2 - 109760*a*b
        + 78890*a - 864*b^5 + 5936*b^4 + 7056*b^3
        - 96040*b^2 + 120050*b - 60025;
end function;

// true in the kill mask means: good marked fiber, but D7 is not in 7J.
function BuildMarkedKillMask(p)
    F := GF(p);
    PF<X> := PolynomialRing(F);
    mask := [ false : i in [1..p^2] ];
    total := p^2;
    good := 0;
    boundary := 0;
    order49 := 0;
    exponent49 := 0;
    marked_div7 := 0;
    exponent_not_marked := 0;

    for ai in [0..p-1] do
        for bi in [0..p-1] do
            idx := ai*p + bi + 1;
            ok,f,h := Contact7PolynomialFinite(F,F!ai,F!bi);
            if not ok then
                boundary +:= 1;
                continue;
            end if;
            good +:= 1;
            J := Jacobian(HyperellipticCurve(f));
            D7 := J![X-1,Evaluate(h,F!1)];
            if Order(D7) ne 7 then
                // Retain an unexpected marked-boundary fiber rather than
                // using it as an obstruction.
                boundary +:= 1;
                good -:= 1;
                continue;
            end if;
            G,phi := AbelianGroup(J);
            invs := [Z!n : n in Invariants(G)];
            if (#G) mod 49 eq 0 then
                order49 +:= 1;
            end if;
            exp49 := HasExponent49(invs);
            if exp49 then
                exponent49 +:= 1;
            end if;
            div7 := MarkedDivisibleBy7(D7,G,phi);
            if div7 then
                marked_div7 +:= 1;
            else
                mask[idx] := true;
                if exp49 then
                    exponent_not_marked +:= 1;
                end if;
            end if;
        end for;
    end for;

    return mask,
        <total,good,boundary,order49,exponent49,marked_div7,
         exponent_not_marked>;
end function;

function ParameterResidues(params,p)
    residues := [];
    F := GF(p);
    for q in params do
        d := Denominator(q);
        if d mod p eq 0 then
            Append(~residues,-1);
        else
            r := F!(Numerator(q) mod p)/F!(d mod p);
            Append(~residues,Z!r);
        end if;
    end for;
    return residues;
end function;

print "Z49 MARKED-DIVISION LOW-MEMORY SIEVE";
print "height",height,"prime_bound",prime_bound,"primes",filter_primes;
print "criterion","D7_bar lies in 7*J(F_p) at every good p != 7";

masks := [];
for p in filter_primes do
    mask,summary := BuildMarkedKillMask(p);
    Append(~masks,mask);
    printf "LOCAL p=%o total=%o good=%o boundary=%o order49=%o exponent49=%o marked_div7=%o exponent_not_marked=%o\n",
        p,summary[1],summary[2],summary[3],summary[4],summary[5],
        summary[6],summary[7];
end for;

params := RationalParametersOfHeight(height);
residue_tables := [ ParameterResidues(params,p) : p in filter_primes ];
nparams := #params;
checked := 0;
mask_survivors := 0;
smooth_survivors := 0;
singular_survivors := 0;
singular_h_component := 0;
singular_q5_component := 0;
singular_both := 0;
singular_unclassified := 0;
marked_failures := 0;
first_kill := [0 : j in [1..#filter_primes]];
denominator_unresolved := [0 : j in [1..#filter_primes]];
smooth_records := [];
unclassified_records := [];

for ia in [1..nparams] do
    a := params[ia];
    for ib in [1..nparams] do
        checked +:= 1;
        if progress_interval gt 0 and checked mod progress_interval eq 0 then
            printf "PROGRESS checked=%o mask_survivors=%o smooth_survivors=%o\n",
                checked,mask_survivors,smooth_survivors;
        end if;

        killed := false;
        for j in [1..#filter_primes] do
            ai := residue_tables[j][ia];
            bi := residue_tables[j][ib];
            if ai lt 0 or bi lt 0 then
                denominator_unresolved[j] +:= 1;
                continue;
            end if;
            p := filter_primes[j];
            if masks[j][ai*p+bi+1] then
                first_kill[j] +:= 1;
                killed := true;
                break;
            end if;
        end for;
        if killed then
            continue;
        end if;

        mask_survivors +:= 1;
        b := params[ib];
        ok,f,h := Contact7Polynomial(a,b);
        if not ok or Degree(f) ne 5 or Discriminant(f) eq 0 then
            singular_survivors +:= 1;
            hz := 2*a+2*b-5 eq 0;
            qz := Q5Component(a,b) eq 0;
            if hz and qz then
                singular_both +:= 1;
            elif hz then
                singular_h_component +:= 1;
            elif qz then
                singular_q5_component +:= 1;
            else
                singular_unclassified +:= 1;
                if #unclassified_records lt max_record then
                    Append(~unclassified_records,<a,b>);
                end if;
            end if;
            continue;
        end if;

        yP := Evaluate(h,Q!1);
        if yP eq 0 then
            marked_failures +:= 1;
            continue;
        end if;
        J := Jacobian(HyperellipticCurve(f));
        D7 := J![x-1,yP];
        if Order(D7) ne 7 then
            marked_failures +:= 1;
            continue;
        end if;
        smooth_survivors +:= 1;
        if #smooth_records lt max_record then
            Append(~smooth_records,<a,b>);
        end if;
    end for;
end for;

print "FIRST_KILL_COUNTS",[<filter_primes[j],first_kill[j]> :
    j in [1..#filter_primes]];
print "DENOMINATOR_UNRESOLVED_VISITS",[<filter_primes[j],denominator_unresolved[j]> :
    j in [1..#filter_primes]];
printf "DONE parameters=%o checked=%o mask_survivors=%o smooth_survivors=%o marked_failures=%o\n",
    nparams,checked,mask_survivors,smooth_survivors,marked_failures;
printf "SINGULAR total=%o h_component=%o q5_component=%o both=%o unclassified=%o\n",
    singular_survivors,singular_h_component,singular_q5_component,
    singular_both,singular_unclassified;
print "SMOOTH_RECORDS",smooth_records;
print "UNCLASSIFIED_SINGULAR_RECORDS",unclassified_records;

quit;
