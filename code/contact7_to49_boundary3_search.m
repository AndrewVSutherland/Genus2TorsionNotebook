//////////////////////////////////////////////////////////////////////
//  Z/49 lane: first blow-up probe for the degenerate h(1)=0
//  p=3 boundary branches of the contact-7 family.
//
//  This is deliberately local: it does not repeat the broad
//  cyclic-49 open search height search.  It uses the chart
//
//      a = a0 + 3*u,
//      b = 5/2 - a + 3^7*s,
//
//  where a0=0 gives the residue branch (a,b)=(0,1), and a0=2 gives
//  (a,b)=(2,2).  Then h(1)=3^7*s.  The first blow-up layer is s a
//  3-adic unit.
//
//  Why exponent 7 is forced: the one-node component group has order
//  n=7*v_3(h(1)).  The marked contact class has order 7 in Z/nZ, and
//  it can lie in 7*(Z/nZ) only when 49 divides n.  Hence the first
//  possible layer is v_3(h(1))=7.
//
//
//  Typical runs:
//      magma -b branch_a0:=0 height:=8 first_layer_only:=1 \
//          code/contact7_to49_boundary3_search.m
//      magma -b branch_a0:=2 height:=8 first_layer_only:=1 \
//          code/contact7_to49_boundary3_search.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned branch_a0 then
    branch_a0 := 0;
elif Type(branch_a0) eq MonStgElt then
    branch_a0 := StringToInteger(branch_a0);
end if;
if not assigned height then
    height := 8;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned first_layer_only then
    first_layer_only := 1;
elif Type(first_layer_only) eq MonStgElt then
    first_layer_only := StringToInteger(first_layer_only);
end if;
if not assigned max_layer then
    max_layer := 1;
elif Type(max_layer) eq MonStgElt then
    max_layer := StringToInteger(max_layer);
end if;
if not assigned progress_interval then
    progress_interval := 20000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;
if not assigned max_contact_checks then
    max_contact_checks := 5;
elif Type(max_contact_checks) eq MonStgElt then
    max_contact_checks := StringToInteger(max_contact_checks);
end if;

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

filter_primes := [5,11,13,17,19,23,29,31,37,41,43];

function RationalParametersOfHeight(B)
    vals := [];
    seen := {};
    for den in [1..B] do
        if Valuation(den, 3) gt 0 then
            continue;
        end if;
        for num in [-B..B] do
            if GCD(num, den) ne 1 then
                continue;
            end if;
            q := Q!num/Q!den;
            key := Sprint(q);
            if key notin seen then
                Include(~seen, key);
                Append(~vals, q);
            end if;
        end for;
    end for;
    return vals;
end function;

function V3(q)
    return Valuation(Numerator(q), 3) - Valuation(Denominator(q), 3);
end function;

function Q5Component(a, b)
    return 432*a^4 - 64*a^3*b^3 + 1008*a^3*b + 3024*a^3
        - 448*a^2*b^3 + 224*a^2*b^2 + 21168*a^2*b - 32536*a^2
        - 2016*a*b^4 + 4480*a*b^3 + 38416*a*b^2 - 109760*a*b
        + 78890*a - 864*b^5 + 5936*b^4 + 7056*b^3
        - 96040*b^2 + 120050*b - 60025;
end function;

function LocalAB(a0, u, s)
    a := Q!a0 + 3*u;
    b := Q!5/2 - a + 3^7*s;
    return a, b;
end function;

function Contact7Polynomial(a, b)
    h := 1 - (Q!7/2)*x + a*x^2 + b*x^3;
    return ExactQuotient(h^2 + (x - 1)^7, x^2), h;
end function;

function GoodHyperellipticPolynomial(f)
    return Degree(f) eq 5 and Discriminant(f) ne 0;
end function;

function Passes49AwayFrom3(f,h)
    used := [];
    for p in filter_primes do
        try
            F := GF(p);
            fp := ChangeRing(f,F);
            hp := ChangeRing(h,F);
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 or Evaluate(hp,F!1) eq 0 then
                continue;
            end if;
            Jp := Jacobian(HyperellipticCurve(fp));
            xp := Parent(fp).1;
            Dp := Jp![xp-1,Evaluate(hp,F!1)];
            Gp,phi := AbelianGroup(Jp);
            gp := Dp @@ phi;
            invs := Invariants(Gp);
            coords := Eltseq(gp);
            div7 := &and[(Z!coords[i]) mod GCD(7,Z!invs[i]) eq 0 : i in [1..#coords]];
            Append(~used,<p,#Jp,invs,coords>);
            if not div7 then
                return false,p,#Jp,used;
            end if;
        catch e
            continue;
        end try;
    end for;
    return true,0,0,used;
end function;

function ContactClassCheck(f, h)
    h1 := Evaluate(h, Q!1);
    if h1 eq 0 then
        return false, false, "h(1)=0";
    end if;
    try
        C := HyperellipticCurve(f);
        J := Jacobian(C);
        D := J![x - 1, P!h1];
        seven_zero := IsZero(7*D);
        exact_order_7 := seven_zero and &and [ not IsZero(i*D) : i in [1..6] ];
        return seven_zero, exact_order_7, "";
    catch e
        return false, false, e`Object;
    end try;
end function;

function Exact49Check(f,h)
    L := LCM([Denominator(Coefficient(f,i)) : i in [0..Degree(f)]]);
    fI := P!(L^2*f);
    JI := Jacobian(HyperellipticCurve(fI));
    D7I := JI![x-1,L*Evaluate(h,Q!1)];
    ok,Q49 := IsDivisibleBy(D7I,7);
    return ok,Q49,Order(D7I),fI;
end function;

procedure Increment(~A, key)
    if not IsDefined(A, key) then
        A[key] := 0;
    end if;
    A[key] +:= 1;
end procedure;

if branch_a0 notin {0,2} then
    error "branch_a0 must be 0 or 2";
end if;

b0 := (branch_a0 eq 0) select 1 else 2;
q5mod3 := GF(3)!Q5Component(Q!branch_a0, Q!b0);

print "Z49 contact-7 p=3 h(1)=0 blow-up probe";
print "branch residue", <branch_a0,b0>;
print "chart", "a=a0+3*u, b=5/2-a+3^7*s, h(1)=3^7*s";
print "first_layer_only", first_layer_only, "max_layer", max_layer;
print "height", height, "filter_primes", filter_primes;
print "Q5(a0,b0) mod 3", q5mod3;
print "h'(1) on chart", "4-a0-3*u+3^8*s; a 3-adic unit on both branches";
print "disc visible factor", "(2a+2b-5)^7=2^7*3^49*s^7";

params := RationalParametersOfHeight(height);
print "z3-integral chart parameters", #params;

checked := 0;
smooth := 0;
layer_counts := AssociativeArray();
first_kill := AssociativeArray();
survivors := 0;
contact_checks := 0;
contact_failures := 0;
exact_tests := 0;
hits49 := 0;

for u in params do
    for s in params do
        if s eq 0 then
            continue;
        end if;
        layer := V3(s);
        if first_layer_only ne 0 and layer ne 0 then
            continue;
        end if;
        if first_layer_only eq 0 and layer gt max_layer then
            continue;
        end if;

        checked +:= 1;
        if progress_interval gt 0 and checked mod progress_interval eq 0 then
            print "progress", "checked", checked, "smooth", smooth,
                  "survivors", survivors, "contact_checks", contact_checks;
        end if;

        a, b := LocalAB(branch_a0, u, s);
        f, h := Contact7Polynomial(a, b);
        if not GoodHyperellipticPolynomial(f) then
            continue;
        end if;
        smooth +:= 1;
        Increment(~layer_counts, layer);

        if contact_checks lt max_contact_checks then
            seven_zero, exact_order_7, msg := ContactClassCheck(f, h);
            contact_checks +:= 1;
            if not exact_order_7 then
                contact_failures +:= 1;
            end if;
            print "CONTACT_SAMPLE", "u", u, "s", s, "a", a, "b", b,
                  "v3_s", layer, "7D=0", seven_zero,
                  "exact_order_7", exact_order_7, "msg", msg;
        end if;

        pass, pbad, nbad, used := Passes49AwayFrom3(f,h);
        if not pass then
            Increment(~first_kill, pbad);
            continue;
        end if;

        survivors +:= 1;
        print "DIV7_SURVIVOR", "u", u, "s", s, "a", a, "b", b,
              "v3_s", layer, "used", used;
        try
            ok49,Q49,ord7,fI := Exact49Check(f,h);
            exact_tests +:= 1;
            print "EXACT_DIV7",ok49,"D7order",ord7;
            if ok49 then
                hits49 +:= 1;
                print "HIT49","Q",Q49,"order",Order(Q49),"f",fI;
            end if;
        catch e
            print "EXACT_ERROR",e`Object;
        end try;
    end for;
end for;

print "DONE";
print "branch", <branch_a0,b0>;
print "checked", checked;
print "smooth", smooth;
print "survivors", survivors;
print "contact_checks", contact_checks;
print "contact_failures", contact_failures;
print "exact_tests",exact_tests;
print "hits49",hits49;

print "LAYER_COUNTS";
for key in Sort([ k : k in Keys(layer_counts) ]) do
    print key, layer_counts[key];
end for;

print "FIRST_KILL_DIV7_AWAY_FROM_3";
for key in Sort([ k : k in Keys(first_kill) ]) do
    print key, first_kill[key];
end for;

quit;
