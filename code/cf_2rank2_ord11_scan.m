// [2,22] hunt, route 2: D_inf of order 11 on the built-in 2-rank-2 family
//   f = (x^2 - 1)(x^2 + a x + b)(x^2 + c x + d),  a,b,c,d integers
// (factor type [1,1,2,2]/[2,2,2]-flavored: 2-rank 2 pointwise, cf.
// two-rank-and-factor-types; the family of code/agent_a2_24_cf_search.m).
// A member with D_inf order divisible by 11 has J(Q) >= (Z/2)^2 x Z/11
// = [2,22] -- no halving needed.  Any such member found is then checked
// with exact TorsionSubgroup, and (separately) certified simple + End=Z.
//
// Funnel (same 3-stage design as code/cf_cyclic_gaps_scan.m):
//   mod-p1 CF order in TARGETS u CALIB -> mod-p2 agree -> exact CF over Q.
// Self-test vectors run first.  Calibration inventory: all odd exact orders
// >= 5 (if the box shows no order 5/7/9, it cannot be trusted about 11).
//
// Symmetries quotiented: (a,b) <-> (c,d) swap and x -> -x (a -> -a, c -> -c):
// enumerate a >= 0 and (a,b) <= (c,d) lexicographically.
//
// Run: magma -b H:=12 code/cf_2rank2_ord11_scan.m

SetColumns(0);
SetSeed(1);
if not assigned H then H := 12; elif Type(H) eq MonStgElt then H := StringToInteger(H); end if;
if not assigned NParts then NParts := 1; elif Type(NParts) eq MonStgElt then NParts := StringToInteger(NParts); end if;
if not assigned Part then Part := 0; elif Type(Part) eq MonStgElt then Part := StringToInteger(Part); end if;
if not assigned MemGB then MemGB := 3; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
if not assigned MaxSteps then MaxSteps := 60; elif Type(MaxSteps) eq MonStgElt then MaxSteps := StringToInteger(MaxSteps); end if;
if not assigned MaxOrd then MaxOrd := 45; elif Type(MaxOrd) eq MonStgElt then MaxOrd := StringToInteger(MaxOrd); end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();
PQ<x> := PolynomialRing(Q);

function SqrtPolyPart(f)
    P := Parent(f); xx := P.1;
    s := xx^3;
    for k in [1..3] do
        d := f - s^2;
        if Degree(d) le 2 then break; end if;
        s := s + (Coefficient(d, 6-k)/(2*Coefficient(s, 3)))*xx^(3-k);
    end for;
    return s;
end function;

function CFOrderB(f, maxsteps, maxord)
    P := Parent(f);
    s := SqrtPolyPart(f);
    Pi := P!0; Qi := P!1; total := 0;
    for i in [0..maxsteps] do
        if Qi eq 0 then return 0; end if;
        ai := (Pi + s) div Qi;
        total +:= Degree(ai);
        if total gt maxord then return 0; end if;
        Pn := ai*Qi - Pi;
        if (f - Pn^2) mod Qi ne 0 then return 0; end if;
        Qn := (f - Pn^2) div Qi;
        Pi := Pn; Qi := Qn;
        if i ge 1 and Degree(Qi) le 0 and Qi ne 0 then return total; end if;
    end for;
    return 0;
end function;

// SELF-TEST
f14 := (x^2+1)*(x^4+5*x^2+4*x+4);
f18 := (x^2-x+1)*(x^4-x^3+9*x^2+8*x-8);
f28t := x^6+2*x^5-5*x^4-14*x^3-3*x^2+24*x+28;
assert CFOrderB(f14, MaxSteps, MaxOrd) eq 14;
assert CFOrderB(f18, MaxSteps, MaxOrd) eq 18;
assert CFOrderB(f28t, MaxSteps, MaxOrd) eq 7;
for p in [101, 103] do
    Pp := PolynomialRing(GF(p));
    assert CFOrderB(Pp!f14, MaxSteps, MaxOrd) eq 14;
end for;
printf "SELFTEST_PASS\n";

primes := [101, 103, 107, 109];
Pps := [PolynomialRing(GF(p)) : p in primes];

// interesting exact orders: odd >= 5 (calibration) and multiples of 11 (hits)
function Interesting(o)
    return (o ge 5 and IsOdd(o)) or (o gt 0 and o mod 11 eq 0);
end function;

counts := AssociativeArray();
tested := 0; s2pass := 0; exactrun := 0; nearmiss := 0; idx := -1;

for a in [0..H] do
 for b in [-H..H] do
  for c in [-H..H] do
   for d in [-H..H] do
    if [a,b] gt [c,d] then continue; end if;   // swap symmetry
    idx +:= 1;
    if idx mod NParts ne Part then continue; end if;
    tested +:= 1;
    if tested mod 200000 eq 0 then
        printf "PROGRESS tested=%o s2pass=%o exact=%o\n", tested, s2pass, exactrun;
    end if;
    // f = (x^2-1)(x^2+ax+b)(x^2+cx+d): coefficients expanded once, cheaply
    // via mod-p polynomial products in stage 1.
    o1 := -1; i1 := 0;
    for i in [1..#primes] do
        Pp := Pps[i]; xp := Pp.1;
        fb := (xp^2-1)*(xp^2+a*xp+b)*(xp^2+c*xp+d);
        if Degree(GCD(fb, Derivative(fb))) eq 0 then
            o1 := CFOrderB(fb, MaxSteps, MaxOrd); i1 := i; break;
        end if;
    end for;
    if o1 le 0 or not Interesting(o1) then continue; end if;
    o2 := -1;
    for i in [i1+1..#primes] do
        Pp := Pps[i]; xp := Pp.1;
        fb := (xp^2-1)*(xp^2+a*xp+b)*(xp^2+c*xp+d);
        if Degree(GCD(fb, Derivative(fb))) eq 0 then
            o2 := CFOrderB(fb, MaxSteps, MaxOrd); break;
        end if;
    end for;
    if o2 ne o1 then continue; end if;
    s2pass +:= 1;
    f := (x^2-1)*(x^2+a*x+b)*(x^2+c*x+d);
    if Degree(GCD(f, Derivative(f))) ne 0 then continue; end if;
    exactrun +:= 1;
    oQ := CFOrderB(f, MaxSteps, MaxOrd);
    if oQ eq 0 then nearmiss +:= 1; continue; end if;
    if not Interesting(oQ) then continue; end if;
    if IsDefined(counts, oQ) then counts[oQ] +:= 1; else counts[oQ] := 1; end if;
    printf "CALIB ord=%o a=%o b=%o c=%o d=%o\n", oQ, a, b, c, d;
    if oQ mod 11 eq 0 then
        printf "HIT11 ord=%o a=%o b=%o c=%o d=%o f=%o\n", oQ, a, b, c, d, f;
    end if;
   end for;
  end for;
 end for;
end for;

printf "CALIBRATION_TABLE (H=%o part %o/%o):\n", H, Part, NParts;
for o in Sort(Setseq(Keys(counts))) do
    printf "  ORDER %o : %o\n", o, counts[o];
end for;
printf "SEARCH_DONE H=%o Part=%o/%o tested=%o s2pass=%o exact=%o nearmiss=%o\n",
       H, Part, NParts, tested, s2pass, exactrun, nearmiss;
quit;
