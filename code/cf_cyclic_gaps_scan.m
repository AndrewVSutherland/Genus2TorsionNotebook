// Lane 2B: CF/Pell exact-order scan for the open cyclic gap orders
// {31, 35, 37, 38} (and calibration inventory 25..40) over raw boxes of
// monic squarefree sextics  f = x^6 + c5 x^5 + ... + c0,  |ci| <= H,
// c5 in [0..3] (translation x->x+e shifts c5 by 6e and x->-x negates it,
// so c5 in [0..3] covers all orbits; other symmetries may double-count,
// which only wastes time, never misses).
//
// D_inf = inf+ - inf- is torsion of exact order N iff the polynomial CF of
// sqrt(f) is quasi-periodic with degree total N (incl. deg a_0 = 3) --
// see .claude/skills/pell-cf-order.  Funnel:
//   stage 1: CF order mod p1 (exact in J(F_p1); N-torsion over Q reduces to
//            exact order N at every good prime, so no true hit is lost);
//   stage 2: same at p2, must agree;
//   stage 3: exact CF order over Q.
// A CFHIT of order N IS a rational point of exact order N on J -- the
// endgame (TorsionSubgroup + simplicity cert) runs separately per hit.
//
// SELF-TEST runs first (f14->14, f18->18, f28-trap->7, over Q and mod p);
// the scan aborts if any vector fails.
//
// Run: magma -b H:=3 NParts:=1 Part:=0 code/cf_cyclic_gaps_scan.m

SetColumns(0);
SetSeed(1);

if not assigned H then H := 3; elif Type(H) eq MonStgElt then H := StringToInteger(H); end if;
if not assigned NParts then NParts := 1; elif Type(NParts) eq MonStgElt then NParts := StringToInteger(NParts); end if;
if not assigned Part then Part := 0; elif Type(Part) eq MonStgElt then Part := StringToInteger(Part); end if;
if not assigned MemGB then MemGB := 3; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
if not assigned MaxSteps then MaxSteps := 60; elif Type(MaxSteps) eq MonStgElt then MaxSteps := StringToInteger(MaxSteps); end if;
if not assigned MaxOrd then MaxOrd := 40; elif Type(MaxOrd) eq MonStgElt then MaxOrd := StringToInteger(MaxOrd); end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();
PQ<x> := PolynomialRing(Q);

// Polynomial part of sqrt(f), monic sextic f, char != 2.
function SqrtPolyPart(f)
    P := Parent(f);
    xx := P.1;
    s := xx^3;
    for k in [1..3] do
        d := f - s^2;
        if Degree(d) le 2 then break; end if;
        s := s + (Coefficient(d, 6-k)/(2*Coefficient(s, 3)))*xx^(3-k);
    end for;
    return s;
end function;

// Exact D_inf order via CF; 0 = no quasi-period within budget (order > maxord,
// non-torsion, or step budget exceeded -- NOT a proof of non-torsion).
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

////////////////////////////////////////////////////////////////////////
// SELF-TEST (mandatory: abort on any failure)
////////////////////////////////////////////////////////////////////////
f14 := (x^2+1)*(x^4+5*x^2+4*x+4);
f18 := (x^2-x+1)*(x^4-x^3+9*x^2+8*x-8);
f28t := x^6+2*x^5-5*x^4-14*x^3-3*x^2+24*x+28;
assert CFOrderB(f14, MaxSteps, MaxOrd) eq 14;
assert CFOrderB(f18, MaxSteps, MaxOrd) eq 18;
assert CFOrderB(f28t, MaxSteps, MaxOrd) eq 7;
for p in [101, 103] do
    Pp := PolynomialRing(GF(p));
    assert CFOrderB(Pp!f14, MaxSteps, MaxOrd) eq 14;
    assert CFOrderB(Pp!f18, MaxSteps, MaxOrd) eq 18;
    assert CFOrderB(Pp!f28t, MaxSteps, MaxOrd) eq 7;
end for;
printf "SELFTEST_PASS (f14->14, f18->18, f28trap->7; over Q and mod 101,103)\n";

////////////////////////////////////////////////////////////////////////
// Scan
////////////////////////////////////////////////////////////////////////
targets := {31, 35, 37, 38};
calibLo := 25; calibHi := MaxOrd;
primes := [101, 103, 107, 109];   // fallback list for bad reduction
Pps := [PolynomialRing(GF(p)) : p in primes];

range := [-H..H];
counts := AssociativeArray();
for N in [calibLo..calibHi] do counts[N] := 0; end for;
tested := 0; s2pass := 0; exactrun := 0; nearmiss := 0;
idx := -1;

for c5 in [0..3] do
 for c4 in range do
  for c3 in range do
   for c2 in range do
    for c1 in range do
     for c0 in range do
        idx +:= 1;
        if idx mod NParts ne Part then continue; end if;
        tested +:= 1;
        if tested mod 100000 eq 0 then
            printf "PROGRESS tested=%o s2pass=%o exact=%o\n", tested, s2pass, exactrun;
        end if;
        coeffs := [c0, c1, c2, c3, c4, c5, 1];
        // stage 1: first good prime
        o1 := -1;
        i1 := 0;
        for i in [1..#primes] do
            fb := Pps[i]!coeffs;
            if Degree(fb) eq 6 and Degree(GCD(fb, Derivative(fb))) eq 0 then
                o1 := CFOrderB(fb, MaxSteps, MaxOrd);
                i1 := i;
                break;
            end if;
        end for;
        if o1 lt calibLo then continue; end if;
        // stage 2: second good prime must agree
        o2 := -1;
        for i in [i1+1..#primes] do
            fb := Pps[i]!coeffs;
            if Degree(fb) eq 6 and Degree(GCD(fb, Derivative(fb))) eq 0 then
                o2 := CFOrderB(fb, MaxSteps, MaxOrd);
                break;
            end if;
        end for;
        if o2 ne o1 then continue; end if;
        s2pass +:= 1;
        // stage 3: exact over Q
        f := PQ!coeffs;
        if Degree(GCD(f, Derivative(f))) ne 0 then continue; end if;
        exactrun +:= 1;
        oQ := CFOrderB(f, MaxSteps, MaxOrd);
        if oQ eq 0 then nearmiss +:= 1; continue; end if;
        if oQ ge calibLo and oQ le calibHi then
            counts[oQ] +:= 1;
            printf "CALIB ord=%o f=%o\n", oQ, f;
        end if;
        if oQ in targets then
            printf "CFHIT ord=%o f=%o\n", oQ, f;
        end if;
     end for;
    end for;
   end for;
  end for;
 end for;
end for;

printf "CALIBRATION_TABLE (exact D_inf orders %o..%o in box H=%o part %o/%o):\n",
       calibLo, calibHi, H, Part, NParts;
for N in [calibLo..calibHi] do
    if counts[N] gt 0 then printf "  ORDER %o : %o\n", N, counts[N]; end if;
end for;
printf "SEARCH_DONE H=%o Part=%o/%o tested=%o s2pass=%o exact=%o nearmiss=%o\n",
       H, Part, NParts, tested, s2pass, exactrun, nearmiss;
quit;
