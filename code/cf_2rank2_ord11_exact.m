// [2,22] hunt: exact stage for C-sieve survivors (cf_2rank2_ord11_sieve.cpp).
// Input TSV lines: a b c d ord_p  (mod-p CF order agreed at >= 2 primes).
// For each: exact CF order of D_inf over Q; any order divisible by 11 gives
// J(Q) >= (Z/2)^2 x Z/11 = [2,22] on f=(x^2-1)(x^2+ax+b)(x^2+cx+d) --
// then exact TorsionSubgroup is computed and printed (HIT22 lines).
// Simplicity/End certification is a separate follow-up per hit.
//
// Run: magma -b InFile:=data/cf_2rank2_ord11_sieve_H100.tsv code/cf_2rank2_ord11_exact.m

SetColumns(0);
SetSeed(1);
if not assigned InFile then InFile := "data/cf_2rank2_ord11_sieve_H100.tsv"; end if;
if not assigned MemGB then MemGB := 3; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
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

assert CFOrderB((x^2+1)*(x^4+5*x^2+4*x+4), 60, 45) eq 14;
assert CFOrderB((x^2-x+1)*(x^4-x^3+9*x^2+8*x-8), 60, 45) eq 18;
printf "SELFTEST_PASS\n";

function ExactTorsionInvs(f)
    C := HyperellipticCurve(f);
    try
        Cm := ReducedMinimalWeierstrassModel(C);
        C := SimplifiedModel(Cm);
    catch e ; end try;
    return Invariants(TorsionSubgroup(Jacobian(C)));
end function;

lines := Split(Read(InFile), "\n");
printf "INPUT_LINES %o\n", #lines;
counts := AssociativeArray();
nhits := 0;
for l in lines do
    parts := Split(l, " ");
    if #parts lt 5 then continue; end if;
    a := StringToInteger(parts[1]); b := StringToInteger(parts[2]);
    c := StringToInteger(parts[3]); d := StringToInteger(parts[4]);
    op := StringToInteger(parts[5]);
    f := (x^2-1)*(x^2+a*x+b)*(x^2+c*x+d);
    if Discriminant(f) eq 0 then continue; end if;
    oQ := CFOrderB(f, 60, 45);
    printf "EXACT a=%o b=%o c=%o d=%o modp=%o overQ=%o\n", a, b, c, d, op, oQ;
    if oQ eq 0 then continue; end if;
    if IsDefined(counts, oQ) then counts[oQ] +:= 1; else counts[oQ] := 1; end if;
    if oQ mod 11 eq 0 then
        invs := ExactTorsionInvs(f);
        nhits +:= 1;
        printf "HIT22 a=%o b=%o c=%o d=%o Dinf_order=%o TORSION=%o f=%o\n",
               a, b, c, d, oQ, invs, f;
    end if;
end for;
printf "EXACT_TABLE:\n";
for o in Sort(Setseq(Keys(counts))) do
    printf "  ORDER %o : %o\n", o, counts[o];
end for;
printf "EXACT_DONE hits=%o\n", nhits;
quit;
