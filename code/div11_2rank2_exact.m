// [2,22] hunt: exact stage for 11-divisibility sieve survivors
// (code/div11_2rank2_sieve.cpp).  Input lines: a b c d.
// For each: exact TorsionSubgroup of Jac(y^2 = (x^2-1)(x^2+ax+b)(x^2+cx+d)).
// Any torsion with 11 | order has, together with the built-in 2-rank 2,
// J(Q) >= [2,22].  HIT22 lines print the full data.
//
// Run: magma -b InFile:=data/div11_2rank2_sieve_H150.tsv code/div11_2rank2_exact.m

SetColumns(0);
SetSeed(1);
if not assigned InFile then InFile := "data/div11_2rank2_sieve_H150.tsv"; end if;
if not assigned MemGB then MemGB := 3; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();
P<x> := PolynomialRing(Q);

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
nhits := 0; nrun := 0;
for l in lines do
    parts := Split(l, " ");
    if #parts lt 4 then continue; end if;
    a := StringToInteger(parts[1]); b := StringToInteger(parts[2]);
    c := StringToInteger(parts[3]); d := StringToInteger(parts[4]);
    f := (x^2-1)*(x^2+a*x+b)*(x^2+c*x+d);
    if Discriminant(f) eq 0 then continue; end if;
    // tier-2 prefilter: 11 | #J(F_p) at further primes kills chance
    // survivors of the C chain at ~1/11 per prime
    disc := Integers()!Discriminant(f);
    ok := true; ngood2 := 0;
    for pp in [31, 37, 41, 43, 47, 53] do
        if disc mod pp eq 0 then continue; end if;
        Jp := Jacobian(HyperellipticCurve(ChangeRing(f, GF(pp))));
        if #Jp mod 11 ne 0 then ok := false; break; end if;
        ngood2 +:= 1;
        if ngood2 ge 3 then break; end if;
    end for;
    if not ok then continue; end if;
    nrun +:= 1;
    invs := ExactTorsionInvs(f);
    printf "EXACT a=%o b=%o c=%o d=%o torsion=%o\n", a, b, c, d, invs;
    if #invs ge 1 and (&*invs) mod 11 eq 0 then
        nhits +:= 1;
        printf "HIT22 a=%o b=%o c=%o d=%o TORSION=%o f=%o\n", a, b, c, d, invs, f;
    end if;
end for;
printf "EXACT_DONE run=%o hits=%o\n", nrun, nhits;
quit;
