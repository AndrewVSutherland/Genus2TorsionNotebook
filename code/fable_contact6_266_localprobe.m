//////////////////////////////////////////////////////////////////////
// fable_contact6_266_localprobe.m  (2026-07-18, Fable session)
//
// Finite-field viability probe for [2,6,6] on the RAW contact-6 chart
// (a,b), using the factorized quintic from contact6_m36.md:
//
//   f = x*((b+3)x^2 + (a-3)x + 2)*(2x^2 + (b-3)x + (a+3)).
//
// The frontier note's [2,6,6] cover (discriminant-adjoined) is empty
// mod 5 and 7.  Here we ask the cruder chart-level question: does ANY
// good chart point over F_q have J(F_q) containing [2,6,6]?  If 0 at
// a good prime, the whole chart is locally obstructed for [2,6,6]
// (stronger conclusion than the cover being empty); if positive, the
// emptiness is specific to that cover's marked structure.
// Canaries: contains [6,6] and [2,6].
//
// Run:  magma -b code/fable_contact6_266_localprobe.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(2*10^9);

function ContainsInv(iv, ta)
    k := #iv; m := #ta;
    if k lt m then return false; end if;
    return &and[ IsDivisibleBy(iv[k-m+j], ta[j]) : j in [1..m] ];
end function;

for q in [5,7,11,13] do
    K := GF(q);
    Pol<x> := PolynomialRing(K);
    tested := 0; good := 0; c26 := 0; c66 := 0; c266 := 0;
    ex := [];
    for av0, bv0 in K do
        a := av0; b := bv0;
        if b eq -3 then continue; end if;           // degree drop
        f := x*((b+3)*x^2 + (a-3)*x + 2)*(2*x^2 + (b-3)*x + (a+3));
        tested +:= 1;
        if Degree(f) ne 5 or not IsSeparable(f) then continue; end if;
        okc := true;
        try
            C := HyperellipticCurve(f);
            G := AbelianGroup(Jacobian(C));
        catch ee okc := false; end try;
        if not okc then continue; end if;
        good +:= 1;
        iv := Invariants(G);
        if ContainsInv(iv, [2,6])   then c26  +:= 1; end if;
        if ContainsInv(iv, [6,6])   then c66  +:= 1; end if;
        if ContainsInv(iv, [2,6,6]) then
            c266 +:= 1;
            if #ex lt 6 then Append(~ex, <a,b,iv>); end if;
        end if;
    end for;
    printf "C6X266 q=%o tested=%o good=%o canary26=%o canary66=%o contains266=%o\n",
        q, tested, good, c26, c66, c266;
    for e in ex do printf "C6HIT q=%o a=%o b=%o inv=%o\n", q, e[1], e[2], e[3]; end for;
end for;

printf "PROBE_DONE\n";
quit;
