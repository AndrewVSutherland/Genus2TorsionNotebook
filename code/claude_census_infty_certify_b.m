// Stragglers: [2,2] via x(x-1)(x-2)(x^3+x+t) (even, 4 factors, rank 2),
// [2,2,2] via x(x-1)(x-2)(x^2+x+t) (odd deg-5, 4 factors, rank 3).
// Strict prime range extended to 600.
SetColumns(0);
SetSeed(1);
SetMemoryLimit(8*10^9);
Q := Rationals();
P<x> := PolynomialRing(Q);
function StrictPrime(C)
    fC, hC := HyperellipticPolynomials(C);
    g := 4*fC + hC^2;
    dsc := Integers()!Numerator(Discriminant(g)/16);
    for p in PrimesInInterval(11, 600) do
        if dsc mod p eq 0 then continue; end if;
        Cp := ChangeRing(SimplifiedModel(C), GF(p));
        Pp := P!Reverse(Coefficients(LPolynomial(Cp)));
        if not IsIrreducible(Pp) then continue; end if;
        K<pi> := NumberField(Pp);
        ok := true;
        for n in [1..12] do
            if Degree(MinimalPolynomial(pi^n)) ne 4 then ok := false; break; end if;
        end for;
        if ok then return p, true; end if;
    end for;
    return 0, false;
end function;
fams := [* <"[2,2]", func< t | x*(x-1)*(x-2)*(x^3 + x + t) >>,
          <"[2,2,2]", func< t | x*(x-1)*(x-2)*(x^2 + x + t) >> *];
tvals := [ Q!3, Q!5, Q!7, Q!-3, Q!11, Q!1/2, Q!2/5, Q!13 ];
for fm in fams do
    G := fm[1];
    exactT := 0; strictT := 0; strictP := 0; invs := {};
    printf "FAMILY %o\n", G;
    for t in tvals do
        f := fm[2](t);
        if not IsSquarefree(f) then continue; end if;
        den := LCM([Denominator(c) : c in Coefficients(f)]);
        f := den^2*f;
        C := HyperellipticCurve(f);
        Include(~invs, G2Invariants(C));
        T := TorsionSubgroup(Jacobian(SimplifiedModel(C)));
        inv := Invariants(T);
        tstr := "[" cat (#inv eq 0 select "" else
            &cat[ IntegerToString(inv[k]) cat (k lt #inv select "," else "") : k in [1..#inv] ]) cat "]";
        printf "  t=%o torsion %o\n", t, tstr;
        if tstr eq G then
            if Type(exactT) eq RngIntElt then exactT := t; end if;
            if strictP eq 0 then
                sp, oks := StrictPrime(C);
                if oks then strictT := t; strictP := sp; end if;
            end if;
        end if;
        if Type(exactT) ne RngIntElt and strictP ne 0 and #invs ge 2 then break; end if;
    end for;
    printf "CERT %o exact_fiber=%o strict_fiber=%o strict_prime=%o distinct_invs=%o\n",
        G, exactT, strictT, strictP, #invs;
end for;
printf "DONE_B\n";
quit;
