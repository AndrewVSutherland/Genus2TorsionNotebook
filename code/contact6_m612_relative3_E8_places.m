SetColumns(0);
Q := Rationals();
K<e> := FunctionField(Q); Kz<zz> := PolynomialRing(K);
mp4Q := Kz ! eval Read("data/contact6_m612_E4_mp4Q.txt");
// mp8(w) = mp4(w^2)
mp8 := Kz ! [ IsEven(i) select Coefficient(mp4Q, i div 2) else 0 : i in [0..8] ];
assert Degree(mp8) eq 8;
F8<w8> := ext<K | mp8>;
printf "E8 genus (recheck) = %o\n", Genus(F8);
for pl in [<"e=0", Zeros(K!e)[1]>, <"e=inf", Poles(K!e)[1]>] do
    printf "== places over %o ==\n", pl[1];
    dec := Decomposition(F8, pl[2]);
    for D in dec do
        P := Type(D) eq Tup select D[1] else D;
        rf := ResidueClassField(P);
        printf "  degree %o, ram %o, residue field: %o\n",
            Degree(P), RamificationIndex(P), rf;
        if Degree(P) eq 1 then
            printf "    ** RATIONAL boundary place (residue field Q) **\n";
        elif Type(rf) eq FldNum then
            // does it embed in Q_2? check 2-adic factorization of its defining poly
            dp := DefiningPolynomial(rf);
            f2 := Factorization(ChangeRing(dp, pAdicField(2, 60)));
            degs := [Degree(g[1]) : g in f2];
            printf "    2-adic completion degrees: %o %o\n", degs,
                1 in degs select "  ** embeds in Q_2: local points exist near boundary **" else "";
        end if;
    end for;
end for;
print "DONE";
quit;
