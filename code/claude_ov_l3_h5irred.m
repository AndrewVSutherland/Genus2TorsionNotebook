// ===========================================================================
// Lane 3 : prove that Z(H5) IS the Humbert surface H_5, not a larger locus.
//
// The derivation gives  Z(H5) ⊇ image(Elkies-Kumar parametrisation) = H_5,
// and Z(H5) is a hypersurface in the weighted P(2,4,6,10), so both are
// 2-dimensional.  They therefore agree as sets AS SOON AS H5 is ABSOLUTELY
// irreducible.  H5 is irreducible over Q (Magma Factorization, one factor of
// multiplicity 1).  For absolute irreducibility we use the standard descent:
//
//   if H5 mod p is absolutely irreducible for ONE prime p of good reduction,
//   then H5 is absolutely irreducible over Q;
//
// and over F_p, if H5 mod p is irreducible over F_p and factors into r
// absolute factors, Frobenius permutes them transitively, so each is defined
// over F_{p^r} with r <= deg H5 = 8.  Hence irreducibility over F_{p^k} for
// k = 1..8 IS absolute irreducibility.
//
// usage: code/claude_magma_slot.sh -b code/claude_ov_l3_h5irred.m
// ===========================================================================
SetColumns(0);
SetMemoryLimit(4*10^9);
Z := Integers();  Q := Rationals();
R<J2,J4,J6,J10> := PolynomialRing(Q, 4);
H5 := R ! eval Read("results/claude_ov_l3_humbert5_eqn.txt");
OUT := "results/claude_ov_l3_h5irred.log";
G := Open(OUT, "w");
procedure Say(G, s) printf "%o\n", s; Puts(G, s); end procedure;

Say(G, Sprintf("H5 #terms=%o totaldeg=%o", #Terms(H5), TotalDegree(H5)));
fac := Factorization(H5);
Say(G, Sprintf("factorisation over Q: %o factor(s), degrees %o, mults %o",
    #fac, [TotalDegree(t[1]) : t in fac], [t[2] : t in fac]));

// integral primitive model
den := LCM([Denominator(c) : c in Coefficients(H5)]);
Hi := R ! (den*H5);
cont := GCD([Z ! c : c in Coefficients(Hi)]);
Hi := R ! (Hi / cont);
Say(G, Sprintf("integral primitive: content=1, #terms=%o", #Terms(Hi)));

for p in [7, 11, 13, 101] do
    Fp := GF(p);
    Rp := PolynomialRing(Fp, 4);
    Hp := Rp ! Hi;
    if TotalDegree(Hp) ne TotalDegree(Hi) or #Terms(Hp) eq 0 then
        Say(G, Sprintf("p=%o : degree drop, skipped", p)); continue;
    end if;
    allirr := true;  firstbad := 0;
    for k in [1..8] do
        Fq := GF(p^k);
        Rq := PolynomialRing(Fq, 4);
        if not IsIrreducible(Rq ! Hi) then
            allirr := false; firstbad := k; break;
        end if;
    end for;
    Say(G, Sprintf("p=%o : irreducible over F_(p^k) for k=1..8 ? %o%o",
        p, allirr, allirr select "" else Sprintf("  (first reducible at k=%o)", firstbad)));
    if allirr then
        Say(G, Sprintf("ABSOLUTELY_IRREDUCIBLE via p=%o  =>  Z(H5) = H_5", p));
        break;
    end if;
end for;

Flush(G); delete G;
printf "IRRED_DONE\n";
quit;
