// sigma-3-congruence test on the five LMFDB quadratic 11-torsion pairs:
// necessary: a_p(E) = a_p(E^sigma) mod 3 at all good primes of K;
// module test: 3-division x-quartics generate matching etale algebras.
SetColumns(0);
SetMemoryLimit(8*10^9);
Q := Rationals();
QT<T> := PolynomialRing(Q);
DATA := [*
  <"2.2.8.1-46.2-b2",   T^2-2,   [[1,1],[-1,-1],[1,0],[-3,1],[3,-2]]>,
  <"2.2.17.1-172.2-f1", T^2-T-4, [[1,0],[-1,-1],[1,1],[-27,11],[455,-179]]>,
  <"2.0.7.1-268.3-b1",  T^2-T+2, [[1,0],[-1,1],[1,1],[-2,-2],[7,-5]]>,
  <"2.2.13.1-828.2-d1", T^2-T-3, [[1,0],[0,0],[0,0],[-104,-297],[-29541,16497]]>,
  <"2.0.8.1-4338.4-b2", T^2+2,   [[1,1],[-1,1],[1,0],[122,62],[-483,286]]>
*];
for rec in DATA do
    lbl := rec[1];
    K<wg> := NumberField(rec[2]);
    OK := Integers(K);
    conj := [ rr[1] : rr in Roots(PolynomialRing(K)!rec[2]) | rr[1] ne wg ][1];
    sig := hom< K -> K | conj >;
    ai := [ K!(c[1] + c[2]*wg) : c in rec[3] ];
    E := EllipticCurve(ai);
    Es := EllipticCurve([ sig(a) : a in ai ]);
    // trace test mod 3
    mm := 0; nt := 0;
    for p in PrimesInInterval(5, 200) do
        for dd in Decomposition(OK, p) do
            pr := dd[1];
            okr := true; t1 := 0; t2 := 0;
            try
                t1 := TraceOfFrobenius(Reduction(E, pr));
                t2 := TraceOfFrobenius(Reduction(Es, pr));
            catch e; okr := false; end try;
            if not okr then continue; end if;
            nt +:= 1;
            if (t1 - t2) mod 3 ne 0 then mm +:= 1; end if;
        end for;
    end for;
    printf "%o: trace-mod-3 mismatches %o/%o", lbl, mm, nt;
    if mm eq 0 then
        // module test: quartic field match
        RK := PolynomialRing(K);
        p3 := RK!DivisionPolynomial(E, 3);
        p3s := RK![ sig(Coefficient(RK!DivisionPolynomial(Es,3), i)) : i in [0..4] ];
        // careful: p3s should be div poly of Es viewed via sigma-matching to E's:
        p3s := RK!DivisionPolynomial(Es, 3);
        f3 := Factorization(p3);
        f3s := Factorization(p3s);
        d3 := Sort([Degree(ff[1]) : ff in f3]);
        d3s := Sort([Degree(ff[1]) : ff in f3s]);
        okm := d3 eq d3s;
        if okm and Degree(f3[#f3][1]) gt 1 then
            L3 := ext< K | f3[#f3][1] >;
            okm := exists{ ff : ff in f3s | Degree(ff[1]) eq Degree(f3[#f3][1]) and #Roots(PolynomialRing(L3)!ff[1]) gt 0 };
        end if;
        printf "  MODULE3 match: %o", okm;
    end if;
    printf "\n";
end for;
printf "LMFDB11_C3_DONE\n";
quit;
