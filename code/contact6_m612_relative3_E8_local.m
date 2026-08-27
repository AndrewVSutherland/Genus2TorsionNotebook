//////////////////////////////////////////////////////////////////////
// LOCAL SOLUBILITY of E8 (genus 4) = { mp4(Y,e)=0, w^2 = Y }.
// E8(Q_p) nonempty iff some (Y0,e0) in C(Q_p) has Y0 a square (or 0).
// Crude-but-informative scan: for disk representatives e0 = u*p^k,
// solve the quartic mp4(Y, e0) over Q_p (Roots), test IsSquare(Y0).
// Also the e=infinity chart via the limiting quartic.
// A prime with NO square-Y point across all disks+depths is a
// candidate local obstruction (then to be made rigorous separately).
//
// Usage: magma -b code/contact6_m612_relative3_E8_local.m
//////////////////////////////////////////////////////////////////////
SetColumns(0); SetSeed(1); SetMemoryLimit(4*10^9);
Q := Rationals();
K<e> := FunctionField(Q); Kz<zz> := PolynomialRing(K);
mp4Q := Kz ! eval Read("data/contact6_m612_E4_mp4Q.txt");

primes := [2,3,5,7,11,13,17,19,23,29,31,37];
prec := 40;

for p in primes do
    Qp := pAdicField(p, prec);
    PQp<T> := PolynomialRing(Qp);
    found := false; witness := <0,0>;
    // unit-times-p^k disk representatives
    m := p le 7 select 3 else (p le 19 select 2 else 1);
    units := [u : u in [1..p^m-1] | u mod p ne 0];
    for k in [-5..5] do
        if found then break; end if;
        for u in units do
            e0 := Qp!u * (Qp!p)^k;
            // quartic coefficients c_i(e0)
            cs := [];
            okc := true;
            for i in [0..4] do
                c := Coefficient(mp4Q, i);
                nn := Numerator(c); dd := Denominator(c);
                nv := &+[Qp | Qp!(Coefficient(nn,h)) * e0^h : h in [0..Degree(nn)]];
                dv := &+[Qp | Qp!(Coefficient(dd,h)) * e0^h : h in [0..Degree(dd)]];
                if Valuation(dv) ge prec-5 then okc := false; break; end if;
                Append(~cs, nv/dv);
            end for;
            if not okc then continue; end if;
            quart := &+[cs[i+1]*T^i : i in [0..4]];
            rts := [];
            try rts := Roots(quart); catch ee rts := []; end try;
            for r in rts do
                Y0 := r[1];
                if Valuation(Y0) ge prec-8 then found := true; witness := <u,k>; break; end if; // Y ~ 0: square
                if IsEven(Valuation(Y0)) then
                    uY := Y0 / (Qp!p)^Valuation(Y0);
                    if IsSquare(uY) then found := true; witness := <u,k>; break; end if;
                end if;
            end for;
            if found then break; end if;
        end for;
    end for;
    // e = infinity chart: limiting quartic Y^4 + 216 Y^2 - 1296 Y - 3888
    if not found then
        quartinf := T^4 + 216*T^2 - 1296*T - 3888;
        rts := [];
        try rts := Roots(quartinf); catch ee rts := []; end try;
        for r in rts do
            Y0 := r[1];
            if IsEven(Valuation(Y0)) and IsSquare(Y0/(Qp!p)^Valuation(Y0)) then
                found := true; witness := <-999,-999>;
            end if;
        end for;
    end if;
    printf "p=%o : E8 local point found = %o %o\n", p, found,
        found select Sprintf("(disk u=%o k=%o)", witness[1], witness[2]) else "  <-- CANDIDATE OBSTRUCTION";
end for;
print "E8_LOCAL_SCAN_DONE";
quit;
