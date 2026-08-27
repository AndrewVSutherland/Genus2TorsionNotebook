// FINAL: rederive the descended rational model, fix the twist, exact torsion.
SetColumns(0);
SetMemoryLimit(12*10^9);
Q := Rationals();
K<w> := QuadraticField(11);
RK<xk> := PolynomialRing(K);
Px<x> := PolynomialRing(Q);
s0 := K!(4/5);
r0 := Roots(xk^2 - (s0^3-3*s0^2+4*s0)*xk + s0)[1][1];
c0 := s0*(r0-1); b0 := r0*c0;
E := EllipticCurve([1-c0, -b0, -b0, 0, 0]);
sig := hom< K -> K | -w >;
Es := EllipticCurve([sig(a) : a in aInvariants(E)]);
f := HyperellipticPolynomials(WeierstrassModel(E));
g := HyperellipticPolynomials(WeierstrassModel(Es));
L := SplittingField(f);
RL := PolynomialRing(L);
al := [ r[1] : r in Roots(RL!f) ];
gm := [ r[1] : r in Roots(RL!g) ];
DfK := Discriminant(f); DgK := Discriminant(g);
GA := Automorphisms(L);
GK := [ gA : gA in GA | gA(L!w) eq L!w ];
Gaction := func< rs | [[ Index(rs, gA(r)) : r in rs ] : gA in GK ] >;
aK := [[ al[i]-al[j] : j in [1..3]] : i in [1..3]];
RtL<tL> := PolynomialRing(L);
S := RK!0;
for sgm in SymmetricGroup(3) do
    beta := [ gm[i^sgm] : i in [1..3] ];
    b := [[ beta[i]-beta[j] : j in [1..3]] : i in [1..3]];
    if al[1]*b[3][2] + al[2]*b[1][3] + al[3]*b[2][1] ne 0 and Gaction(al) eq Gaction(beta) then
        a1 := aK[3][2]^2/b[3][2]+aK[2][1]^2/b[2][1]+aK[1][3]^2/b[1][3];
        a2 := al[1]*b[3][2]+al[2]*b[1][3]+al[3]*b[2][1];
        b1 := b[3][2]^2/aK[3][2]+b[2][1]^2/aK[2][1]+b[1][3]^2/aK[1][3];
        b2 := beta[1]*aK[3][2]+beta[2]*aK[1][3]+beta[3]*aK[2][1];
        A := DgK*a1/a2;  B := DfK*b1/b2;
        sext := -(A*aK[2][1]*aK[1][3]*tL^2+B*b[2][1]*b[1][3])
                *(A*aK[3][2]*aK[2][1]*tL^2+B*b[3][2]*b[2][1])
                *(A*aK[1][3]*aK[3][2]*tL^2+B*b[1][3]*b[3][2]);
        S := RK![ K!c : c in Coefficients(sext) ];
        break;
    end if;
end for;
// descent (as validated): M = [0,1;1,0], A = [1,w;0,1]
MKK := Matrix(K,2,2,[0,1,1,0]);
AM := Matrix(K,2,2,[K|1,w,0,1]);
AMs := Matrix(K,2,2,[sig(e) : e in Eltseq(AM)]);
N := AM + MKK*AMs;
nu := N[1,1]*xk + N[1,2]; de := N[2,1]*xk + N[2,2];
T6 := &+[ Coefficient(S,i)*nu^i*de^(6-i) : i in [0..6] ];
dg := Max([ i : i in [0..6] | Coefficient(T6,i) ne 0 ]);
T6s := RK![ sig(co) : co in Coefficients(T6) ];
lam2 := Coefficient(T6s,dg)/Coefficient(T6,dg);
assert T6s eq lam2*T6;
FQ := Px!0;
for cc in [K| 1, w, 1+w, 2+w, 1-2*w, 3+w, 1+2*w ] do
    mu := cc + lam2*sig(cc);
    if mu eq 0 then continue; end if;
    TT := mu*T6;
    if RK![ sig(co) : co in Coefficients(TT) ] eq TT then
        FQ := Px![ Q!co : co in Coefficients(TT) ];
        break;
    end if;
end for;
assert FQ ne 0;

// integralize + strip square content
denl := LCM([ Denominator(co) : co in Coefficients(FQ) ]);
fZ := Px!(FQ*denl^2);
cont := GCD([ Integers()!co : co in Coefficients(fZ) ]);
sq := 1;
for pr in Factorization(cont) do sq *:= pr[1]^(2*(pr[2] div 2)); end for;
fZ := Px![ co div sq : co in [Integers()!c : c in Coefficients(fZ)] ];
printf "integral model: y^2 = %o\n", fZ;
dsc := Integers()!Discriminant(fZ);
printf "disc trial division: %o\n", TrialDivision(AbsoluteValue(dsc), 10^6);

// twist sieve: 11 | #J_d(F_p) for all good p (split and inert)
goodp := [ p : p in PrimesInInterval(13, 250) | p ne 11 and dsc mod p ne 0 ];
OK := Integers(K);
winners := [];
for dsupp in Subsets({2,5,11}) do
    d0 := &*[ Integers() | p : p in dsupp ];
    for d in [d0, -d0] do
        ok := true; ntest := 0;
        for p in goodp do
            if ntest ge 14 then break; end if;
            if d mod p eq 0 then continue; end if;
            Fp := GF(p);
            fp := PolynomialRing(Fp)!(d*fZ);
            if Degree(fp) lt 5 or not IsSquarefree(fp) then continue; end if;
            ntest +:= 1;
            if #Jacobian(HyperellipticCurve(fp)) mod 11 ne 0 then ok := false; break; end if;
        end for;
        if ok and ntest ge 10 then
            printf "TWIST d=%o passes at %o primes\n", d, ntest;
            Append(~winners, d);
        end if;
    end for;
end for;
for d in winners do
    fT := d*fZ;
    dn2 := LCM([Denominator(c) : c in Coefficients(fT)]);
    fT := Px!(fT*dn2^2);
    C := HyperellipticCurve(fT);
    J := Jacobian(C);
    I := Invariants(TorsionSubgroup(J));
    printf "EXACT TORSION (twist %o): %o\n", d, I;
    // cross-check vs Res(E) at one split and one inert prime
    for p in [19, 29] do
        if dsc mod p eq 0 or d mod p eq 0 then continue; end if;
        Fp := GF(p);
        fp := PolynomialRing(Fp)!fT;
        nj := #Jacobian(HyperellipticCurve(fp));
        dec := Decomposition(OK, p);
        if #dec eq 2 then
            pred := #Reduction(E, dec[1][1]) * #Reduction(Es, dec[2][1]);
            printf "   p=%o split: #J=%o pred(E x Es)=%o\n", p, nj, pred;
        else
            pred := #Reduction(E, dec[1][1]);
            printf "   p=%o inert: #J=%o pred(#E(Fp2))=%o\n", p, nj, pred;
        end if;
    end for;
end for;
printf "S45FINAL_DONE winners=%o\n", winners;
quit;
