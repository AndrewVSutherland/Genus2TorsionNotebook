// Weil descent of the K-glue CK (s0=4/5, K=Q(sqrt 11)) to Q.
// Step 1: is the sextic already sigma-proportional?  If yes, scalar H90.
// Step 2: else find a Mobius M/K with S^sigma ~ S o M, solve N = A + M A^s
//         (matrix H90), substitute, rescale.
// Step 3: final rational model up to Q-twist; sieve twists by matching
//         #J(F_p) with the Res(E) prediction at split AND inert primes;
//         exact TorsionSubgroup of the winner.
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
S := RK!0; found := false;
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
        found := true;
        break;
    end if;
end for;
error if not found, "no matching";
Ss := RK![ sig(c) : c in Coefficients(S) ];
printf "deg S = %o (even sextic? coeffs of x^5,x^3,x: %o %o %o)\n",
    Degree(S), Coefficient(S,5), Coefficient(S,3), Coefficient(S,1);

// Step 1: sigma-proportional?
prop := false;
if Degree(S) eq Degree(Ss) then
    lam := Coefficient(Ss, Degree(Ss))/Coefficient(S, Degree(S));
    if Ss eq lam*S then
        prop := true;
        printf "S^sigma = lam * S with lam = %o, Norm(lam) = %o\n", lam, Norm(lam);
        // want mu with mu^sigma * lam = mu  =>  mu = a + lam-solving; try
        // mu = 1 + lam-corrected (Hilbert 90): mu := c + lam*c^sigma for generic c
        for cc in [K| 1, w, 1+w, 2+w, 1-2*w ] do
            mu := cc + lam*sig(cc);
            if mu eq 0 then continue; end if;
            T := mu*S;
            Ts := RK![ sig(co) : co in Coefficients(T) ];
            if Ts eq T then
                FQ := Px![ Q!co : co in Coefficients(T) ];
                printf "RATIONAL MODEL: y^2 = %o\n", FQ;
                break;
            end if;
        end for;
    end if;
end if;
if not prop then
    printf "not sigma-proportional; even-sextic Mobius descent\n";
    a6 := Coefficient(S,6); a4 := Coefficient(S,4); a2 := Coefficient(S,2); a0 := Coefficient(S,0);
    printf "coeffs a6,a4,a2,a0 nonzero: %o %o %o %o\n", a6 ne 0, a4 ne 0, a2 ne 0, a0 ne 0;
    // shape 1: S^sigma(x) = mu * x^6 S(kap/x):  a6^s = mu a0, a4^s = mu a2 k2,
    //          a2^s = mu a4 k2^2, a0^s = mu a6 k2^3   (k2 = kap^2)
    mu1 := sig(a6)/a0;
    k2 := sig(a4)/(mu1*a2);
    ok1 := sig(a2) eq mu1*a4*k2^2 and sig(a0) eq mu1*a6*k2^3;
    printf "shape x->kap/x: consistent %o, k2 = %o\n", ok1, k2;
    // shape 2: S^sigma(x) = mu * S(lam x)
    mu2 := sig(a6)/a6;
    l2 := sig(a4)/(mu2*a4);
    ok2 := a2 ne 0 and sig(a2) eq mu2*a2*l2^2 and sig(a0) eq mu2*a0*l2^3;
    printf "shape x->lam x: consistent %o\n", ok2;
    MKK := 0; haveM := false;
    if ok1 then
        okk, kap := IsSquare(k2);
        printf "k2 square in K: %o\n", okk;
        if okk then MKK := Matrix(K,2,2,[0,kap,1,0]); haveM := true;
        else
            // use M with M(x) = k2/ (x)?  matrix [0,k2,1,0]: M(x)=k2/x: then
            // x^6 S(k2/x) differs from S(kap/x)-scaling by even powers: still fine:
            MKK := Matrix(K,2,2,[0,k2,1,0]); haveM := true;
        end if;
    elif ok2 then
        okk, lam1 := IsSquare(l2);
        if okk then MKK := Matrix(K,2,2,[lam1,0,0,1]); haveM := true; end if;
    end if;
    if haveM then
        Mks := Matrix(K, 2,2, [ sig(e) : e in Eltseq(MKK) ]);
        CC := Mks*MKK;
        issc := CC[1,2] eq 0 and CC[2,1] eq 0 and CC[1,1] eq CC[2,2];
        printf "M = %o\ncocycle scalar: %o (%o)\n", MKK, issc, CC[1,1];
        done := false;
        if issc then
            for Aseq in [[K|1,0,0,1],[K|1,1,0,1],[K|1,w,0,1],[K|1,0,w,1],[K|2,w,1,1],[K|1,2,w,1]] do
                AM := Matrix(K,2,2,Aseq);
                AMs := Matrix(K,2,2,[sig(e) : e in Eltseq(AM)]);
                N := AM + MKK*AMs;
                if Determinant(N) eq 0 then continue; end if;
                nu := N[1,1]*xk + N[1,2]; de := N[2,1]*xk + N[2,2];
                T6 := &+[ Coefficient(S,i)*nu^i*de^(6-i) : i in [0..6] ];
                dg := Max([ i : i in [0..6] | Coefficient(T6,i) ne 0 ]);
                T6s := RK![ sig(co) : co in Coefficients(T6) ];
                lam2 := sig(Coefficient(T6,dg))/Coefficient(T6,dg);
                lam2 := Coefficient(T6s,dg)/Coefficient(T6,dg);
                if T6s eq lam2*T6 then
                    printf "substitution A=%o: sigma-proportional, Norm(lam)=%o\n", Aseq, Norm(lam2);
                    for cc in [K| 1, w, 1+w, 2+w, 1-2*w, 3+w, 1+2*w ] do
                        mu := cc + lam2*sig(cc);
                        if mu eq 0 then continue; end if;
                        TT := mu*T6;
                        TTs := RK![ sig(co) : co in Coefficients(TT) ];
                        if TTs eq TT then
                            FQ := Px![ Q!co : co in Coefficients(TT) ];
                            printf "RATIONAL MODEL: y^2 = %o\n", FQ;
                            done := true;
                            break;
                        end if;
                    end for;
                else
                    printf "substitution A=%o: NOT proportional\n", Aseq;
                end if;
                if done then break; end if;
            end for;
        end if;
    end if;
end if;
printf "S45DESCEND_DONE\n";
quit;
