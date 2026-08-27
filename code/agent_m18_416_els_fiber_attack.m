
//////////////////////////////////////////////////////////////////////
//  Attack on the everywhere-locally-solvable (ELS) fibers of the
//  P_R-halving problem in M_1(8,4).
//
//  agent_m18_416_fiber_local_scan.m: fibers R = -8, -25/4, -29/8 are
//  ELS up to p <= 23 for the full halving fiber product.  This script
//  attacks each fiber:
//   1. builds the genus-1 model  y^2 = Gamma(m)  of the C2 curve
//      (Gamma = numerator of G(R,m)*(m^2-K)^2, an even quartic in m);
//   2. deep m-sweep: all rational m of height <= sweepH; every C2 point
//      (G = square) gets the exact second-stage descent test;
//   3. Mordell-Weil: converts to an elliptic curve via a known point,
//      computes rank bounds / generators, enumerates MW combinations,
//      maps back to m, and runs the exact second-stage test on each.
//  Any full pass ==> rational half of P_R ==> order-16 point (>= [2,16],
//  and [4,16] if the first tangent cover also has a point).
//
//  Usage: magma -b sweepH:=800 mwN:=12 agent_m18_416_els_fiber_attack.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();

if not assigned sweepH then sweepH := 500;
elif Type(sweepH) eq MonStgElt then sweepH := StringToInteger(sweepH); end if;
if not assigned mwN then mwN := 10;
elif Type(mwN) eq MonStgElt then mwN := StringToInteger(mwN); end if;

fibers := [<Q!-8, Q!-28>, <Q!-25/4, Q!-25/2>, <Q!-29/8, Q!-29/4>];

P<xq> := PolynomialRing(Q);

function IsSquareQ(qv)
    qv := Q!qv;
    if qv le 0 then return false, 0; end if;
    okn, sn := IsSquare(Numerator(qv));
    okd, sd := IsSquare(Denominator(qv));
    if okn and okd then return true, sn/sd; end if;
    return false, 0;
end function;

function FamilyData(Rv, wv)
    tv := (2*Rv^2 + (1-wv^2)*Rv - 2*wv^2)/(4*(wv^2-1));
    c4v := Rv + 2 + 4*tv;
    Av := xq^2 + (Rv^3 + 4*Rv^2*tv + Rv - 8*Rv*tv + 4*tv)*xq + Rv^4;
    Bv := c4v*xq^2 + (Rv^2 + 4*Rv + 1 + 8*tv)*xq + (2*Rv^2 + Rv + 4*tv);
    return xq*Av*Bv, Av, Bv, c4v;
end function;

function SecondStagePass(Rv, wv)
    fq, Av, Bv, c4v := FamilyData(Rv, wv);
    if Degree(fq) ne 5 or Discriminant(fq) eq 0 then return false; end if;
    XR := -c4v*Rv;
    Atq := P![Q!co : co in Coefficients(c4v^2*Evaluate(Av, xq/c4v))];
    Btq := P![Q!co : co in Coefficients(c4v*Evaluate(Bv, xq/c4v))];
    for gq in [Atq, Btq] do
        dsc := Discriminant(gq);
        if dsc eq 0 then return false; end if;
        okd, _ := IsSquareQ(dsc);
        if okd then
            for rt in Roots(gq) do
                val := XR - rt[1];
                if val eq 0 then return false; end if;
                okv, _ := IsSquareQ(val);
                if not okv then return false; end if;
            end for;
        else
            Kf<th> := NumberField(gq);
            if not IsSquare(XR - th) then return false; end if;
        end if;
    end for;
    return true;
end function;

procedure CertifyAndReport(Rv, mv, wv)
    printf "!!!! FULL DESCENT PASS R=%o m=%o w=%o -- certifying\n", Rv, mv, wv;
    fq, Av, Bv, c4v := FamilyData(Rv, wv);
    L := 1;
    for i in [0..Degree(fq)] do L := LCM(L, Denominator(Coefficient(fq, i))); end for;
    fI := P!(L^2*fq);
    J := Jacobian(HyperellipticCurve(fI));
    Qf := Rv^2 - (Q!1/2)*Rv*wv^2 + (Q!1/2)*Rv - wv^2;
    YRv := -2*Rv*(Rv-1)^2*Qf/(wv^2-1);
    PR := J![xq + Rv, P!(L*YRv)];
    ok, half := IsDivisibleBy(PR, 2);
    printf "IsDivisibleBy(P_R,2) = %o\n", ok;
    if ok then
        printf "half order = %o\n", Order(half);
        invs := Invariants(TorsionSubgroup(J));
        printf "TORSION = %o\n", invs;
        printf "f = %o\n", fI;
    end if;
end procedure;

for fib in fibers do
    Rv := fib[1]; mex := fib[2];
    Kv := -2*Rv*(Rv^2-1);
    printf "\n================ FIBER R = %o ================\n", Rv;
    // Gamma(m): numerator of G*(m^2-K)^2
    Pm<mm> := PolynomialRing(Q);
    Gam := 2*(Rv^2-1)*(Rv*(2*Rv+1)*(mm^2-Kv)^2 - (Rv+2)*(mm^2+Kv)^2);
    printf "Gamma(m) = %o\n", Gam;
    printf "disc = %o\n", Discriminant(Gam);

    // ---- deep m-sweep ----
    c2pts := [];
    tested := 0;
    for den in [1..sweepH] do
        for num in [-sweepH..sweepH] do
            if GCD(num, den) ne 1 then continue; end if;
            mv := Q!num/den;
            if mv eq 0 then continue; end if;
            tested +:= 1;
            Gv := Evaluate(Gam, mv);
            okG, _ := IsSquareQ(Gv);
            if not okG then continue; end if;
            den0 := mv^2 - Kv;
            if den0 eq 0 then continue; end if;
            wv := (mv^2 + Kv)/den0;
            if wv in {Q!0, Q!1, Q!-1} then continue; end if;
            Append(~c2pts, <mv, wv>);
            if SecondStagePass(Rv, wv) then
                CertifyAndReport(Rv, mv, wv);
            end if;
        end for;
    end for;
    printf "sweep height %o: tested=%o C2 points=%o (second stage reported above if any)\n",
        sweepH, tested, #c2pts;
    printf "C2 m-values: %o\n", [pt[1] : pt in c2pts];

    // ---- Mordell-Weil attack ----
    try
        C := HyperellipticCurve(Gam);
        // known rational point from mex
        Gv0 := Evaluate(Gam, mex);
        okg, g0 := IsSquareQ(Gv0);
        error if not okg, "example point not on curve";
        pt0 := C![mex, g0];
        E, phi := EllipticCurve(C, pt0);
        Emin, mpmin := MinimalModel(E);
        printf "elliptic fiber: %o\n", aInvariants(Emin);
        r1, r2 := RankBounds(Emin);
        printf "rank bounds: [%o, %o]\n", r1, r2;
        T, mT := TorsionSubgroup(Emin);
        printf "torsion: %o\n", Invariants(T);
        if r2 eq 0 and #Invariants(T) le 1 and #T le 4 then
            printf "rank 0: fiber has only finitely many C2 points (torsion); enumerating those via sweep is sufficient.\n";
        end if;
        if r1 ge 1 then
            gens := Generators(Emin);
            printf "generators found: %o\n", #gens;
            // enumerate MW combos and map back to m
            phiInv := Inverse(phi);
            mpminInv := Inverse(mpmin);
            torspts := [mT(tt) : tt in T];
            count := 0;
            for n1 in [-mwN..mwN] do
                base := n1*gens[1];
                for tt in torspts do
                    Ept := base + tt;
                    if Ept eq Emin!0 then continue; end if;
                    try
                        Cpt := phiInv(mpminInv(Ept));
                        mv := Cpt[1]/Cpt[3];
                        den0 := mv^2 - Kv;
                        if den0 eq 0 then continue; end if;
                        wv := (mv^2 + Kv)/den0;
                        if wv in {Q!0, Q!1, Q!-1} then continue; end if;
                        count +:= 1;
                        if SecondStagePass(Rv, wv) then
                            CertifyAndReport(Rv, mv, wv);
                        end if;
                    catch e2
                        continue;
                    end try;
                end for;
            end for;
            printf "MW enumeration: %o points tested through n=%o\n", count, mwN;
        end if;
    catch err
        printf "MW attack failed: %o\n", err`Object;
    end try;
end for;
print "DONE";
quit;
