
//////////////////////////////////////////////////////////////////////
//  Second-stage obstruction diagnostics for P_R halving in M_1(8,4).
//
//  On the C1&C2 locus (Pell-parametrized, C2 = square), the remaining
//  halving conditions are
//    S_A: u_A = X_R - theta_A  square in K_A = Q[T]/At,
//    S_B: u_B = X_R - theta_B  square in K_B = Q[T]/Bt.
//  For a quadratic field K = Q(sqrt(d)) and u = alpha + beta*sqrt(d)
//  with N(u) = n^2 (n rational, sign choice +-):
//    u in K*^2  <=>  2*(alpha + n) = square  or  2*(alpha - n) = square
//  (equivalently: the quaternion algebra (2*(alpha+n), d) splits).
//
//  Empirically (h50 run) S_A/S_B kill 36/36 C1&C2 points.  This script
//  re-enumerates the C1&C2 points and prints, for each:
//    d_A, d_B (discriminant squareclasses),
//    the four candidate square arguments 2*(alpha +- n) for A and B,
//    and the ramified primes of the corresponding Hilbert classes,
//  hunting a CONSTANT nontrivial Brauer class (a provable obstruction).
//
//  Usage: magma -b hR:=50 hM:=50 agent_m18_416_sstage_obstruction.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();

if not assigned hR then hR := 50;
elif Type(hR) eq MonStgElt then hR := StringToInteger(hR); end if;
if not assigned hM then hM := 50;
elif Type(hM) eq MonStgElt then hM := StringToInteger(hM); end if;

P<xq> := PolynomialRing(Q);

function IsSquareQ(qv)
    qv := Q!qv;
    if qv le 0 then return false; end if;
    return IsSquare(Numerator(qv)) and IsSquare(Denominator(qv));
end function;

function SqfreePart(qv)   // squarefree integer representative of the class
    qv := Q!qv;
    if qv eq 0 then return 0; end if;
    n := Z!(Numerator(qv)*Denominator(qv));
    s := Sign(n); n := Abs(n);
    sf := 1;
    for pe in Factorization(n) do
        if IsOdd(pe[2]) then sf *:= pe[1]; end if;
    end for;
    return s*sf;
end function;

function RamPlacesQA(aa, bb)  // of quaternion algebra (aa,bb) over Q
    a0 := SqfreePart(aa); b0 := SqfreePart(bb);
    if a0 eq 1 or b0 eq 1 then return []; end if;
    A := QuaternionAlgebra<Q | a0, b0>;
    return RamifiedPlaces(A);  // intrinsic
end function;

function FamilyData(Rv, wv)
    tv := (2*Rv^2 + (1-wv^2)*Rv - 2*wv^2)/(4*(wv^2-1));
    c4v := Rv + 2 + 4*tv;
    Av := xq^2 + (Rv^3 + 4*Rv^2*tv + Rv - 8*Rv*tv + 4*tv)*xq + Rv^4;
    Bv := c4v*xq^2 + (Rv^2 + 4*Rv + 1 + 8*tv)*xq + (2*Rv^2 + Rv + 4*tv);
    return xq*Av*Bv, Av, Bv, c4v;
end function;

function RationalParametersOfHeight(Bd)
    vals := [];
    for den in [1..Bd] do for num in [-Bd..Bd] do
        if GCD(num, den) eq 1 then Append(~vals, Q!num/den); end if;
    end for; end for;
    return Sort(Setseq(Seqset(vals)));
end function;

printf "S-STAGE OBSTRUCTION DIAGNOSTICS hR=%o hM=%o\n", hR, hM;
Rparams := RationalParametersOfHeight(hR);
Mparams := RationalParametersOfHeight(hM);
found := 0;
seenW := {};

for Rv in Rparams do
    if Rv in {Q!0, Q!1, Q!-1} then continue; end if;
    K := -2*Rv*(Rv^2 - 1);
    for mv in Mparams do
        den := mv^2 - K;
        if den eq 0 then continue; end if;
        wv := (mv^2 + K)/den;
        if wv in {Q!0, Q!1, Q!-1} then continue; end if;
        if <Rv,wv> in seenW then continue; end if;
        Include(~seenW, <Rv,wv>);
        Wv := wv^2;
        G := 2*(Rv^2-1)*(Rv*(2*Rv+1) - Wv*(Rv+2));
        if not IsSquareQ(G) then continue; end if;
        found +:= 1;
        // second-stage data
        fq, Av, Bv, c4v := FamilyData(Rv, wv);
        if Degree(fq) ne 5 or Discriminant(fq) eq 0 then continue; end if;
        XR := -c4v*Rv;
        printf "\nPOINT %o: R=%o w=%o (m=%o)\n", found, Rv, wv, mv;
        for pair in [<"A", c4v^2*Evaluate(Av, xq/c4v)>, <"B", c4v*Evaluate(Bv, xq/c4v)>] do
            nm := pair[1];
            gq := P![Q!co : co in Coefficients(pair[2])];
            dsc := Discriminant(gq);
            dcls := SqfreePart(dsc);
            // u = XR - theta: alpha = XR - (-g1/2) = XR + g1/2 where
            // gq = T^2 + g1*T + g0;  N(u) = gq(XR).
            g1 := Coefficient(gq,1);
            alpha := XR + g1/2;
            Nu := Evaluate(gq, XR);
            okN, n := IsSquare(Q!Nu);
            printf "  %o: d=%o  N(u)=%o square=%o", nm, dcls, Nu, okN;
            if dcls eq 1 then
                // split component: report the two values XR - e_i
                rts := Roots(gq);
                printf "  SPLIT vals=%o squares=%o\n",
                    [XR - r[1] : r in rts],
                    [IsSquareQ(XR - r[1]) : r in rts];
                continue;
            end if;
            if not okN then printf "  (norm not square: C-level fail)\n"; continue; end if;
            v1 := 2*(alpha + n); v2 := 2*(alpha - n);
            s1 := IsSquareQ(v1); s2 := IsSquareQ(v2);
            printf "  2(a+n)=%o sq=%o | 2(a-n)=%o sq=%o\n", v1, s1, v2, s2;
            // Hilbert class of the failing torsor: (v1, d) ~ (v2, d)
            if not (s1 or s2) and v1 ne 0 then
                ram := RamPlacesQA(v1, dsc);
                printf "    quaternion (2(a+n), d) ramified at: %o\n",
                    [Sprint(pl) : pl in ram];
            end if;
        end for;
    end for;
end for;
printf "\ntotal C1&C2 points examined: %o\n", found;
print "DONE";
quit;
