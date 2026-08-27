//////////////////////////////////////////////////////////////////////
// Exact A/B second-stage diagnostics for selected (R,w) points.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
P<x> := PolynomialRing(Q);

function IsSquareQ(qv)
    qv := Q!qv;
    if qv lt 0 then return false, 0; end if;
    okn, sn := IsSquare(Numerator(qv));
    okd, sd := IsSquare(Denominator(qv));
    if okn and okd then return true, sn/sd; end if;
    return false, 0;
end function;

function FamilyData(Rv, wv)
    tv := (2*Rv^2 + (1-wv^2)*Rv - 2*wv^2)/(4*(wv^2-1));
    c4v := Rv + 2 + 4*tv;
    Av := x^2 + (Rv^3 + 4*Rv^2*tv + Rv - 8*Rv*tv + 4*tv)*x + Rv^4;
    Bv := c4v*x^2 + (Rv^2 + 4*Rv + 1 + 8*tv)*x + (2*Rv^2 + Rv + 4*tv);
    XR := -c4v*Rv;
    Atq := P![Q!co : co in Coefficients(c4v^2*Evaluate(Av, x/c4v))];
    Btq := P![Q!co : co in Coefficients(c4v*Evaluate(Bv, x/c4v))];
    return tv, c4v, XR, Atq, Btq;
end function;

procedure CheckComponent(name, gq, XR)
    printf "  %o polynomial=%o\n", name, gq;
    dsc := Discriminant(gq);
    printf "  %o discriminant=%o\n", name, dsc;
    if dsc eq 0 then
        printf "  %o FAIL: zero discriminant\n", name;
        return;
    end if;
    okd, sd := IsSquareQ(dsc);
    printf "  %o discriminant square over Q=%o", name, okd;
    if okd then printf " sqrt=%o", sd; end if;
    printf "\n";
    if okd then
        pass := true;
        for rt in Roots(gq) do
            val := XR - rt[1];
            okv, sv := IsSquareQ(val);
            printf "    root=%o XR-root=%o square=%o", rt[1], val, okv;
            if okv then printf " sqrt=%o", sv; end if;
            printf "\n";
            if not okv then pass := false; end if;
        end for;
        printf "  %o %o\n", name, pass select "PASS" else "FAIL";
    else
        K<th> := NumberField(gq);
        val := K!XR - th;
        ok := IsSquare(val);
        printf "  %o quadratic-field square=%o\n", name, ok;
        if not ok then
            printf "  %o FAIL: XR-theta is nonsquare in quadratic algebra\n", name;
        else
            printf "  %o PASS\n", name;
        end if;
    end if;
end procedure;

tests := [
    <Q!-8, Q!-8, "R=-8, lambda=-5/7, m=+/-28, w=-8">,
    <Q!-8, Q!8, "R=-8, lambda=-11/9, m=+/-36, w=8">,
    <Q!-25/4, Q!71/13, "R=-25/4 branch m=+/-105/4, w=71/13">,
    <Q!-25/4, Q!-71/13, "R=-25/4 branch m=+/-145/8, w=-71/13">
];

for T in tests do
    Rv := T[1];
    wv := T[2];
    printf "\n================ %o ================\n", T[3];
    tv, c4v, XR, Atq, Btq := FamilyData(Rv, wv);
    printf "t=%o c4=%o XR=%o\n", tv, c4v, XR;
    CheckComponent("A", Atq, XR);
    CheckComponent("B", Btq, XR);
end for;

quit;
