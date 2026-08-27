// Lane 8 (2026-07-25, third session): exact verification of contact-9 / contact-7
// survivors of the Z/63 sieve (code/claude_ov_lane8_c63sieve.c).
//
//   contact-9 (1 param a): h9 = 1 - 9/2 x + 63/8 x^2 - 105/16 x^3 + a x^4
//                          f  = (h9^2 + (x-1)^9) / x^4         (monic quintic)
//   contact-7 (2 params) : h7 = 1 - 7/2 x + a x^2 + b x^3
//                          f  = (h7^2 + (x-1)^7) / x^2         (monic quintic)
// The marked class is [P - infty] with P = (1, h(1)); it has order 9 resp. 7.
//
// TRAP (repo-wide): TorsionSubgroup needs an INTEGRAL model.  We scale
// x -> X/m^2, y -> Y/m^5 so that m^10 f(X/m^2) is integral.
//
// usage: code/claude_magma_slot.sh -b CHART:=9 AS:="201289/2665" code/claude_ov_lane8_c63verify.m
SetColumns(0);
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
if not assigned CHART then CHART := 9; elif Type(CHART) eq MonStgElt then CHART := StringToInteger(CHART); end if;

QQ := Rationals();
R<x> := PolynomialRing(QQ);

function c9poly(a)
    h := 1 - 9/2*x + 63/8*x^2 - 105/16*x^3 + a*x^4;
    g := h^2 + (x-1)^9;
    return &+[ Coefficient(g,i)*x^(i-4) : i in [4..9] ], h;
end function;

function c7poly(a,b)
    h := 1 - 7/2*x + a*x^2 + b*x^3;
    g := h^2 + (x-1)^7;
    return &+[ Coefficient(g,i)*x^(i-2) : i in [2..7] ], h;
end function;

// integral model of y^2 = f (deg 5): X = m^2 x, Y = m^5 y
function IntegralQuintic(f)
    d := LCM([Denominator(c) : c in Coefficients(f)]);
    m := 1;
    for t in Factorization(d) do
        m *:= t[1]^(Ceiling(t[2]/2));
    end for;
    repeat
        g := R![ m^(10-2*i)*Coefficient(f,i) : i in [0..5] ];
        ok := &and[ IsIntegral(c) : c in Coefficients(g) ];
        if not ok then m *:= 2; end if;
    until ok;
    return R![Integers()!c : c in Coefficients(g)], m;
end function;

as := [];
if assigned AS then as := [ StringToRational(s) : s in Split(AS, ",") ]; end if;
bs := [];
if assigned BS then bs := [ StringToRational(s) : s in Split(BS, ",") ]; end if;

n := CHART eq 9 select #as else Min(#as,#bs);
printf "CHART %o : %o survivor(s) to verify\n", CHART, n;

for i in [1..n] do
    a := as[i];
    if CHART eq 9 then
        f, h := c9poly(a);
        printf "\n=========== contact-9  a = %o\n", a;
    else
        b := bs[i];
        f, h := c7poly(a,b);
        printf "\n=========== contact-7  a = %o  b = %o\n", a, b;
    end if;
    printf "  f = %o\n", f;
    if Degree(f) ne 5 or not IsSquarefree(f) then
        printf "  DEGENERATE (deg %o, squarefree %o) -- skip\n", Degree(f), IsSquarefree(f);
        continue;
    end if;
    g, m := IntegralQuintic(f);
    printf "  integral model (X = %o^2 x): %o\n", m, g;
    C := HyperellipticCurve(g);
    printf "  genus %o   disc factored %o\n", Genus(C), Factorization(Integers()!Discriminant(g));
    J := Jacobian(C);
    T := TorsionSubgroup(J);
    printf "  TORSION invariants: %o   order %o\n", Invariants(T), #T;
    if #T mod 63 eq 0 then
        printf "  *** ORDER 63 DIVIDES THE TORSION ***\n";
    else
        printf "  --- 63 does NOT divide the torsion: FALSE POSITIVE of the sieve\n";
    end if;
    // marked class order (P = (1, h(1)) on the unscaled model -> (m^2, m^5 h(1)))
    v := Evaluate(h,1);
    ok := true;
    try
        P := C ! [m^2, m^5*v, 1];
        D := J ! (P - PointsAtInfinity(C)[1]);
        printf "  marked class [P-infty] order: %o\n", Order(D);
    catch e
        printf "  marked class: could not place P (%o)\n", e`Object;
    end try;
end for;
printf "SEARCH_DONE c63verify\n";
quit;
