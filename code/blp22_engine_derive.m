// BLP [2,22] engine, part A: derivation + validation of the order-11
// closure conditions t0, t1 for C_{a,b,c,d}: y^2 = F = R^2 - 4c^2 S^2,
// R = x^3 - x^2 + a x + b, S = x^2 + d  (Bernard-Leprevost-Pohst,
// Experiment. Math. 18 (2009), 65-70).
//
// D_inf has order 11 iff there is psi = yU - V (U monic cubic, V monic
// sextic) with F U^2 - V^2 = t S.  KEY LINEARIZATION (this session):
// with sqrtF = x^3 * sum_j s_j x^-j (s_0 = 1) the sqrt series of F at
// infinity, and V := polynomial part of U*sqrtF, one has
//   F U^2 - V^2 = rho (2V + rho),   rho := U*sqrtF - V = sum_m r_m x^-m,
// so deg <= 5 automatically, and killing the x^5, x^4, x^3 coefficients
// is EXACTLY r_1 = r_2 = r_3 = 0 with
//   r_m = s_{6+m} + u_2 s_{5+m} + u_1 s_{4+m} + u_0 s_{3+m}
// -- LINEAR in (u_0, u_1, u_2).  Then with G = F U^2 - V^2:
//   t := g_2,   t_1 := g_1,   t_0 := g_0 - t*d,
// and the order-11 locus is t_0 = t_1 = 0.
//
// VALIDATION: t_0 = t_1 = 0 at all valid BLP Table-1 tuples; both nonzero
// at perturbed tuples.  (Row C4 has a suspected typo in the paper's b and
// is EXPECTED to fail; row tC6 duplicates tC3.)
//
// Run: magma -b code/blp22_engine_derive.m > results/blp22_engine_derive.log

SetColumns(0);
SetSeed(1);
SetMemoryLimit(3*10^9);

// Derive <t0, t1, ok> over the field containing a,b,c,d.
// ok=false when the 3x3 linear system is singular (degenerate locus).
function DeriveT(a, b, c, d)
    K := Parent(a);
    P := PolynomialRing(K); x := P.1;
    R := x^3 - x^2 + a*x + b;
    S := x^2 + d;
    F := R^2 - 4*c^2*S^2;
    if Degree(F) ne 6 then return K!0, K!0, false; end if;
    fc := [Coefficient(F, i) : i in [0..6]];   // fc[i+1] = coeff of x^i
    // sqrt series: (1 + sum_{k>=1} s_k X^k)^2 = 1 + sum Ft_k X^k,
    // X = 1/x, Ft_k = coeff of x^{6-k} of F (monic), 0 for k > 6.
    s := [K| 1];   // s[j+1] = s_j
    for k in [1..9] do
        Ftk := k le 6 select fc[6-k+1] else K!0;
        acc := Ftk;
        for i in [1..k-1] do
            acc -:= s[i+1]*s[k-i+1];
        end for;
        Append(~s, acc/2);
    end for;
    // linear system r_m = 0, m = 1..3, unknowns u = (u0, u1, u2):
    // r_m = s_{6+m} + u2 s_{5+m} + u1 s_{4+m} + u0 s_{3+m}
    M := Matrix(K, 3, 3, [ [ s[3+m+1], s[4+m+1], s[5+m+1] ] : m in [1..3] ]);
    rhs := Vector(K, [ -s[6+m+1] : m in [1..3] ]);
    if Determinant(M) eq 0 then return K!0, K!0, false; end if;
    u := Solution(Transpose(M), rhs);   // u[1]=u0, u[2]=u1, u[3]=u2
    U := x^3 + u[3]*x^2 + u[2]*x + u[1];
    // V = polynomial part of U*sqrtF: V_n = sum_i u_i s_{i+3-n}, n=0..6
    uco := [u[1], u[2], u[3], K!1];   // u_0..u_3
    V := P!0;
    for n in [0..6] do
        vn := K!0;
        for i in [Max(0, n-3)..3] do
            j := i + 3 - n;
            if j ge 0 and j le 9 then vn +:= uco[i+1]*s[j+1]; end if;
        end for;
        V +:= vn*x^n;
    end for;
    G := F*U^2 - V^2;
    error if Degree(G) gt 5, "internal: G degree > 5";
    assert Coefficient(G,5) eq 0 and Coefficient(G,4) eq 0 and Coefficient(G,3) eq 0;
    t := Coefficient(G, 2);
    return Coefficient(G, 0) - t*d, Coefficient(G, 1), true;   // t0, t1
end function;

Q := Rationals();

rows := [
    <"C1",  -101/48,   -61/48,     1/4,   -5/12>,
    <"C2",  473/147,   -4013/343,  6/7,   207/49>,
    <"C3",  8/49,      -134/49,    3/7,   47/49>,
    <"C4",  1159/81,   261607/2187, 40/9, 13/27>,   // suspected b typo
    <"C5",  -1/13,     -191/2197,  8/13,  15/169>,
    <"C6",  -28/169,   103/2197,   3/13,  -4/169>,
    <"C7",  594/1805,  13348/34295, 8/19, -64/361>,
    <"C8",  208/867,   1338/4913,  5/17,  -39/289>,
    <"C9",  415/1089,  -2207/1089, 8/33,  119/121>,
    <"C10", 4989/2500, -13599/12500, 27/50, -81/250>,
    <"tC1", -3,        59,         4,     -7>,
    <"tC2", -163/1215, -367/3645,  2/3,   13/243>,
    <"tC3", -13/18,    71/6,       5/3,   -13/3>,
    <"tC4", -2287/27,  -1171/9,    10/3,  -323/3>,
    <"tC5", 121/147,   -141/343,   2/7,   15/49>,
    <"tC7", -1494/847, 19480/9317, 2/11,  -256/121>,
    <"tC8", 125/121,   -223/1331,  6/11,  29/121>,
    <"tC9", 187/361,   -649/6859,  6/19,  23/361>
];

npass := 0; nfail := 0;
for row in rows do
    t0, t1, ok := DeriveT(Q!row[2], Q!row[3], Q!row[4], Q!row[5]);
    vanish := ok and t0 eq 0 and t1 eq 0;
    printf "VALIDATE %o ok=%o t0_zero=%o t1_zero=%o\n",
           row[1], ok, ok and t0 eq 0, ok and t1 eq 0;
    if vanish then npass +:= 1; else nfail +:= 1;
        printf "  NONVANISHING %o: t0=%o t1=%o\n", row[1], t0, t1;
    end if;
end for;
printf "VANISH_PASS %o / %o (C4 expected to fail if its printed b is a typo)\n",
       npass, #rows;

// negative controls: perturbed tuples must NOT satisfy t0 = t1 = 0
nneg := 0;
for row in [ <"tC1p", -3, 59, 4, -8>, <"C5p", -1/13, -191/2197, 8/13, 16/169>,
             <"rand", 2/3, 1/5, 1/2, -3> ] do
    t0, t1, ok := DeriveT(Q!row[2], Q!row[3], Q!row[4], Q!row[5]);
    bad := ok and (t0 ne 0 or t1 ne 0);
    printf "NEGCTRL %o nonvanishing=%o\n", row[1], bad;
    if bad then nneg +:= 1; end if;
end for;
assert nneg eq 3;
assert npass ge 16;   // 17 valid distinct tuples (18 rows incl dup); C4 may fail
printf "BLP22_ENGINE_DERIVATION_VALIDATED\n";

// If C4 failed: solve for the correct b at C4's (a,c,d) as a rational root
// of t0(b), t1(b) over Q(b) -- recovers the paper's intended value.
Kb<bb> := RationalFunctionField(Q);
t0b, t1b, okb := DeriveT(Kb!(1159/81), bb, Kb!(40/9), Kb!(13/27));
if okb then
    n0 := Numerator(t0b); n1 := Numerator(t1b);
    g := GCD(PolynomialRing(Q)!n0, PolynomialRing(Q)!n1);
    printf "C4_B_RECOVERY gcd_degree=%o roots=%o\n", Degree(g), Roots(g);
end if;
quit;
