/*
Exact symbolic models for the first orbit-12 radicand on the complete
Clebsch--Klein chart, and the fixed marked-ratio fibration.

For the marked class {r1^2,r2^2}, put

  G0 = -(r1^2-r3^2)(r1^2-r4^2)(r1^2-r5^2).

The script verifies the chart pullback, derives the two square equations on
a fixed q=r2/r1 fiber, and specializes them at q=-2.  It also prints the
local expansion used in the Q_23 obstruction on that special fiber.
*/

Qtm<t,m> := PolynomialRing(Rationals(), 2);
R := [
    1 + t*(t+2)*m,
    t*m*(m-t-2),
    -1 + m + t*(t+1)*m^2,
    1 + t - m - t*m^2,
    -(1+t)*(1+t*m^2)
];

assert &+R eq 0;
assert &+[z^3 : z in R] eq 0;

G0 := -&*[R[1]^2-R[k]^2 : k in [3..5]];
printf "ELKIES22210_ORBIT12_FIRST_COVER_FIBRATION\n";
printf "chart_G0_factorization %o\n", Factorization(G0);

/* Fix r1=1, r2=q, r3=x.  The two remaining roots u,v have sum S and
   product P.  Their discriminant is a square exactly when the CK tuple
   completes over Q.  The first radicand is symmetric in u,v. */
Qqx<q,x> := PolynomialRing(Rationals(), 2);
S := -(1+q+x);
P := (q+x)*(1+q)*(1+x)/(1+q+x);
Delta := S^2 - 4*P;
First := -(1-x^2)*((1+P)^2-S^2);

/* Multiplying numerator by denominator preserves rational squareclass. */
DeltaModel := Numerator(Delta)*Denominator(Delta);
FirstModel := Numerator(First)*Denominator(First);
printf "fixed_q_ck_discriminant_model %o\n", Factorization(DeltaModel);
printf "fixed_q_first_radicand_model %o\n", Factorization(FirstModel);

Dminus2 := Evaluate(DeltaModel, q, -2);
Fminus2 := Evaluate(FirstModel, q, -2);
printf "q_minus_2_ck_model %o\n", Factorization(Dminus2);
printf "q_minus_2_first_model %o\n", Factorization(Fminus2);

/* After deleting square factors, the two q=-2 quartics are

     V^2 = (x-1)(x^3+x^2-x-9),
     U^2 = -x(x-2)(x-1)(x+1).
*/
assert Dminus2 eq (x-1)*(x^3+x^2-x-9);
assert Fminus2 eq 4*(x-1)^2*(-x*(x-2)*(x-1)*(x+1));

/* Every mod-23 point satisfying all four square conditions is one of the
   six permutations of the same boundary triple. */
F := GF(23);
base_solutions := {};
for xf,uf in F do
    vf := 1-xf-uf;
    rf := [F!1,F!-2,xf,uf,vf];
    af := [z^2 : z in rf];
    gf := [
        -(af[1]-af[3])*(af[1]-af[4])*(af[1]-af[5]),
        (af[3]-af[2])*(af[1]-af[4])*(af[1]-af[5]),
        (af[4]-af[2])*(af[1]-af[3])*(af[1]-af[5]),
        (af[5]-af[2])*(af[1]-af[3])*(af[1]-af[4])
    ];
    if &+[z^3:z in rf] eq 0 and &and[IsSquare(z):z in gf] then
        Include(~base_solutions,<xf,uf,vf>);
    end if;
end for;
expected_solutions := {
    <F!0,F!2,F!-1>, <F!0,F!-1,F!2>,
    <F!2,F!0,F!-1>, <F!2,F!-1,F!0>,
    <F!-1,F!0,F!2>, <F!-1,F!2,F!0>
};
assert base_solutions eq expected_solutions;
printf "q_minus_2_mod23_cover_reductions %o\n", #base_solutions;

/* Local expansion at a forced F_23 boundary point.  Write
   (r3,r4,r5)=(A,2+B,-1-A-B). */
QAB<A,B> := PolynomialRing(Rationals(), 2);
xx := A;
uu := 2+B;
vv := -1-A-B;
CK := 1+(-2)^3+xx^3+uu^3+vv^3;
G_u := (uu^2-4)*(1-xx^2)*(1-vv^2);
printf "q_minus_2_boundary_ck %o\n", CK;
printf "q_minus_2_boundary_Gu %o\n", G_u;
assert MonomialCoefficient(CK,A) eq -3;
assert MonomialCoefficient(CK,B) eq 9;
assert MonomialCoefficient(G_u,A*B) eq -8;
assert MonomialCoefficient(G_u,B^2) eq -8;
/* A/B=3 in the first nonzero 23-adic digit, so G_u/B^2=-32 mod 23. */
assert (-32 mod 23) eq 14;
assert LegendreSymbol(Integers()!14, 23) eq -1;
printf "q_minus_2_Q23_leading_unit %o nonsquare true\n", -32 mod 23;
printf "DONE\n";
