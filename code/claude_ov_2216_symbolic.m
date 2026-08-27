//////////////////////////////////////////////////////////////////////
//  claude_ov_2216_symbolic.m
//
//  Symbolic derivation of the x-T (delta) criteria for the marked
//  order-8 class Q and for its seven twists Q+T, T in J[2](Q)\{0},
//  on the odd M_1(8,2,2) chart, over the function field Q(u,v).
//
//  Odd model:  f = ((1-u)X+1)((1-v)X+1)((u+v+2)X+1) * qt(X)
//              qt = -X^2 + B X - C,
//              B = u^2 v - u^2 + u v^2 - u - v^2 - v - 2,
//              C = u^2 + u v + v^2 + u + v + 1.
//  Rational roots  r1 = 1/(u-1), r2 = 1/(v-1), r3 = -1/(u+v+2)
//  correspond to the sextic-chart roots  u, v, w = -(u+v+1)
//  under X = 1/(x-1)   (the fourth sextic root x=1 goes to infinity).
//
//  delta' conventions are those of code/claude_ov_2216_delta.m
//  (validated 372/372 against IsDivisibleBy).
//
//  Halving criterion (both mod Q*, so only ratios matter):
//      Q+T is 2-divisible  <=>  exists d in Q* with
//         delta'_1, delta'_2, delta'_3  all in d*(Q*)^2   and
//         delta'_K in d*(K*)^2 .
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(4*10^9);

FF<u,v> := RationalFunctionField(Rationals(), 2);
R := PolynomialRing(Rationals(), 2);
PX<X> := PolynomialRing(FF);

B := u^2*v - u^2 + u*v^2 - u - v^2 - v - 2;
C := u^2 + u*v + v^2 + u + v + 1;
qt  := -X^2 + B*X - C;
qtm := X^2 - B*X + C;                     // monic
f   := ((1-u)*X + 1)*((1-v)*X + 1)*((u+v+2)*X + 1)*qt;
c   := LeadingCoefficient(f);

r := [ 1/(u-1), 1/(v-1), -1/(u+v+2) ];
lin := [ X - r[i] : i in [1..3] ];

// sanity: f = c * prod(lin) * qtm
assert f eq c*&*lin*qtm;

// ---- squarefree normal form of an element of Q(u,v), modulo squares
function SqFree(g)
    if g eq 0 then return R!0; end if;
    nu := R!Numerator(g);  de := R!Denominator(g);
    h := nu*de;                                    // same square class
    fa := Factorization(h);
    unit := h;
    for t in fa do unit := unit div t[1]^t[2]; end for;   // dropped constant
    out := R!1;
    for t in fa do if IsOdd(t[2]) then out *:= t[1]; end if; end for;
    // put the constant back, reduced mod squares
    cst := Rationals()!unit;
    n := Numerator(cst); d := Denominator(cst);
    m := n*d; s := Sign(m); m := AbsoluteValue(m); rr := 1;
    for pe in Factorization(m) do
        if IsOdd(pe[2]) then rr *:= pe[1]; end if;
    end for;
    return (s*rr)*out;
end function;

// Print the squarefree class as  const * (p1)(p2)...  with every factor
// made PRIMITIVE-INTEGRAL (leading coeff > 0) and the leftover rational
// constant shown explicitly.  Magma's Factorization drops the unit, so
// we recompute it (lab rule #1) -- a dropped 2 changes the square class.
function ShowQ(g)
    sf := SqFree(g);
    if sf eq 0 then return "0"; end if;
    fa := Factorization(sf);
    unit := sf;
    for t in fa do unit := unit div t[1]^t[2]; end for;
    cst := Rationals()!unit;
    s := "";
    for t in fa do
        p := t[1];
        // make p primitive integral with positive leading coefficient
        den := LCM([ Denominator(cc) : cc in Coefficients(p) ]);
        pp := den*p;
        cont := GCD([ Numerator(cc) : cc in Coefficients(pp) ]);
        pp := pp div cont;
        if LeadingCoefficient(pp) lt 0 then pp := -pp; cst *:= -1; end if;
        cst *:= (Rationals()!cont/den)^t[2];
        for e in [1..t[2]] do s := s cat "(" cat Sprint(pp) cat ")"; end for;
    end for;
    // reduce the constant mod squares
    n := Numerator(cst); d := Denominator(cst);
    m := n*d; sg := Sign(m); m := AbsoluteValue(m); rr := 1;
    for pe in Factorization(m) do
        if IsOdd(pe[2]) then rr *:= pe[1]; end if;
    end for;
    cs := sg*rr;
    if cs eq 1 and s ne "" then return s; end if;
    return Sprint(cs) cat s;
end function;

// ---- K-component: reduce a polynomial mod qtm to a + b*X
function ModQ(g)
    return g mod qtm;
end function;

// delta' of a coprime class with Mumford u-polynomial uD
function DeltaCop(uD)
    return [* Evaluate(uD, r[1]), Evaluate(uD, r[2]), Evaluate(uD, r[3]), ModQ(uD) *];
end function;

// delta' of the 2-torsion class of the even finite subset encoded by
// the monic divisor uS | f:   a_j notin S -> uS(a_j),  a_j in S -> -(f/uS)(a_j)
function DeltaTwo(uS)
    cof := f div uS;
    assert uS*cof eq f;
    out := [* *];
    for i in [1..3] do
        if Evaluate(uS, r[i]) eq 0 then Append(~out, -Evaluate(cof, r[i]));
        else Append(~out, Evaluate(uS, r[i])); end if;
    end for;
    vS := ModQ(uS);
    if vS eq 0 then Append(~out, -ModQ(cof)); else Append(~out, vS); end if;
    return out;
end function;

function DMul(a, b)
    return [* a[1]*b[1], a[2]*b[2], a[3]*b[3], ModQ(a[4]*b[4]) *];
end function;

// Norm of  a + b*theta  from K = Q(u,v)[X]/(X^2 - B X + C) down to Q(u,v):
//   N(a+b*th) = a^2 + a*b*(th+th') + b^2*th*th' = a^2 + a*b*B + b^2*C.
// The field condition  delta'_1 * delta'_K in (K*)^2  forces
//   N(delta'_1 * delta'_K) = delta'_1^2 * N(delta'_K) in (Q*)^2,
// i.e. the purely RATIONAL necessary condition   N(delta'_K) = square.
function NormK(g)
    h := ModQ(g);
    a := Coefficient(h, 0); b := Coefficient(h, 1);
    return a^2 + a*b*B + b^2*C;
end function;

// ---------------------------------------------------------------
print "==== delta'(Q),  u_Q = X+1 ====";
dQ := DeltaCop(X+1);
for i in [1..3] do
    printf "  delta'_%o(Q) = %o     [sqfree] %o\n", i, dQ[i], ShowQ(dQ[i]);
end for;
printf "  delta'_K(Q) = %o\n", dQ[4];

// the seven nonzero 2-torsion classes
labels := [];  subs := [];
for i in [1..3] do for j in [i+1..3] do
    Append(~subs, lin[i]*lin[j]);
    Append(~labels, Sprintf("A_r%or%o  (S={r%o,r%o})", i, j, i, j));
end for; end for;
Append(~subs, qtm);  Append(~labels, "B_quad  (S={th,th'})");
for k in [1..3] do
    Append(~subs, qtm * &*[ lin[j] : j in [1..3] | j ne k ]);
    Append(~labels, Sprintf("C_r%o    (S=compl{r%o,inf}, class [(r%o,0)-inf])", k, k, k));
end for;

print "";
print "==== delta'(Q+T) for the seven twists  (square classes as factorizations) ====";
for k in [1..#subs] do
    dT := DeltaTwo(subs[k]);
    dQT := DMul(dQ, dT);
    printf "\n---- T = %o\n", labels[k];
    for i in [1..3] do
        printf "   d%o : %o\n", i, ShowQ(dQT[i]);
    end for;
    // ratios that must be squares (mod Q*): d1*d2, d1*d3
    printf "   d1*d2 : %o\n", ShowQ(dQT[1]*dQT[2]);
    printf "   d1*d3 : %o\n", ShowQ(dQT[1]*dQT[3]);
    printf "   dK   : %o\n", dQT[4];
    printf "   N(dK) [MUST BE SQUARE] : %o\n", ShowQ(NormK(dQT[4]));
end for;

print "";
print "==== untwisted ratios ====";
printf "   d1*d2 : %o\n", ShowQ(dQ[1]*dQ[2]);
printf "   d1*d3 : %o\n", ShowQ(dQ[1]*dQ[3]);
printf "   N(dK) [MUST BE SQUARE] : %o\n", ShowQ(NormK(dQ[4]));

print "SEARCH_DONE";
quit;
