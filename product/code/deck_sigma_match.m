// deck_sigma_match.m — match each Mobius deck map delta of the E26 family to
// the sigma-systems it satisfies identically: for each delta and sigma test
// whether all three class conditions A_p(delta(u)) * A_sigma(p)(u) are
// squares in Q(u).  A full match = the graph t = delta(u) is a SECTION of the
// sigma-surface.  Pure symbolic, seconds.
// Usage: cd product/code && magma -b deck_sigma_match.m > ../logs/deck_sigma_match.log
SetColumns(0);
SetMemoryLimit(2*10^9);

RQ := Rationals();
QU<u> := FunctionField(RQ);
Pu<U> := PolynomialRing(RQ);

function SqClassFF(q)  // squarefree class rep of nonzero element of Q(u)
    P := Numerator(q)*Denominator(q);
    P := Pu!P;
    if P eq 0 then error "zero"; end if;
    fac := Factorization(P);
    unit := P div &*[ f[1]^f[2] : f in fac ];
    error if Degree(unit) ne 0, "unit not constant";
    c := RQ!unit;
    sc := Sign(Numerator(c)*Denominator(c));
    n := Abs(Numerator(c)*Denominator(c));
    sq := SquarefreeFactorization(n);
    r := Pu!(sc*sq);
    for f in fac do
        if IsOdd(f[2]) then r *:= f[1]; end if;
    end for;
    return r;
end function;

// ordered pair-class value A(a,b)(x) over Q(u)
function AV(a, b, x)
    if a lt b then
        if a eq 1 then
            return b eq 2 select (x+3)*(x-5) else 2*(x-3);
        else
            return -(x-1)*(x-9);
        end if;
    else
        return -AV(b, a, x);
    end if;
end function;

SIGMAS := [ [1,2,3],[2,1,3],[3,2,1],[1,3,2],[2,3,1],[3,1,2] ];
PAIRS  := [ [1,2], [1,3], [2,3] ];

DECKS := [* <"id: t=u", QU!u>,
           <"t=6-u", 6-u>,
           <"t=(5u-9)/(u-5)", (5*u-9)/(u-5)>,
           <"t=(u-21)/(u-5)", (u-21)/(u-5)>,
           <"t=(5u-21)/(u-1)", (5*u-21)/(u-1)>,
           <"t=(u+15)/(u-1)", (u+15)/(u-1)> *];

for D in DECKS do
    name := D[1]; del := D[2];
    for si in [1..6] do
        sg := SIGMAS[si];
        allok := true;
        conds := [];
        for j in [1..3] do
            p := PAIRS[j];
            cls := SqClassFF(AV(p[1], p[2], del) * AV(sg[p[1]], sg[p[2]], QU!u));
            Append(~conds, cls eq Pu!1);
            if cls ne Pu!1 then allok := false; end if;
        end for;
        if allok then
            printf "SECTION deck %o IS a section of sigma=%o (si=%o)\n", name, sg, si;
        elif &or conds then
            printf "partial: deck %o sigma=%o conds=%o\n", name, sg, conds;
        end if;
    end for;
end for;
printf "DECK_SIGMA_MATCH_DONE\n";
quit;
