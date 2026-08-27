// lane_goodbeta.m — scan X1(N) (N in {8,10,12}, class-(b) points) for curves
// whose beta = x(g1) - e_conj is a SQUARE in K = Q(sqrt(d)).  Such curves are
// one-sided halving carriers for MIXED gluings (validated delta machinery,
// see lane_1212.m):
//   t in GOOD10  +  any Delta-matched X1(12)-partner  ->  [2,60]-candidate
//        (order 120; psi~delta_E(T_rat) = delta_F(T'_rat) with the N=12 side
//        trivial and delta(T_rat^{(10)}) == delta(g1) == (alpha, beta))
//   u in GOOD12  +  any Delta-matched class-(b) X1(8)-partner -> [4,24]-type
//        candidate (order 96; halving im(2g,3g') with delta(2g) trivial)
//   two GOOD12 curves with equal d -> immediate [12,12] pair.
// Usage: cd product/code && magma -b lane_goodbeta.m > ../logs/lane_goodbeta.log
//   optional: H:=<int> (default 150)
SetColumns(0);
if not assigned H then H := 150; elif Type(H) eq MonStgElt then H := StringToInteger(H); end if;
SetMemoryLimit(4*10^9);

RQ := Rationals();
QT<t> := FunctionField(RQ);
RQx<x> := PolynomialRing(RQ);

BC := AssociativeArray();
BC[8]  := [* (2*t-1)*(t-1), (2*t-1)*(t-1)/t *];
BC[10] := [* t^3*(t-1)*(2*t-1)/(t^2-3*t+1)^2, -t*(t-1)*(2*t-1)/(t^2-3*t+1) *];
BC[12] := [* t*(2*t-1)*(3*t^2-3*t+1)*(2*t^2-2*t+1)/(t-1)^4, -t*(2*t-1)*(3*t^2-3*t+1)/(t-1)^3 *];

function SFrat(xv)
    n := Numerator(xv)*Denominator(xv);
    s := Sign(n); n := Abs(n);
    a := SquarefreeFactorization(n);
    return s*a;
end function;

function BetaData(N, tv, pair)
    ok := true; bv := 0; cv := 0;
    try bv := Evaluate(pair[1], tv); cv := Evaluate(pair[2], tv); catch e ok := false; end try;
    if not ok or bv eq 0 then return false, 0; end if;
    E := 0;
    try E := EllipticCurve([1-cv, -bv, -bv, 0, 0]); catch e return false, 0; end try;
    P := 0;
    try P := E![0,0]; catch e return false, 0; end try;
    if Order(P) ne N then return false, 0; end if;
    E2, phi := WeierstrassModel(E);
    f := RQx!(HyperellipticPolynomials(E2));
    rts := Roots(f);
    if #rts ne 1 then return false, 0; end if;
    erat := rts[1][1];
    q := f div (RQx.1 - erat);
    q1 := Coefficient(q,1); q0 := Coefficient(q,0);
    D := q1^2 - 4*q0;
    d := SFrat(D);
    okw, w := IsSquare(D/d);
    if not okw then return false, 0; end if;
    xg := (phi(P))[1];
    return true, < tv, d, xg + q1/2, w/2 >;   // beta = A + B*sqrt(d)
end function;

function HeightRats(H)
    S := [];
    for a in [1..H], b in [1..H] do
        if GCD(a,b) eq 1 then Append(~S, a/b); Append(~S, -a/b); end if;
    end for;
    return S;
end function;

for N in [8,10,12] do
    ngood := 0; ntot := 0;
    for tv in HeightRats(H) do
        ok, rec := BetaData(N, tv, BC[N]);
        if not ok then continue; end if;
        ntot +:= 1;
        d := rec[2];
        K := QuadraticField(d);
        if IsSquare(K!rec[3] + rec[4]*K.1) then
            ngood +:= 1;
            printf "GOODBETA N=%o t=%o d=%o\n", N, rec[1], d;
        end if;
    end for;
    printf "GOODSCAN N=%o: %o/%o good\n", N, ngood, ntot;
end for;
printf "GOODBETA_DONE\n";
quit;
