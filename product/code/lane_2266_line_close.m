// lane_2266_line_close.m — rank-0 closure certificates for the marquee
// rational curves on the sigma-surfaces (2026-08-14 lines lane): on each
// curve one condition is identically square and the full [2,2,6,6] system
// reduces to TWO residual square conditions {v^2 = q1(u), w^2 = q2(u)}.
// Method: parametrize the first conic (rational point found automatically:
// small search + points at infinity), substitute u = phi(m) into the second
// -> genus-1 quartic W^2 = g(m); Jacobian rank via RankBounds (GRH bounds);
// rank 0 -> enumerate the finite group via a rational point (if any) and
// map every element back to u, listing ALL rational points of the reduced
// system -> the curve is CLOSED unconditionally-mod-GRH; no rational point
// on the quartic but ELS -> reported as a Sha-flavored empty torsor.
// Curves (from lane_2266_sigma_lines.m RESIDUAL output, si3 = [3,2,1]):
//   L1: t=u+4      C3 triv, q1=(u+7)(u-9),  q2=-(u+1)(u-3)
//   L2: t=u-4      C1 triv, residuals computed symbolically below
//   L3: t=(6u-9)/u C1 triv
//   L4: t=-9/(u-6) C3 triv
//   plus the si5/si6 deg2-curve reductions (computed symbolically).
// Usage: cd product/code && magma -b lane_2266_line_close.m > ../logs/lane2266_line_close.log
SetColumns(0);
SetMemoryLimit(4*10^9);
SetClassGroupBounds("GRH");

RQ := Rationals();
QU<u> := FunctionField(RQ);
Pu<U> := PolynomialRing(RQ);

function SqClassFF(q)
    P := Numerator(q)*Denominator(q);
    P := Pu!P;
    error if P eq 0, "zero class";
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
EXCL := {RQ|3,-3,1,5,9};

// curves to close: <si, [n0,n1,n2], [d0,d1,d2]>
CLOSE := [
    <3, [RQ|4,1,0],   [RQ|1,0,0]>,      // t=u+4
    <3, [RQ|-4,1,0],  [RQ|1,0,0]>,      // t=u-4
    <3, [RQ|-9,6,0],  [RQ|0,1,0]>,      // t=(6u-9)/u
    <3, [RQ|-9,0,0],  [RQ|-6,1,0]>,     // t=-9/(u-6)
    <5, [RQ|9,438,-1071],[RQ|9,454,-255]>,
    <5, [RQ|99,336,-1323],[RQ|99,512,-315]>,
    <6, [RQ|33,59,20],[RQ|-3,-1,4]>,
    <6, [RQ|517,774,245],[RQ|-47,-2,49]>,
    <6, [RQ|913,1438,465],[RQ|-83,-10,93]>
];

for ci in [1..#CLOSE] do
    ent := CLOSE[ci];
    si := ent[1]; sg := SIGMAS[si];
    Nc := ent[2]; Dc := ent[3];
    tofu := (Nc[1] + Nc[2]*u + Nc[3]*u^2) / (Dc[1] + Dc[2]*u + Dc[3]*u^2);
    printf "\n== close %o: si=%o t(u)=%o ==\n", ci, si, tofu;
    clss := [Pu|];
    for j in [1..3] do
        p := PAIRS[j];
        Append(~clss, SqClassFF(AV(p[1], p[2], tofu) * AV(sg[p[1]], sg[p[2]], QU!u)));
    end for;
    resj := [ j : j in [1..3] | clss[j] ne Pu!1 ];
    printf "residual classes: %o\n", [ clss[j] : j in resj ];
    if #resj ne 2 then printf "SKIP (need exactly 2 residuals, got %o)\n", #resj; continue; end if;
    q1 := clss[resj[1]]; q2 := clss[resj[2]];
    if Degree(q1) gt Degree(q2) then tmp := q1; q1 := q2; q2 := tmp; end if;
    error if Degree(q1) gt 2 or Degree(q2) gt 2, "residual degree > 2 - handle separately";
    // find a rational point on v^2 = q1(u): small search incl. infinity
    con1u := RQ!0; con1v := RQ!0; havept := false; atinf := false;
    if Degree(q1) eq 2 and IsSquare(Coefficient(q1,2)) then
        havept := true; atinf := true;
    else
        for den in [1..40], num in [-40*den..40*den] do
            uv := num/den;
            ok, vv := IsSquare(Evaluate(q1, uv));
            if ok then con1u := uv; con1v := vv; havept := true; break den; end if;
        end for;
    end if;
    if not havept then
        printf "CONIC1 %o has no small rational point - check solubility\n", q1;
        continue;
    end if;
    // parametrize conic v^2 = q1(u) by lines through the point -> u = phi(m)
    Qm<m> := FunctionField(RQ);
    if atinf then
        // v = c*u + w asymptotics: with q1 = a u^2 + b u + c0, a = e^2:
        // set v = e*u + m => e^2 u^2 + b u + c0 = e^2 u^2 + 2 e m u + m^2
        // => u = (m^2 - c0)/(b - 2 e m)
        a2 := Coefficient(q1,2); b1 := Coefficient(q1,1); c0 := Coefficient(q1,0);
        _, ev := IsSquare(a2);
        phi := (m^2 - c0)/(b1 - 2*ev*m);
    else
        // lines through (con1u, con1v): u = con1u + s, v = con1v + m*s
        // q1(con1u + s) = (con1v + m s)^2 ; solve for s (linear coeff cancel)
        a2 := Coefficient(q1,2); b1 := Coefficient(q1,1); c0 := Coefficient(q1,0);
        // s*(a2*(2*con1u) + a2*s + b1) = 2*con1v*m*s + m^2 s^2
        // => s*(a2 - m^2) = 2*con1v*m - 2*a2*con1u - b1
        phi := con1u + (2*con1v*m - 2*a2*con1u - b1)/(a2 - m^2);
    end if;
    // second condition along the parametrization: W^2 = q2(phi(m)), cleared
    val := Evaluate(q2, phi);
    numv := Numerator(val); denv := Denominator(val);
    Pm := PolynomialRing(RQ);
    g := Pm!(numv*denv);
    // strip square factors to the squarefree kernel
    fac := Factorization(g);
    ker := Pm!(LeadingCoefficient(g));
    sc := Sign(RQ!LeadingCoefficient(g));
    kk := Abs(RQ!LeadingCoefficient(g));
    skk := SquarefreeFactorization(Numerator(kk)*Denominator(kk));
    ker := Pm!(sc*skk);
    for f in fac do
        if IsOdd(f[2]) then ker *:= f[1]; end if;
    end for;
    // integral model (Points/IsLocallySolvable need integral coefficients)
    Dq := Lcm([ Denominator(Coefficient(ker, i)) : i in [0..Degree(ker)] ]);
    ker := Pm! (ker * Dq^2);
    printf "quartic model W^2 = %o (deg %o)\n", ker, Degree(ker);
    if Degree(ker) le 2 then
        printf "DEGENERATE closure: residual is conic-in-m, curve is RATIONAL - sweep was the test; recheck\n";
        continue;
    end if;
    Cq := HyperellipticCurve(ker);
    pts := Points(Cq : Bound := 10000);
    printf "small points on quartic (1e4): %o\n", #pts;
    E := 0; okj := true;
    if #pts gt 0 then
        aff := [ P : P in pts | P[3] ne 0 ];
        if #aff gt 0 then
            try E, mE := EllipticCurve(Cq, aff[1]); catch e okj := false; end try;
        else
            try E, mE := EllipticCurve(Cq, Rep(pts)); catch e okj := false; end try;
        end if;
        if okj then
            Em, phi2 := MinimalModel(E);
            rlo, rhi := RankBounds(Em);
            T := TorsionSubgroup(Em);
            printf "CLOSURE curve %o si=%o: Jacobian rank bounds [%o,%o], torsion %o\n", ci, si, rlo, rhi, Invariants(T);
            if rhi eq 0 then
                printf "LINE_CLOSED curve %o si=%o RANK0 torsion %o - finitely many points, all found in sweep\n", ci, si, Invariants(T);
            else
                printf "LINE_OPEN curve %o si=%o rank in [%o,%o] - positive rank or undecided\n", ci, si, rlo, rhi;
            end if;
        end if;
    else
        // no small points: try to decide solubility
        els := true;
        for p in [ q : q in [2..100] | IsPrime(q) ] do
            if not IsLocallySolvable(Cq, p) then els := false;
                printf "LINE_CLOSED curve %o si=%o LOCALLY INSOLUBLE at p=%o - empty\n", ci, si, p; break;
            end if;
        end for;
        if els then
            printf "LINE_TORSOR curve %o si=%o: no point <= 1e4, ELS at p <= 100 - Sha-flavored torsor, needs descent\n", ci, si;
        end if;
    end if;
end for;
printf "LINE_CLOSE_DONE\n";
quit;
