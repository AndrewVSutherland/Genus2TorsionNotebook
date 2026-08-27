//////////////////////////////////////////////////////////////////////
// fable_z24_richelot_probe.m  (2026-07-18, Fable session)
//
// Depth-1 rational Richelot probe of the two known geometrically
// simple Z/24 curves (A(8) chart seeds from agent_a2_24_composite8x3):
//   curve A: r=5,   p=-5/2, t=-9/2
//   curve B: r=1/3, p=-1/9, t=-1
// Each codomain is (2,2)-isogenous to a simple Jacobian, hence simple.
// We look for 2-primary redistribution: [2,24], [4,12], [2,2,12], ...
// The odd part 3 survives every (2,2)-isogeny, so every codomain has
// 3 | #torsion automatically.
//
// Run:  magma -b code/fable_z24_richelot_probe.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(2*10^9);

Q := Rationals();
P<x> := PolynomialRing(Q);

// A(8) chart, copied VERBATIM from named-charts-reference /
// code/agent_a2_24_composite8x3.m (m=1 gauge).
function A8f(rv, pv, tv)
    e      := tv^2 - 2*pv*tv/rv;
    d      := e + 2*pv - rv^2;
    lambda := rv/tv;
    aa := rv^2 - lambda;
    bb := 2*rv*pv - 2*lambda*(pv + rv*tv) + 2*rv*lambda;
    cc := pv^2 + 2*pv*rv^2 - rv^4 - rv^3*tv - rv*pv^2/tv
          - lambda*(rv^2 + e) + 2*lambda*(rv*pv + rv^2*tv - 3*pv*tv + rv*tv^2);
    Qq := x^2 + d;
    qq := aa*x^2 + bb*x + cc;
    return qq*(Qq^2 + qq);
end function;

function IntegralSquareScale(f)
    den := LCM([Denominator(Coefficient(f,i)) : i in [0..Degree(f)]]);
    return P!(den^2*f), den;
end function;

procedure ProbeSeed(name, rv, pv, tv)
    f := A8f(rv, pv, tv);
    fI, den := IntegralSquareScale(f);
    C := HyperellipticCurve(fI);
    J := Jacobian(C);
    T, mp := TorsionSubgroup(J);
    printf "SEED %o r=%o p=%o t=%o torsion=%o\n", name, rv, pv, tv, Invariants(T);
    Rs := RichelotIsogenousSurfaces(J);
    printf "SEED %o codomains=%o\n", name, #Rs;
    for i in [1..#Rs] do
        S := Rs[i];
        tS := Type(S);
        if tS eq JacHyp or tS eq CrvHyp then
            D := tS eq JacHyp select Curve(S) else S;
            fD, hD := HyperellipticPolynomials(D);
            FD := hD eq 0 select P!fD else P!(hD^2 + 4*fD);
            FI, den := IntegralSquareScale(FD);
            TS := TorsionSubgroup(Jacobian(HyperellipticCurve(FI)));
            printf "COD %o.%o inv=%o\n", name, i, Invariants(TS);
        else
            // product of elliptic curves would contradict simplicity
            printf "COD %o.%o NONJACOBIAN type=%o\n", name, i, tS;
        end if;
    end for;
end procedure;

ProbeSeed("A", 5, -5/2, -9/2);
ProbeSeed("B", 1/3, -1/9, -1);
printf "PROBE_DONE\n";
quit;
