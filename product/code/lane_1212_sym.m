// lane_1212_sym.m — symbolic analysis of the diagonal (self-gluing) halving
// condition for the X1(N) families, N in {8,10,12}:
//   alpha(t) = x(g1) - e_rat on the y^2 = cubic model;  Delta(t) = disc of the
//   conjugate quadratic factor.  Diagonal condition: alpha == 1 or Delta mod
//   squares.  This script computes alpha, Delta mod squares over Q(t), builds
//   the two condition curves  w^2 = cls(alpha)  and  w^2 = cls(alpha*Delta),
//   reports genus / rank / points, and enumerates MW points where possible
//   (candidate t-values printed for the numeric funnel).
// Usage: cd product/code && magma -b lane_1212_sym.m > ../logs/lane1212_sym.log
SetColumns(0);
SetMemoryLimit(4*10^9);

RQ := Rationals();
QT<t> := FunctionField(RQ);
Pt<T> := PolynomialRing(RQ);
Rx<x> := PolynomialRing(QT);

function SFint(v)
    n := Numerator(v)*Denominator(v);
    s := Sign(n); n := Abs(n);
    a := SquarefreeFactorization(n);
    return s*a;
end function;

function SqClassFF(q)
    error if q eq 0, "SqClassFF: zero input";
    N := Numerator(q); D := Denominator(q);
    P := Pt!N * Pt!D;
    fac := Factorization(P);
    if #fac eq 0 then unit := P;
    else unit := P div &*[ f[1]^f[2] : f in fac ]; end if;
    error if Degree(unit) ne 0, "unit not constant";
    c := RQ!Coefficient(unit, 0);
    cls := Pt!SFint(c);
    for f in fac do if IsOdd(f[2]) then cls *:= f[1]; end if; end for;
    return cls;
end function;

// verified Kubert variants (from the runtime verification in lane_1212.m)
BC := AssociativeArray();
BC[8]  := [* (2*t-1)*(t-1), (2*t-1)*(t-1)/t *];
BC[10] := [* t^3*(t-1)*(2*t-1)/(t^2-3*t+1)^2, -t*(t-1)*(2*t-1)/(t^2-3*t+1) *];
BC[12] := [* t*(2*t-1)*(3*t^2-3*t+1)*(2*t^2-2*t+1)/(t-1)^4, -t*(2*t-1)*(3*t^2-3*t+1)/(t-1)^3 *];

for N in [8,10,12] do
    printf "===== N = %o =====\n", N;
    bv := BC[N][1]; cv := BC[N][2];
    E := EllipticCurve([QT| 1-cv, -bv, -bv, 0, 0]);
    P := E![0,0];
    // sanity: N*P = 0
    error if not (N*P eq E!0), "N*P != 0";
    T2 := (N div 2)*P;
    E2, phi := WeierstrassModel(E);
    f := Rx![ Coefficient(HyperellipticPolynomials(E2), i) : i in [0..3] ];
    rts := Roots(f);
    error if #rts ne 1, Sprintf("expected 1 rational root over Q(t), got %o", #rts);
    erat := rts[1][1];
    q := f div (Rx.1 - erat);
    q1 := Coefficient(q,1); q0 := Coefficient(q,0);
    Delta := q1^2 - 4*q0;
    xg := (phi(P))[1];
    alpha := xg - erat;
    // consistency: x(6P)/x(5P)... the rational 2-torsion x-coord equals erat
    xT := (phi(T2))[1];
    error if xT ne erat, "e_rat mismatch with (N/2)g1";
    clsA  := SqClassFF(alpha);
    clsD  := SqClassFF(Delta);
    clsAD := SqClassFF(alpha*Delta);
    printf "alpha      = %o\n", alpha;
    printf "cls(alpha) = %o\n", clsA;
    printf "cls(Delta) = %o\n", clsD;
    printf "cls(alpha*Delta) = %o\n", clsAD;
    // checksum: for 4|N, kappa-class: delta(T_rat) trivial was verified
    // numerically on the pools; here only alpha matters.
    for CASE in [ <"untwisted", clsA>, <"twisted", clsAD> ] do
        cls := CASE[2];
        den := LCM([ Denominator(c) : c in Coefficients(cls) ]);
        cls := Pt!(cls*den^2);   // same condition curve, integral model
        printf "-- %o condition curve: w^2 = %o\n", CASE[1], cls;
        dg := Degree(cls);
        if dg eq 0 then
            printf "   constant class: %o (condition %o)\n", cls,
                IsSquare(RQ!Coefficient(cls,0)) select "ALWAYS" else "NEVER";
            continue;
        end if;
        if dg le 2 then
            printf "   degree <= 2 (conic-like): %o -- handle manually\n", cls;
            continue;
        end if;
        C := HyperellipticCurve(cls);
        g := Genus(C);
        printf "   genus %o\n", g;
        pts := Points(C : Bound := 2000);
        printf "   small points (bound 2000): %o\n", pts;
        if g eq 1 then
            if #pts eq 0 then
                printf "   genus 1, no small points: possible nontrivial torsor\n";
                continue;
            end if;
            EE, mE := EllipticCurve(C, Rep(pts));
            EEm, mmin := MinimalModel(EE);
            printf "   E = %o, conductor %o\n", aInvariants(EEm), Conductor(EEm);
            G, mw := MordellWeilGroup(EEm);
            printf "   MW %o\n", Invariants(G);
            // enumerate
            mEi := Inverse(mE);
            mmi := Inverse(mmin);
            ng := Ngens(G);
            NB := 10;
            ranges := [ Order(G.i) eq 0 select [-NB..NB] else [0..Order(G.i)-1] : i in [1..ng] ];
            CP := ng eq 0 select [ [Integers()|] ] else [ [c[i] : i in [1..ng]] : c in CartesianProduct(ranges) ];
            for tup in CP do
                gg := G!0;
                for i in [1..ng] do gg +:= tup[i]*G.i; end for;
                Q := mw(gg);
                okq := true; PP := 0;
                try PP := mEi(mmi(Q)); catch e okq := false; end try;
                if not okq or PP[3] eq 0 then continue; end if;
                tv := PP[1]/PP[3];
                if Max(Abs(Numerator(tv)),Abs(Denominator(tv))) gt 10^18 then continue; end if;
                printf "   CANDT N=%o case=%o t=%o\n", N, CASE[1], tv;
            end for;
        elif g ge 2 then
            printf "   genus %o >= 2: finitely many points; listed above\n", g;
            for PP in pts do
                if PP[3] eq 0 then continue; end if;
                printf "   CANDT N=%o case=%o t=%o\n", N, CASE[1], PP[1]/PP[3];
            end for;
        end if;
    end for;
end for;
printf "SYM_DONE\n";
quit;
