// lane_2266_sections.m — [2,2,6,6] structured attack, symbolic version:
// on HLP's (6)-surface (the C1&C2-locus, fibered over the u-line with
// section P0 = (y,z)=(-1,1) of infinite order), compute for small section
// combinations P = a*P1 + b*P2 the y-coordinate y_{a,b}(u) in Q(u), the
// induced t_{a,b}(u) = 3 + (u-3) y^2, and the second-gain condition
//   X_{a,b}:  w^2 = cls( (t-1)(t-9)(u-1)(u-9) )  over Q(u);
// then analyze each X (genus, small points, rank if genus 1) and report
// candidate u's.  Any nondegenerate point gives a [2,2,6,6] candidate pair.
// Usage: cd product/code && magma -b lane_2266_sections.m > ../logs/lane2266_sections.log
SetColumns(0);
SetMemoryLimit(4*10^9);

RQ := Rationals();
QU<u> := FunctionField(RQ);
Pt<T> := PolynomialRing(RQ);
Ru<Y> := PolynomialRing(QU);

function SFint(x)
    n := Numerator(x)*Denominator(x);
    s := Sign(n); n := Abs(n);
    a := SquarefreeFactorization(n);
    return s*a;
end function;

function SqClassFF(q)   // over Q(u), class rep in Pt
    error if q eq 0, "zero";
    N := Numerator(q); D := Denominator(q);
    P := Pt!N * Pt!D;
    fac := Factorization(P);
    if #fac eq 0 then unit := P; else unit := P div &*[ f[1]^f[2] : f in fac ]; end if;
    error if Degree(unit) ne 0, "unit not constant";
    c := RQ!Coefficient(unit, 0);
    cls := Pt!SFint(c);
    for f in fac do if IsOdd(f[2]) then cls *:= f[1]; end if; end for;
    return cls;
end function;

quart := (u+3)*(u-5)*((u-3)*Y^2+6)*((u-3)*Y^2-2);
C := HyperellipticCurve(quart);
D0 := (u+3)*(u-5);
base := C![1, D0];
E, mE := EllipticCurve(C, base);
printf "generic fiber E/Q(u) built\n";
P1 := mE(C![-1, D0]);
P2 := mE(C![1, -D0]);
P3 := mE(C![-1, -D0]);
printf "section points mapped\n";
mEi := Inverse(mE);

// enumerate small combinations
cand := [* *];
for a in [0..4] do
    for b in [-2..2] do
        if a eq 0 and b le 0 then continue; end if;
        P := a*P1 + b*P2;
        if P eq E!0 then continue; end if;
        Q := 0; ok := true;
        try Q := mEi(P); catch e ok := false; end try;
        if not ok then continue; end if;
        if Q[3] eq 0 then continue; end if;
        yv := Q[1]/Q[3];
        // t(u) and the condition class
        tv := 3 + (u-3)*yv^2;
        num := (tv-1)*(tv-9)*(u-1)*(u-9);
        if num eq 0 then continue; end if;
        cls := SqClassFF(num);
        dg := Degree(cls);
        printf "== section (a,b)=(%o,%o): deg(y)=%o, condition class degree %o ==\n",
            a, b, Degree(Numerator(yv)), dg;
        if dg eq 0 then
            printf "   CONSTANT class %o: condition %o IDENTICALLY!\n", cls,
                IsSquare(RQ!Coefficient(cls,0)) select "HOLDS" else "fails";
            if IsSquare(RQ!Coefficient(cls,0)) then
                printf "   PARAMETRIC FAMILY: t(u) = %o\n", tv;
            end if;
            continue;
        end if;
        if dg le 2 then
            printf "   low degree class: %o (conic-like; note)\n", cls;
            continue;
        end if;
        X := HyperellipticCurve(cls);
        gX := Genus(X);
        pts := Points(X : Bound := 1000);
        printf "   X: w^2 = %o, genus %o, small pts %o\n", cls, gX, #pts;
        for PP in pts do
            if PP[3] eq 0 then continue; end if;
            uv := PP[1]/PP[3];
            if uv in {RQ|3,-3,1,5,9} then continue; end if;
            tvv := Evaluate(Numerator(tv), uv)/Evaluate(Denominator(tv), uv);
            if tvv in {RQ|3,-3,1,5,9} or tvv eq uv then continue; end if;
            // verify all three conditions numerically
            c1 := (tvv-3)*(uv-3); c2 := (tvv+3)*(tvv-5)*(uv+3)*(uv-5); c3 := (tvv-1)*(tvv-9)*(uv-1)*(uv-9);
            if c1 eq 0 or c2 eq 0 or c3 eq 0 then continue; end if;
            if IsSquare(c1) and IsSquare(c2) and IsSquare(c3) then
                printf "   SECSURV t=%o u=%o (a,b)=(%o,%o)\n", tvv, uv, a, b;
            end if;
        end for;
        if gX eq 1 and #pts gt 0 then
            EX, mX := EllipticCurve(X, Rep(pts));
            EXm, mmin := MinimalModel(EX);
            GX, mwX := MordellWeilGroup(EXm);
            printf "   X rank data: MW %o (E=%o)\n", Invariants(GX), aInvariants(EXm);
            // enumerate MW points for candidates
            mXi := Inverse(mX); mmi := Inverse(mmin);
            ng := Ngens(GX);
            ranges := [ Order(GX.i) eq 0 select [-8..8] else [0..Order(GX.i)-1] : i in [1..ng] ];
            CP := ng eq 0 select [ [Integers()|] ] else [ [c[i] : i in [1..ng]] : c in CartesianProduct(ranges) ];
            for tup in CP do
                g := GX!0;
                for i in [1..ng] do g +:= tup[i]*GX.i; end for;
                PP := 0; okq := true;
                try PP := mXi(mmi(mwX(g))); catch e okq := false; end try;
                if not okq or PP[3] eq 0 then continue; end if;
                uv := PP[1]/PP[3];
                if uv in {RQ|3,-3,1,5,9} then continue; end if;
                if Max(Abs(Numerator(uv)),Abs(Denominator(uv))) gt 10^30 then continue; end if;
                tvv := Evaluate(Numerator(tv), uv)/Evaluate(Denominator(tv), uv);
                if tvv in {RQ|3,-3,1,5,9} or tvv eq uv then continue; end if;
                c1 := (tvv-3)*(uv-3); c2 := (tvv+3)*(tvv-5)*(uv+3)*(uv-5); c3 := (tvv-1)*(tvv-9)*(uv-1)*(uv-9);
                if c1 eq 0 or c2 eq 0 or c3 eq 0 then continue; end if;
                if IsSquare(c1) and IsSquare(c2) and IsSquare(c3) then
                    printf "   SECSURV t=%o u=%o (a,b)=(%o,%o) [mw]\n", tvv, uv, a, b;
                end if;
            end for;
        end if;
    end for;
end for;
printf "SECTIONS_DONE\n";
quit;
