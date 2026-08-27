// lane_2266_sigma_sections.m — symbolic section-multiple machinery on the
// five non-identity sigma-surfaces, using the UNIVERSAL DECK SECTIONS found
// 2026-08-14 (deck_derive.m + deck_sigma_match.m): each sigma-system is
// satisfied identically by exactly one Mobius deck map,
//   si2=[2,1,3]:(12)  t=(u+15)/(u-1)     si3=[3,2,1]:(13)  t=6-u
//   si4=[1,3,2]:(23)  t=(5u-9)/(u-5)     si5=[2,3,1]:(123) t=(u-21)/(u-5)
//   si6=[3,1,2]:(132) t=(5u-21)/(u-1)
// (the non-id analogue of HLP's (y,z)=(-1,1) universal section for sigma=id).
// Over Q(u) the C2-parametrized fiber is Y^2 = quart(S) with
//   d := 2*A(sg[1],sg[3])(u),  quart := cA*(S^2+6d)(S^2-2d),  cA := A(sg[1],sg[2])(u)
// (Org=1 organization; the test condition is C3 with partner cT := A(sg[2],sg[3])(u)).
// The deck section gives a base point (s0,Y0) in C(Q(u)); the sign-flips
// (-s0,Y0),(s0,-Y0) give independent-looking sections P1,P2 of E(Q(u)).
// For each multiple a*P1+b*P2 pull back, t(u) := 3 + S^2/d, and compute the
// class of the test condition in Q(u)*/squares:
//   degree 0 class 1  -> PARAMETRIC FAMILY of full [2,2,6,6] candidates
//                        (unless t(u) is itself one of the six deck maps ->
//                         DEGENERATE)
//   else              -> condition curve w^2 = cls(u): genus, small points,
//                        each rational point -> full-candidate test + report
// Mirrors lane_2266_sections.m (which closed the sigma=id case).
// Usage: cd product/code && magma -b lane_2266_sigma_sections.m > ../logs/lane2266_sigma_sections.log
//   optional: AMAX:=4 BMAX:=2 PtsB:=1000 MemGB:=4
SetColumns(0);
if not assigned AMAX then AMAX := 4; elif Type(AMAX) eq MonStgElt then AMAX := StringToInteger(AMAX); end if;
if not assigned BMAX then BMAX := 2; elif Type(BMAX) eq MonStgElt then BMAX := StringToInteger(BMAX); end if;
if not assigned PtsB then PtsB := 1000; elif Type(PtsB) eq MonStgElt then PtsB := StringToInteger(PtsB); end if;
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

RQ := Rationals();
QU<u> := FunctionField(RQ);
Pu<U> := PolynomialRing(RQ);
Ru<S> := PolynomialRing(QU);

function SqClassFF(q)  // squarefree class rep in Pu of nonzero q in Q(u)
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
DECKMAP := [ QU!u, (u+15)/(u-1), 6-u, (5*u-9)/(u-5), (u-21)/(u-5), (5*u-21)/(u-1) ];
EXCL := {RQ|3,-3,1,5,9};
HeightOf := func<q | Max(Abs(Numerator(q)), Abs(Denominator(q)))>;

for si in [2..6] do
    sg := SIGMAS[si];
    del := DECKMAP[si];
    printf "\n== si=%o sigma=%o deck section t=%o ==\n", si, sg, del;
    d  := 2*AV(sg[1], sg[3], QU!u);
    cA := AV(sg[1], sg[2], QU!u);
    cT := AV(sg[2], sg[3], QU!u);
    quart := cA*(S^2+6*d)*(S^2-2*d);
    // base point from the deck section
    oks, s0 := IsSquare((del-3)*d);
    error if not oks, "deck section fails C2 symbolically";
    oky, Y0 := IsSquare(Evaluate(quart, s0));
    error if not oky, "deck section fails C1 symbolically";
    printf "base s0=%o Y0=%o\n", s0, Y0;
    C := HyperellipticCurve(quart);
    base := C![s0, Y0];
    E, mE := EllipticCurve(C, base);
    P1 := mE(C![-s0, Y0]);
    P2 := mE(C![s0, -Y0]);
    mEi := Inverse(mE);
    printf "section orders: P1 %o, P2 %o\n", Order(P1), Order(P2);
    seenT := {};
    for a in [0..AMAX], b in [-BMAX..BMAX] do
        if a eq 0 and b le 0 then continue; end if;
        P := a*P1 + b*P2;
        okq := true; Q := 0;
        try Q := mEi(P); catch e okq := false; end try;
        if not okq or Q[3] eq 0 then continue; end if;
        sv := Q[1]/Q[3];
        tv := 3 + sv^2/d;
        if tv in seenT then continue; end if;
        Include(~seenT, tv);
        cls := SqClassFF(AV(2, 3, tv) * cT);  // test condition C3: pair {2,3}
        dg := Degree(cls);
        if dg eq 0 then
            if cls eq Pu!1 then
                isdeck := &or[ tv eq DECKMAP[k] : k in [1..6] ];
                if isdeck then
                    printf "(a,b)=(%o,%o): DEGENERATE deck family t(u)=%o\n", a, b, tv;
                else
                    printf "(a,b)=(%o,%o): PARAMETRIC FAMILY t(u) = %o  <<<<< FULL-CANDIDATE FAMILY\n", a, b, tv;
                end if;
            else
                printf "(a,b)=(%o,%o): CONSTANT nontrivial class %o (t(u)=%o)\n", a, b, cls, tv;
            end if;
            continue;
        end if;
        if dg le 2 then
            printf "(a,b)=(%o,%o): conic-like class deg %o: %o (t(u) deg %o)\n", a, b, dg, cls, Degree(Numerator(tv));
            // still search small points on w^2 = cls(u)
        end if;
        X := HyperellipticCurve(cls);
        gX := Genus(X);
        pts := [];
        try pts := Points(X : Bound := PtsB); catch e pts := []; end try;
        printf "(a,b)=(%o,%o): condition curve genus %o, deg %o, pts(%o)=%o\n", a, b, gX, dg, PtsB, #pts;
        for pt in pts do
            if pt[3] eq 0 then continue; end if;
            uv := pt[1]/pt[3];
            if uv in EXCL then continue; end if;
            // evaluate t(u) at uv; guard poles
            okev := true; tvv := RQ!0;
            try tvv := Evaluate(tv, uv); catch e okev := false; end try;
            if not okev or tvv in EXCL then continue; end if;
            // full triple test
            c1 := IsSquare(AV(1,2,tvv) * AV(sg[1],sg[2], uv));
            c2 := IsSquare(AV(1,3,tvv) * AV(sg[1],sg[3], uv));
            c3 := IsSquare(AV(2,3,tvv) * AV(sg[2],sg[3], uv));
            if c1 and c2 and c3 then
                printf "SECSURV si=%o t=%o u=%o (a,b)=(%o,%o) ht=%o\n", si, tvv, uv, a, b, HeightOf(tvv);
            end if;
        end for;
    end for;
end for;
printf "SIGMA_SECTIONS_DONE\n";
quit;
