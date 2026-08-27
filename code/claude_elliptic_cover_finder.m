// claude_elliptic_cover_finder.m -- numerical elliptic-quotient map finder
// used to DISCOVER the explicit cover certificates embedded in
// paper/scripts_and_data/verify_split_certificates.m (2026-08-24): analytic
// Jacobian + endomorphism idempotents -> elliptic sublattice -> g2/g3 and
// Weierstrass-P via q-series -> exact algebraic recognition (deg <= DEG)
// -> exact rational-function fit over K -> exact identity verification.
// The resulting certificates are independent of this finder.
// Usage: magma -b GC:=<sextic coeffs> [PREC:=200] [DEG:=2] [SMALLX:=1] claude_elliptic_cover_finder.m
// K-capable elliptic-quotient map finder: builds E from the lattice's own
// g2, g3 (recognized in a number field K of degree <= 2), evaluates the map
// via Weierstrass-P q-series, fits exactly over K, verifies exactly over K.
// Params: GC (sextic coeffs), PREC (default 200)
SetColumns(0);
SetMemoryLimit(8*10^9);
prec := assigned PREC select StringToInteger(PREC) else 200;
CC := ComplexField(prec);
RR := RealField(prec);
QQ := Rationals(); Z := Integers();
PxQ<X> := PolynomialRing(QQ);
Px<x> := PolynomialRing(CC);
gc := assigned GC select [ StringToInteger(t) : t in Split(GC, ",") ]
                  else [-12,4,-19,4,-6,-4,1];
gQ := PxQ!gc;
g := Px![ CC!c : c in Coefficients(gQ) ];

function EisWP(w1, w2, zs)
    // returns g2, g3 of Z w1 + Z w2 (Im(w1/w2) > 0) and [ wp(z) : z in zs ]
    CCl := Parent(w1); il := CCl.1;
    pr := Precision(CCl);
    tau := w1/w2;
    q := Exp(2*Pi(RealField(pr))*il*tau);
    N := Ceiling(pr/AbsoluteValue(Log(10, AbsoluteValue(q)))) + 5;
    E4 := CCl!1; E6 := CCl!1;
    for n in [1..N] do
        t := q^n/(1-q^n);
        E4 +:= 240*n^3*t;
        E6 -:= 504*n^5*t;
    end for;
    c := 2*Pi(RealField(pr))*il/w2;
    g2 := c^4*E4/12;
    g3 := -c^6*E6/216;
    wps := [ CCl | ];
    for z in zs do
        u := Exp(2*Pi(RealField(pr))*il*z/w2);
        S := CCl!1/12 + u/(1-u)^2;
        for n in [1..N] do
            S +:= u*q^n/(1-u*q^n)^2 + (q^n/u)/(1-q^n/u)^2 - 2*q^n/(1-q^n)^2;
        end for;
        Append(~wps, c^2*S);
    end for;
    return g2, g3, wps;
end function;

// algebraic recognition of degree <= 2 with verification
MAXDEG := assigned DEG select StringToInteger(DEG) else 2;
function RecAlg2(c, prec)
    for dg in [1..MAXDEG] do
        mp := PowerRelation(c, dg);
        if Degree(mp) ne dg then continue; end if;
        ht := Maximum([ AbsoluteValue(co) : co in Coefficients(mp) ]);
        if ht gt 10^(prec div 5) then continue; end if;
        val := AbsoluteValue(Evaluate(PolynomialRing(Parent(c))!mp, c));
        if val lt ht*Maximum(RealField(prec)!1, AbsoluteValue(c))^dg*10^(-(2*prec) div 3) then
            return true, PxQ!mp;
        end if;
    end for;
    return false, PxQ!0;
end function;

A := AnalyticJacobian(g);
M, fs := EndomorphismRing(A);
MQ := MatrixAlgebra(QQ,4);
Bm := Basis(M);
cands0 := [ MQ!b : b in Bm ];
for i in [1..#Bm], j in [i+1..#Bm] do
    Append(~cands0, MQ!Bm[i] + MQ!Bm[j]);
    Append(~cands0, MQ!Bm[i] - MQ!Bm[j]);
end for;
idems := [ MQ | ];
for t in cands0 do
    if IsScalar(t) then continue; end if;
    fac := Factorization(PxQ!MinimalPolynomial(t));
    lin := [ u[1] : u in fac | Degree(u[1]) eq 1 ];
    if #lin lt 2 then continue; end if;
    for a in [1..#lin], b in [1..#lin] do
        if a eq b then continue; end if;
        r1 := -Coefficient(lin[a],0); r2 := -Coefficient(lin[b],0);
        ee := (t - r2*One(MQ))/(r1-r2);
        if ee^2 ne ee or Rank(ee) ne 2 then continue; end if;
        if ee notin idems then Append(~idems, ee); end if;
    end for;
end for;
printf "%o candidate idempotents\n", #idems;
Pi0 := BigPeriodMatrix(A);

if assigned SMALLX then
    xs := [ QQ | (IsEven(k) select 1 else -1) * (QQ!k/29 + QQ!1/(k+30)) : k in [1..26] ];
else
    xs := [ QQ | (IsEven(k) select 1 else -1) * (QQ!k/7 + QQ!3/(5+k)) : k in [1..40] ];
end if;
xs := [ xi : xi in xs | Evaluate(gQ, xi) ne 0 ];

for e in idems do
    eC := Matrix(CC,4,4,[ CC!c : c in Eltseq(e) ]);
    PiE := Pi0*eC;
    S12 := Matrix(CC,2,2,[Pi0[1][1],Pi0[1][2],Pi0[2][1],Pi0[2][2]]);
    T12 := Matrix(CC,2,2,[PiE[1][1],PiE[1][2],PiE[2][1],PiE[2][2]]);
    eA := T12*S12^-1;
    if &+[ AbsoluteValue(c) : c in Eltseq(eA*Pi0 - PiE) ] gt 10^(-prec div 2) then continue; end if;
    den := LCM([ Denominator(c) : c in Eltseq(e) ]);
    EI := Matrix(Z,4,4,[ Z!(den*c) : c in Eltseq(e) ]);
    H := HermiteForm(Transpose(EI));
    if Rank(H) ne 2 then continue; end if;
    bi := [ [ H[r][i]/den : i in [1..4] ] : r in [1..2] ];
    col1 := Pi0*Matrix(CC,4,1,[ CC!bi[1][i] : i in [1..4] ]);
    col2 := Pi0*Matrix(CC,4,1,[ CC!bi[2][i] : i in [1..4] ]);
    idx := AbsoluteValue(col1[1][1]) ge AbsoluteValue(col1[2][1]) select 1 else 2;
    w1 := col1[idx][1]; w2 := col2[idx][1];
    for it in [1..300] do
        if AbsoluteValue(w1) lt AbsoluteValue(w2) then tt := w1; w1 := w2; w2 := tt; end if;
        nn := Round(Real(w1/w2));
        w1n := w1 - nn*w2;
        if nn eq 0 and AbsoluteValue(w1n) ge AbsoluteValue(w2) then break; end if;
        w1 := w1n;
    end for;
    if Imaginary(w1/w2) lt 0 then w1 := -w1; end if;
    // g2, g3 and their field
    g2c, g3c, _ := EisWP(w1, w2, [ CC | ]);
    ok2, mp2 := RecAlg2(g2c, prec);
    ok3, mp3 := RecAlg2(g3c, prec);
    if not (ok2 and ok3) then printf "g2/g3 not deg<=2 algebraic, next idempotent\n"; continue; end if;
    if Degree(mp2) eq 1 and Degree(mp3) eq 1 then
        K := RationalsAsNumberField();
        emb := InfinitePlaces(K)[1]; cnj := false;
    else
        useg2 := Degree(mp2) ge Degree(mp3);
        mpK := useg2 select mp2 else mp3;
        K := NumberField(mpK/LeadingCoefficient(mpK));
        tst := useg2 select g2c else g3c;
        emb := 0; cnj := false;
        for pl in InfinitePlaces(K) do
            ev := CC!Evaluate(K.1, pl : Precision := prec);
            if AbsoluteValue(ev - tst) lt 10^(-prec div 2) then
                emb := pl; cnj := false; break;
            elif AbsoluteValue(ComplexConjugate(ev) - tst) lt 10^(-prec div 2) then
                emb := pl; cnj := true; break;
            end if;
        end for;
        if emb cmpeq 0 then printf "no matching embedding\n"; continue; end if;
    end if;
    // recognize an element of K from its complex value
    function RecK(c, K, emb, cnj, prec)
        okc, mp := RecAlg2(c, prec);
        if not okc then return false, K!0; end if;
        rts := Roots(mp, K);
        for r in rts do
            ev := CC!Evaluate(r[1], emb : Precision := prec);
            if cnj then ev := ComplexConjugate(ev); end if;
            if AbsoluteValue(ev - c) lt 10^(-prec div 2) then
                return true, r[1];
            end if;
        end for;
        return false, K!0;
    end function;
    okg2, g2K := RecK(g2c, K, emb, cnj, prec);
    okg3, g3K := RecK(g3c, K, emb, cnj, prec);
    if not (okg2 and okg3) then printf "g2/g3 not in common K\n"; continue; end if;
    printf "K def by %o; g2 = %o; g3 = %o\n", DefiningPolynomial(K), g2K, g3K;
    // sample points
    P0x := CC!xs[1]; P0y := Sqrt(Evaluate(g, P0x));
    s := (eA*ToAnalyticJacobian(P0x, P0y, A))[idx][1]
       + (eA*ToAnalyticJacobian(P0x, -P0y, A))[idx][1];
    phis := [ CC | ];
    for xi in xs do
        yi := Sqrt(Evaluate(g, CC!xi));
        zi := ToAnalyticJacobian(CC!xi, yi, A);
        Append(~phis, (eA*zi)[idx][1] - s/2);
    end for;
    PK<XK> := PolynomialRing(K);
    for d1 in [0,1], d2 in [0,1] do
        delta := (d1*w1 + d2*w2)/2;
        // prescreen
        _, _, wp3 := EisWP(w1, w2, [ phis[i] + delta : i in [2,5,9] ]);
        pre := true;
        for w in wp3 do
            okp, _ := RecK(w, K, emb, cnj, prec);
            if not okp then pre := false; break; end if;
        end for;
        printf "  shift (%o,%o): prescreen %o\n", d1, d2, pre;
        if not pre then continue; end if;
        _, _, wpall := EisWP(w1, w2, [ ph + delta : ph in phis ]);
        uvals := [ K | ]; xgood := [ QQ | ];
        for i in [1..#wpall] do
            okp, uv := RecK(wpall[i], K, emb, cnj, prec);
            if okp then Append(~uvals, uv); Append(~xgood, xs[i]); end if;
        end for;
        printf "  recognized %o of %o samples\n", #uvals, #wpall;
        if #uvals lt 13 then continue; end if;
        for n in [1..((#uvals - 4) div 2)] do
            nun := 2*n + 1;
            if #uvals lt nun + 3 then break; end if;
            rows := [];
            rhs := [ K | ];
            for i in [1..nun] do
                xi := K!xgood[i]; ui := uvals[i];
                Append(~rows, [ xi^k : k in [0..n] ] cat [ -ui*xi^k : k in [0..n-1] ]);
                Append(~rhs, ui*xi^n);
            end for;
            ok, sol := IsConsistent(Transpose(Matrix(K, nun, nun, rows)),
                                    Matrix(K, 1, nun, rhs));
            if not ok then continue; end if;
            co := Eltseq(sol);
            pQ := PK![ co[k] : k in [1..n+1] ];
            qQ := XK^n + PK![ co[n+1+k] : k in [1..n] ];
            good := true;
            for i in [nun+1..#uvals] do
                if Evaluate(qQ, K!xgood[i]) eq 0 then continue; end if;
                if Evaluate(pQ, K!xgood[i]) ne uvals[i]*Evaluate(qQ, K!xgood[i]) then
                    good := false; break;
                end if;
            end for;
            if not good then continue; end if;
            // exact verification over K: (4p^3 - g2 p q^2 - g3 q^3) q = g h^2
            NN := 4*pQ^3 - g2K*pQ*qQ^2 - g3K*qQ^3;
            num := NN*qQ;
            if num eq 0 then continue; end if;
            okd, quo := IsDivisibleBy(num, PK!gQ);
            if not okd then printf "  n=%o fit, not divisible\n", n; continue; end if;
            oksq, hQ := IsSquare(quo);
            if not oksq then printf "  n=%o fit, quotient not square\n", n; continue; end if;
            printf "CERTIFIEDK|%o|%o|%o|%o|%o|%o|%o\n",
                DefiningPolynomial(K), n, g2K, g3K,
                Coefficients(pQ), Coefficients(qQ), Coefficients(hQ);
            quit;
        end for;
    end for;
end for;
print "NOCOVER";
quit;
