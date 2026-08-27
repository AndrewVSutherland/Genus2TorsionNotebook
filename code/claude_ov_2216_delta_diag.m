//////////////////////////////////////////////////////////////////////
//  claude_ov_2216_delta_diag.m
//
//  Diagnose the correct x-T formula for the seven rational 2-torsion
//  classes of the odd M_1(8,2,2) model.
//
//  "Truth" for delta(T) is obtained WITHOUT any 2-torsion formula:
//  pick any D in J(Q) with both u_D and u_{D+T} coprime to f; then
//  delta(T) = delta(D) * delta(D+T)  (delta is a homomorphism and
//  delta(D)^2 = 1).
//
//  We then print the ratio (truth / naive-formula) for each T so the
//  correction can be read off.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned height then height := 12;
elif Type(height) eq MonStgElt then height := StringToInteger(height); end if;
if not assigned nsample then nsample := 5;
elif Type(nsample) eq MonStgElt then nsample := StringToInteger(nsample); end if;
SetMemoryLimit(8*10^9);

Qq := Rationals();
P<X> := PolynomialRing(Qq);

function OddPolynomial(u, v)
    qt := -X^2
        + (u^2*v - u^2 + u*v^2 - u - v^2 - v - 2)*X
        - (u^2 + u*v + v^2 + u + v + 1);
    return ((1-u)*X + 1)*((1-v)*X + 1)*((u+v+2)*X + 1)*qt, qt;
end function;

function SqClass(a)
    n := Numerator(a); d := Denominator(a);
    m := n*d; s := Sign(m); m := AbsoluteValue(m); r := 1;
    for pe in Factorization(m) do
        if IsOdd(pe[2]) then r *:= pe[1]; end if;
    end for;
    return s*r;
end function;

function EvalK(g, KK, th)
    val := KK!0; cs := Coefficients(g); pw := KK!1;
    for i in [1..#cs] do val +:= (KK!cs[i])*pw; pw *:= th; end for;
    return val;
end function;

function DeltaCoprime(uD, roots, KK, th)
    d := [];
    for r in roots do
        val := Evaluate(uD, r);
        if val eq 0 then return false, _; end if;
        Append(~d, SqClass(val));
    end for;
    valK := EvalK(uD, KK, th);
    if valK eq 0 then return false, _; end if;
    return true, <d[1], d[2], d[3], valK>;
end function;

function DeltaMul(a, b)
    return <SqClass(a[1]*b[1]), SqClass(a[2]*b[2]), SqClass(a[3]*b[3]), a[4]*b[4]>;
end function;

// naive formula
function DeltaNaive(uS, f, roots, KK, th)
    cofac := f div uS;
    if uS*cofac ne f then return false, _; end if;
    d := [];
    for r in roots do
        if Evaluate(uS, r) eq 0 then val := Evaluate(cofac, r);
        else val := Evaluate(uS, r); end if;
        Append(~d, SqClass(val));
    end for;
    vS := EvalK(uS, KK, th);
    if vS eq 0 then valK := EvalK(cofac, KK, th); else valK := vS; end if;
    return true, <d[1], d[2], d[3], valK>;
end function;

hs := [];
for a in [-height..height] do
    for b in [1..height] do
        if GCD(AbsoluteValue(a), b) eq 1 then Append(~hs, Qq!a/b); end if;
    end for;
end for;
hs := Sort(Setseq(Seqset(hs)));

SetSeed(20260725);
tested := 0;
for iter in [1..200000] do
    if tested ge nsample then break; end if;
    u := Random(hs); v := Random(hs);
    if u eq 0 or v eq 0 or u eq v or u eq 1 or v eq 1 then continue; end if;
    if u+v+1 eq 0 or u+v+2 eq 0 then continue; end if;
    f0, qt := OddPolynomial(u, v);
    if Degree(f0) ne 5 or Discriminant(f0) eq 0 then continue; end if;
    if not IsIrreducible(qt) then continue; end if;
    y0 := u*v*(u+v+1);
    if y0 eq 0 or y0^2 ne Evaluate(f0, Qq!-1) then continue; end if;
    Lden := 1;
    for i in [0..Degree(f0)] do Lden := LCM(Lden, Denominator(Coefficient(f0,i))); end for;
    f := P!(Lden^2*f0);
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    ok := true;
    try DQ := J![X+1, Qq!(Lden*y0)]; catch e ok := false; end try;
    if not ok then continue; end if;
    if Order(DQ) ne 8 then continue; end if;
    roots := [ Qq!rt[1] : rt in Roots(f) ];
    if #roots ne 3 then continue; end if;
    qtm := qt / LeadingCoefficient(qt);
    KK<th> := NumberField(qtm);
    c := LeadingCoefficient(f);
    tested +:= 1;

    lin := [ X - r : r in roots ];
    subsets := [ qtm ];
    labels  := [ "quadpair" ];
    for i in [1..3] do for j in [i+1..3] do
        Append(~subsets, lin[i]*lin[j]);
        Append(~labels, Sprintf("pair_%o%o", i, j));
    end for; end for;
    for i in [1..3] do
        Append(~subsets, qtm * &*[ lin[j] : j in [1..3] | j ne i ]);
        Append(~labels, Sprintf("compl_%o", i));
    end for;

    G, phi := TorsionSubgroup(J);
    pool := [ phi(g) : g in G ];

    printf "== u=%o v=%o  lc(f)=%o  sqclass(lc)=%o\n", u, v, c, SqClass(c);
    for k in [1..#subsets] do
        uS := subsets[k];
        T := J![uS, P!0];
        okn, dn := DeltaNaive(uS, f, roots, KK, th);
        // truth via a helper D
        dtrue := 0; found := false;
        for D in pool do
            if D eq J!0 then continue; end if;
            if GCD(D[1], f) ne 1 then continue; end if;
            D2 := D + T;
            if D2 eq J!0 then continue; end if;
            if GCD(D2[1], f) ne 1 then continue; end if;
            g1, d1 := DeltaCoprime(D[1], roots, KK, th);
            g2, d2 := DeltaCoprime(D2[1], roots, KK, th);
            if not (g1 and g2) then continue; end if;
            dtrue := DeltaMul(d1, d2);
            found := true; break;
        end for;
        if not found then
            printf "   %-10o  NO HELPER FOUND\n", labels[k];
            continue;
        end if;
        rat := DeltaMul(dtrue, dn);   // truth/naive  (mod squares)
        sqK, _ := IsSquare(rat[4]);
        printf "   %-10o naive=(%o,%o,%o|K) true=(%o,%o,%o|K) ratio=(%o,%o,%o|Ksq=%o) trueTrivial=%o\n",
               labels[k], dn[1],dn[2],dn[3], dtrue[1],dtrue[2],dtrue[3],
               rat[1],rat[2],rat[3], sqK, (dtrue[1] eq 1 and dtrue[2] eq 1 and dtrue[3] eq 1 and IsSquare(dtrue[4]));
    end for;
end for;
print "SEARCH_DONE";
quit;
