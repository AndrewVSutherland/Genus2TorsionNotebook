// BLP [2,22] engine, part B: scan of the 2-rank-2 slice.
//
// Slice construction (this session; see notes/nonrm_222_hunt_2026_07_31.md
// route 5): f = (R-2cS)(R+2cS) =: g*h with h = g + 4cS, S = x^2+d.
// Prescribe g = (x-r)(x^2+px+q) (rational root r) and h(s) = 0 (rational
// root s != r).  The two normalizations
//   (i)  R = g + 2cS has x^2-coefficient -1  =>  c = (r-1-p)/2,
//   (ii) h(s) = 0  =>  2(r-1-p)(s^2+d) + (s-r)(s^2+ps+q) = 0,
// give p LINEARLY:  p = [2(r-1)(s^2+d) + (s-r)(s^2+q)] / [2(s^2+d) - s(s-r)].
// Then f has factor type [1,1,2,2] (or finer) => 2-rank 2 automatically,
// and the order-11 closure t0 = t1 = 0 (DeriveT, validated in
// code/blp22_engine_derive.m) cuts the [2,22] locus.  For fixed (r,s,d)
// the conditions are univariate in q: common rational roots of
// gcd(Numer t0, Numer t1) give complete [2,22] CANDIDATE curves.
// Each candidate: CFOrder must return 11 (independent D_inf check), then
// exact TorsionSubgroup (expect >= [2,22]), printed as HIT22.
//
// Run: magma -b Hr:=4 NParts:=1 Part:=0 code/blp22_engine_scan.m

SetColumns(0);
SetSeed(1);
if not assigned Hr then Hr := 4; elif Type(Hr) eq MonStgElt then Hr := StringToInteger(Hr); end if;
if not assigned NParts then NParts := 1; elif Type(NParts) eq MonStgElt then NParts := StringToInteger(NParts); end if;
if not assigned Part then Part := 0; elif Type(Part) eq MonStgElt then Part := StringToInteger(Part); end if;
if not assigned MemGB then MemGB := 3; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();
PQ<X> := PolynomialRing(Q);

// ---- DeriveT (verbatim from code/blp22_engine_derive.m; keep in sync) ----
function DeriveT(a, b, c, d)
    K := Parent(a);
    P := PolynomialRing(K); x := P.1;
    R := x^3 - x^2 + a*x + b;
    S := x^2 + d;
    F := R^2 - 4*c^2*S^2;
    if Degree(F) ne 6 then return K!0, K!0, false; end if;
    fc := [Coefficient(F, i) : i in [0..6]];
    s := [K| 1];
    for k in [1..9] do
        Ftk := k le 6 select fc[6-k+1] else K!0;
        acc := Ftk;
        for i in [1..k-1] do
            acc -:= s[i+1]*s[k-i+1];
        end for;
        Append(~s, acc/2);
    end for;
    M := Matrix(K, 3, 3, [ [ s[3+m+1], s[4+m+1], s[5+m+1] ] : m in [1..3] ]);
    rhs := Vector(K, [ -s[6+m+1] : m in [1..3] ]);
    if Determinant(M) eq 0 then return K!0, K!0, false; end if;
    u := Solution(Transpose(M), rhs);
    U := x^3 + u[3]*x^2 + u[2]*x + u[1];
    uco := [u[1], u[2], u[3], K!1];
    V := P!0;
    for n in [0..6] do
        vn := K!0;
        for i in [Max(0, n-3)..3] do
            j := i + 3 - n;
            if j ge 0 and j le 9 then vn +:= uco[i+1]*s[j+1]; end if;
        end for;
        V +:= vn*x^n;
    end for;
    G := F*U^2 - V^2;
    if Degree(G) gt 5 then return K!0, K!0, false; end if;
    t := Coefficient(G, 2);
    return Coefficient(G, 0) - t*d, Coefficient(G, 1), true;
end function;

// ---- CFOrder (pell-cf-order; independent exact D_inf order check) ----
function SqrtPolyPart6(f)
    P := Parent(f); xx := P.1;
    sp := xx^3;
    for k in [1..3] do
        dd := f - sp^2;
        if Degree(dd) le 2 then break; end if;
        sp := sp + (Coefficient(dd, 6-k)/(2*Coefficient(sp, 3)))*xx^(3-k);
    end for;
    return sp;
end function;
function CFOrderB(f, maxsteps, maxord)
    P := Parent(f);
    sp := SqrtPolyPart6(f);
    Pi := P!0; Qi := P!1; total := 0;
    for i in [0..maxsteps] do
        if Qi eq 0 then return 0; end if;
        ai := (Pi + sp) div Qi;
        total +:= Degree(ai);
        if total gt maxord then return 0; end if;
        Pn := ai*Qi - Pi;
        if (f - Pn^2) mod Qi ne 0 then return 0; end if;
        Qn := (f - Pn^2) div Qi;
        Pi := Pn; Qi := Qn;
        if i ge 1 and Degree(Qi) le 0 and Qi ne 0 then return total; end if;
    end for;
    return 0;
end function;

function ExactTorsionInvs(f)
    C := HyperellipticCurve(f);
    try
        Cm := ReducedMinimalWeierstrassModel(C);
        C := SimplifiedModel(Cm);
    catch e ; end try;
    return Invariants(TorsionSubgroup(Jacobian(C)));
end function;

// rational values of height <= Hr (numerator and denominator), with 0
vals := [Q!0];
for nn in [1..Hr] do for mm in [1..Hr] do
    if GCD(nn, mm) eq 1 then vals cat:= [Q!nn/mm, -Q!nn/mm]; end if;
end for; end for;
printf "PARAM_VALUES %o (Hr=%o)\n", #vals, Hr;

Kq<qv> := RationalFunctionField(Q);
Pq := PolynomialRing(Q);

tested := 0; syst := 0; hits := 0; idx := -1;
for r in vals do
 for sv in vals do
  if sv eq r then continue; end if;
  for d in vals do
    idx +:= 1;
    if idx mod NParts ne Part then continue; end if;
    if sv^2 + d eq 0 then continue; end if;
    den := 2*(sv^2 + d) - sv*(sv - r);
    if den eq 0 then continue; end if;
    tested +:= 1;
    if tested mod 500 eq 0 then
        printf "PROGRESS tested=%o syst=%o hits=%o\n", tested, syst, hits;
    end if;
    // p(q), c(q), a(q), b(q) over Kq
    p := (2*(r-1)*(sv^2+d) + (sv-r)*(sv^2+qv)) / den;
    c := (r - 1 - p)/2;
    // g = (x-r)(x^2+px+q) = x^3 + (p-r)x^2 + (q-rp)x - rq;
    // R = g + 2c(x^2+d): x^2 coeff p-r+2c = -1 by construction, so
    // a = q - rp and b = -rq + 2cd.
    aq := qv - r*p;
    bq := -r*qv + 2*c*d;
    t0, t1, ok := DeriveT(aq, bq, c, d);
    if not ok then continue; end if;
    if t0 eq 0 and t1 eq 0 then
        printf "IDENTIC r=%o s=%o d=%o (t0 = t1 = 0 identically in q!)\n", r, sv, d;
        continue;
    end if;
    syst +:= 1;
    N0 := Pq!Numerator(t0); N1 := Pq!Numerator(t1);
    if N0 eq 0 or N1 eq 0 then
        gpol := N0 eq 0 select N1 else N0;
    else
        gpol := GCD(N0, N1);
    end if;
    if Degree(gpol) lt 1 then continue; end if;
    for rt in Roots(gpol) do
        q0 := rt[1];
        // rebuild exactly over Q
        p0 := Evaluate(p, q0); c0 := Evaluate(c, q0);
        if c0 eq 0 then continue; end if;
        g0 := (X - r)*(X^2 + p0*X + q0);
        h0 := g0 + 4*c0*(X^2 + d);
        f := g0*h0;
        if Degree(f) ne 6 or Discriminant(f) eq 0 then continue; end if;
        assert Evaluate(g0, r) eq 0 and Evaluate(h0, sv) eq 0;
        // independent D_inf order check via CF (f monic sextic)
        cford := CFOrderB(f, 60, 45);
        dtype := {* Degree(fp[1])^^fp[2] : fp in Factorization(f) *};
        printf "CAND r=%o s=%o d=%o q=%o cford=%o type=%o\n",
               r, sv, d, q0, cford, dtype;
        if cford ne 0 and cford mod 11 ne 0 then continue; end if;
        invs := ExactTorsionInvs(f);
        printf "EXACT r=%o s=%o d=%o q=%o TORSION=%o\n", r, sv, d, q0, invs;
        if #invs ge 2 and invs[#invs] mod 11 eq 0 and invs[#invs-1] mod 2 eq 0 then
            hits +:= 1;
            printf "HIT22 r=%o s=%o d=%o q=%o TORSION=%o f=%o\n",
                   r, sv, d, q0, invs, f;
        end if;
    end for;
  end for;
 end for;
end for;
printf "SEARCH_DONE Hr=%o Part=%o/%o tested=%o syst=%o hits=%o\n",
       Hr, Part, NParts, tested, syst, hits;
quit;
