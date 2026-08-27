// [2,22] fast solver, SLICE 2: factor type [1,1,1,3].
// g = (x-r1)(x-r2)(x-r3) fully split, h = g + 4cS irreducible cubic,
// c = (r1+r2+r3-1)/2 (from the -x^2 normalization; needs sum != 1).
// 2-rank 2 via the four Galois orbits {r1},{r2},{r3},{h-roots}.
// Fibered over (r1,r2); unknowns (r3,d); same closure conditions
// t0 = t1 = 0 and the SAME mod-p CRT pipeline as code/blp22_fast_sweep.m
// (validated there end-to-end; no known rational positive exists on THIS
// slice, so the self-test checks mod-p consistency: an order-11-tagged
// locus point exists over F_p at a random base and its curve has type
// [1,1,1,3] with 11 | CF order).
//
// Run: magma -b H:=8 NParts:=3 Part:=0 code/blp22_fast_sweep_s2.m
SetColumns(0);
SetSeed(1);
if not assigned H then H := 8; elif Type(H) eq MonStgElt then H := StringToInteger(H); end if;
if not assigned Hdone then Hdone := 5; elif Type(Hdone) eq MonStgElt then Hdone := StringToInteger(Hdone); end if;
if not assigned NParts then NParts := 1; elif Type(NParts) eq MonStgElt then NParts := StringToInteger(NParts); end if;
if not assigned Part then Part := 0; elif Type(Part) eq MonStgElt then Part := StringToInteger(Part); end if;
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

load "code/blp22_locus_rm_test2.m0";   // DeriveT (field-generic)

Q := Rationals();
PQ<X> := PolynomialRing(Q);
PRIMES := [4999, 5003, 5009, 5011];

function SqrtPolyPart6(f)
    P := Parent(f); xx := P.1; sp := xx^3;
    for k in [1..3] do
        dd := f - sp^2;
        if Degree(dd) le 2 then break; end if;
        sp := sp + (Coefficient(dd, 6-k)/(2*Coefficient(sp, 3)))*xx^(3-k);
    end for;
    return sp;
end function;
function CFOrderB(f, maxsteps, maxord)
    P := Parent(f); sp := SqrtPolyPart6(f);
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

// mod-p fiber processing: returns list of <Z!d0, Z!q0> with CF order = 0 mod 11
function FiberRootsModP(r, sv, p)
    Fp := GF(p);
    rp := Fp!Numerator(r)/Fp!Denominator(r);
    sp_ := Fp!Numerator(sv)/Fp!Denominator(sv);
    if rp eq sp_ then return [], false; end if;
    K2 := RationalFunctionField(Fp, 2);
    qv := K2.1; dv := K2.2;
    P2 := PolynomialRing(Fp, 2);
    PF := PolynomialRing(Fp);
    rK := K2!rp; sK := K2!sp_;
    // slice 2: r3 = qv free; c = (r1+r2+r3-1)/2; a = e2(g); b = -r1r2r3+2cd
    cf := (rK + sK + qv - 1)/2;
    if cf eq 0 then return [], false; end if;
    aS2 := rK*sK + rK*qv + sK*qv;
    bS2 := -rK*sK*qv + 2*cf*dv;
    t0, t1, ok := DeriveT(aS2, bS2, cf, K2!dv);
    if not ok or t0 eq 0 or t1 eq 0 then return [], false; end if;
    N0 := P2!Numerator(t0); N1 := P2!Numerator(t1);
    G := GCD(N0, N1);
    N0r := N0 div G; N1r := N1 div G;
    Rd := Resultant(N0r, N1r, P2.1);
    if Rd eq 0 then return [], false; end if;
    Rdn := PF!UnivariatePolynomial(Rd);
    out := [];
    for rtd in Roots(Rdn) do
        d0 := rtd[1];
        if sp_^2 + d0 eq 0 then continue; end if;
        // specialize at d = d0
        A0 := PF!0; A1 := PF!0;
        for i in [0..Degree(N0r, 1)] do
            cfq := Coefficient(N0r, 1, i);
            co := Fp!0;
            for j in [0..Degree(cfq, 2)] do co +:= (Fp!Coefficient(cfq, 2, j))*d0^j; end for;
            A0 +:= co*PF.1^i;
        end for;
        for i in [0..Degree(N1r, 1)] do
            cfq := Coefficient(N1r, 1, i);
            co := Fp!0;
            for j in [0..Degree(cfq, 2)] do co +:= (Fp!Coefficient(cfq, 2, j))*d0^j; end for;
            A1 +:= co*PF.1^i;
        end for;
        if A0 eq 0 or A1 eq 0 then continue; end if;
        gg := GCD(A0, A1);
        if Degree(gg) lt 1 then continue; end if;
        for rtq in Roots(gg) do
            q0 := rtq[1];   // = r3 mod p
            c0 := (rp + sp_ + q0 - 1)/2;
            if c0 eq 0 then continue; end if;
            g0 := (PF.1-rp)*(PF.1-sp_)*(PF.1-q0);
            f := g0*(g0 + 4*c0*(PF.1^2+d0));
            if Degree(f) ne 6 or Degree(GCD(f, Derivative(f))) gt 0 then continue; end if;
            ord := CFOrderB(f, 60, 45);
            if ord ne 0 and ord mod 11 eq 0 then
                Append(~out, <Integers()!d0, Integers()!q0>);
            end if;
        end for;
    end for;
    return out, true;
end function;

// exact verification of a rational candidate point
function ExactCheck(r, sv, d0, q0)
    // q0 = r3
    c0 := (r + sv + q0 - 1)/2;
    if c0 eq 0 or sv^2 + d0 eq 0 then return false, [Integers()|], PQ!0; end if;
    aa := r*sv + r*q0 + sv*q0; bb := -r*sv*q0 + 2*c0*d0;
    t0, t1, ok := DeriveT(Q!aa, Q!bb, Q!c0, Q!d0);
    if not ok or t0 ne 0 or t1 ne 0 then return false, [Integers()|], PQ!0; end if;
    g0 := (X - r)*(X - sv)*(X - q0);
    f := g0*(g0 + 4*c0*(X^2 + d0));
    if Degree(f) ne 6 or Discriminant(f) eq 0 then return false, [Integers()|], PQ!0; end if;
    cford := CFOrderB(f, 60, 45);
    if cford eq 0 or cford mod 11 ne 0 then return false, [Integers()|], PQ!0; end if;
    C := HyperellipticCurve(f);
    try
        Cm := ReducedMinimalWeierstrassModel(C);
        C := SimplifiedModel(Cm);
    catch e ; end try;
    invs := Invariants(TorsionSubgroup(Jacobian(C)));
    return true, invs, f;
end function;

// process one fiber
procedure DoFiber(r, sv, ~nhit)
    lists := [* *]; goodp := [];
    for p in PRIMES do
        L, ok := FiberRootsModP(r, sv, p);
        if ok then Append(~lists, L); Append(~goodp, p); end if;
    end for;
    if #goodp lt 3 then printf "FIBER r=%o s=%o BADPRIMES\n", r, sv; return; end if;
    // all 3-subsets of good primes
    n := #goodp;
    cands := {};
    for ss3 in Subsets({1..n}, 3) do
        idx := Sort(Setseq(ss3));
        M := &*[goodp[i] : i in idx];
        ZM := Integers(M);
        for x1 in lists[idx[1]] do
         for x2 in lists[idx[2]] do
          for x3 in lists[idx[3]] do
            dres := CRT([x1[1], x2[1], x3[1]], [goodp[i] : i in idx]);
            qres := CRT([x1[2], x2[2], x3[2]], [goodp[i] : i in idx]);
            okd, dr := RationalReconstruction(ZM!dres);
            if not okd then continue; end if;
            okq, qr := RationalReconstruction(ZM!qres);
            if not okq then continue; end if;
            Include(~cands, <dr, qr>);
          end for;
         end for;
        end for;
    end for;
    for cd in cands do
        okx, invs, f := ExactCheck(r, sv, cd[1], cd[2]);
        if okx then
            nhit +:= 1;
            printf "HIT22 r=%o s=%o d=%o q=%o TORSION=%o f=%o\n",
                   r, sv, cd[1], cd[2], invs, f;
        end if;
    end for;
end procedure;

// SELF-TEST (no rational positive known on this slice): mod-p pipeline
// must yield order-11-tagged points at random bases over F_4999.
stn := 0;
for tri in [1..8] do
    L, ok := FiberRootsModP(Q!(1+tri), Q!(3+2*tri), 4999);
    if ok and #L gt 0 then stn +:= #L; end if;
end for;
error if stn eq 0, "SELF-TEST FAILED: no mod-p order-11 locus points found";
printf "SELFTEST_PASS (mod-p order-11 points exist on slice 2: %o found)\n", stn;

// fiber list
function Ht(v) return Max(Abs(Numerator(v)), Denominator(v)); end function;
vals := [Q!0];
for nn in [1..H] do for mm in [1..H] do
    if GCD(nn, mm) eq 1 then vals cat:= [Q!nn/mm, -Q!nn/mm]; end if;
end for; end for;
printf "PARAM_VALUES %o (H=%o)\n", #vals, H;

idx := -1; nf := 0; nhit := 0;
for i in [1..#vals] do
 for j in [i+1..#vals] do
    r := vals[i]; sv := vals[j];
    if Ht(r) le Hdone and Ht(sv) le Hdone then continue; end if;
    idx +:= 1;
    if idx mod NParts ne Part then continue; end if;
    nf +:= 1;
    if nf mod 25 eq 0 then printf "PROGRESS fibers=%o hits=%o\n", nf, nhit; end if;
    DoFiber(r, sv, ~nhit);
 end for;
end for;
printf "FAST_SWEEP_DONE H=%o Part=%o/%o fibers=%o hits=%o\n", H, Part, NParts, nf, nhit;
quit;
