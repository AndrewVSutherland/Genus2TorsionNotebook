// lane_c33.m — (3,3)-gluing hunt for the full-product gap groups via the
// CONTRAVARIANT PENCIL (unparks the 2026-08-12 z3 lane, whose blind Z(3,2)
// surface sweep was torsion-thin by construction).
//
// Engine: for E1 = E1(t) in a genus-0 torsion family, let U be a degree-3
// genus-one model of E1 and (P,Q) := Contravariants(U) (Fisher).  The pencil
//     E2(lam) := Jacobian(P + lam*Q),   lam in P^1(Q),
// is the family of ANTI-3-congruent partners of E1 with the correct twists
// built in (runtime-verified in the 2026-08-26 session: members are
// 3-congruent to E1 and Genus2Elliptic3 glues them; the covariant Hesse
// pencil U + lam*H is the symplectic family and never glues).  We sweep lam,
// prefilter E2(lam) for the target torsion T2 by exact divisibility of
// #E2(F_p) (no twist ambiguity), check exact elliptic torsion on survivors,
// glue with Genus2Elliptic3 (BHLS, genus2.m) and funnel (split_lab, early
// abort vs KNOWN before any exact genus-2 torsion computation).
//
// For a (3,3)-gluing with gcd(#T1 * #T2, 3) = 1 the product T1 x T2 injects
// into J(Q), so any target-pair hit realizes the full product:
//   [10]x[10] -> [10,10]     [8]x[10] -> [2,40]    [7]x[8]  -> [56]
//   [10]x[2,4] -> [2,2,20]   [7]x[2,4] -> [2,28]   [8]x[2,8] -> [2,8,8]
//   [10]x[2,8] -> [2,2,40]   [7]x[2,8] -> [2,56]   [2,8]x[2,8] -> [2,2,8,8]
//
// Usage (from product/code/):
//   magma -b N1:=10 T2:=10 TH:=20 SH:=40 Part:=1 NParts:=1 lane_c33.m > ../logs/c33_10x10_p1.log
// N1 in {7,8,10,24,28} (24 = the [2,4] family, 28 = HLP E_{2,8});
// T2 in {7,8,10,24,28} (torsion the partner must contain).
SetColumns(0);
if not assigned MemGB then MemGB := 6; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
if not assigned N1 then N1 := 10; elif Type(N1) eq MonStgElt then N1 := StringToInteger(N1); end if;
if not assigned T2 then T2 := 10; elif Type(T2) eq MonStgElt then T2 := StringToInteger(T2); end if;
if not assigned TH then TH := 20; elif Type(TH) eq MonStgElt then TH := StringToInteger(TH); end if;
if not assigned SH then SH := 40; elif Type(SH) eq MonStgElt then SH := StringToInteger(SH); end if;
if not assigned Part then Part := 1; elif Type(Part) eq MonStgElt then Part := StringToInteger(Part); end if;
if not assigned NParts then NParts := 1; elif Type(NParts) eq MonStgElt then NParts := StringToInteger(NParts); end if;

load "split_lab.m";  // attaches Drew's spec + ../../genus2.m, defines Funnel/KNOWN

Q := Rationals();
Px<x> := PolynomialRing(Q);

// ---------- torsion families (runtime-verified below) ----------
function FamilyCurve(N, tv)
    if N eq 7 then
        bv := tv^3-tv^2; cv := tv^2-tv;
    elif N eq 8 then
        if tv eq 0 then return false, 0; end if;
        bv := (2*tv-1)*(tv-1); cv := (2*tv-1)*(tv-1)/tv;
    elif N eq 10 then
        de := tv^2-3*tv+1;
        if de eq 0 then return false, 0; end if;
        bv := tv^3*(tv-1)*(2*tv-1)/de^2; cv := -tv*(tv-1)*(2*tv-1)/de;
    elif N eq 24 then
        E := 0;
        try E := EllipticCurve([0, 1+tv^2, 0, tv^2, 0]); catch e; return false, 0; end try;
        return true, E;
    elif N eq 28 then
        u := tv;
        if u eq 0 or 8*u^2-1 eq 0 or 4*u+1 eq 0 then return false, 0; end if;
        e1 := (16*u^3+12*u^2+2*u)/(8*u^2-1)^2;
        e2 := (32*u^3+24*u^2+8*u+1)/(16*u^2*(8*u^2-1));
        e3 := (-32*u^4-32*u^3-12*u^2-2*u)/((4*u+1)^2*(8*u^2-1));
        E := 0;
        try E := EllipticCurve((x-e1)*(x-e2)*(x-e3)); catch e; return false, 0; end try;
        return true, E;
    else
        error "unknown family", N;
    end if;
    if bv eq 0 then return false, 0; end if;
    E := 0;
    try E := EllipticCurve([1-cv, -bv, -bv, 0, 0]); catch e; return false, 0; end try;
    return true, E;
end function;

function FamTorsOK(N, I)
    if N eq 24 then return #I eq 2 and I[1] mod 2 eq 0 and I[2] mod 4 eq 0;
    elif N eq 28 then return I eq [2,8];
    else return I eq [N];
    end if;
end function;

nv := 0;
for tv in [Q| 3, 5/2, -7/3 ] do
    ok, E := FamilyCurve(N1, tv);
    if ok and FamTorsOK(N1, Invariants(TorsionSubgroup(E))) then nv +:= 1; end if;
end for;
error if nv lt 2, "family", N1, "failed runtime verification";
printf "FAMILY N1=%o verified %o/3\n", N1, nv;

// target spec: ordm = full order of T2 (divisibility gate), tors test on invs
if T2 eq 7 then ordm := 7; elif T2 eq 8 then ordm := 8; elif T2 eq 10 then ordm := 10;
elif T2 eq 24 then ordm := 8; elif T2 eq 28 then ordm := 16;
else error "unknown T2", T2; end if;

function T2OK(T2, I)
    if T2 in {7,8,10} then return I eq [T2];
    elif T2 eq 24 then return #I eq 2 and I[1] mod 2 eq 0 and I[2] mod 4 eq 0;
    elif T2 eq 28 then return #I eq 2 and I[1] mod 2 eq 0 and I[2] mod 8 eq 0;
    end if;
    return false;
end function;

PRE_P := [211, 223, 227, 229, 233, 239];

// exact divisibility gate on the explicit curve E2 (no twist ambiguity)
function ModpGate(E2, ordm)
    a := aInvariants(E2);
    ntest := 0;
    for p in PRE_P do
        Fp := GF(p);
        bad := false;
        for c in a do
            if Denominator(c) mod p eq 0 then bad := true; break; end if;
        end for;
        if bad then continue; end if;
        Ep := 0; okp := true;
        try Ep := EllipticCurve([Fp!Numerator(c)/Fp!Denominator(c) : c in a]); catch e; okp := false; end try;
        if not okp then continue; end if;
        ntest +:= 1;
        if #Ep mod ordm ne 0 then return false; end if;
    end for;
    return ntest ge 4;
end function;

function HeightRats(H)
    S := [];
    for a in [1..H], b in [1..H] do
        if GCD(a,b) eq 1 then Append(~S, Q!a/b); Append(~S, Q!-a/b); end if;
    end for;
    return S;
end function;

tvals := HeightRats(TH);
lamvals := [Q|0] cat HeightRats(SH);
nt := 0; nlam := 0; npre := 0; ncand := 0; nglue := 0; nhit := 0;
t0c := Cputime();

for ti in [1..#tvals] do
    if ti mod NParts ne (Part - 1) then continue; end if;
    tv := tvals[ti];
    ok, E1 := FamilyCurve(N1, tv);
    if not ok then continue; end if;
    if not FamTorsOK(N1, Invariants(TorsionSubgroup(E1))) then continue; end if;
    nt +:= 1;
    U := 0; Pc := 0; Qc := 0; okU := true;
    try
        U := GenusOneModel(3, E1);
        Pc, Qc := Contravariants(U);
    catch e;
        okU := false;
    end try;
    if not okU then printf "SKIPT t=%o (model)\n", tv; continue; end if;
    j1 := jInvariant(E1);
    seenj := {Q| j1};
    for li in [0..#lamvals] do
        E2 := 0; okm := true;
        try
            M := li eq 0 select Qc else Pc + lamvals[li]*Qc;
            E2 := Jacobian(M);
        catch e;
            okm := false;
        end try;
        if not okm then continue; end if;
        nlam +:= 1;
        j2 := jInvariant(E2);
        if j2 in seenj then continue; end if;
        Include(~seenj, j2);
        if not ModpGate(E2, ordm) then continue; end if;
        npre +:= 1;
        I2 := Invariants(TorsionSubgroup(E2));
        lamstr := li eq 0 select "oo" else Sprintf("%o", lamvals[li]);
        printf "PRECAND t=%o lam=%o tors2=%o\n", tv, lamstr, I2;
        if not T2OK(T2, I2) then continue; end if;
        ncand +:= 1;
        glues := [];
        try glues := Genus2Elliptic3(E1, E2); catch e; printf "GLUEFAIL t=%o lam=%o: %o\n", tv, lamstr, e`Object; end try;
        for C in glues do
            nglue +:= 1;
            res := Funnel(C, Sprintf("c33|%o(t=%o)x%o(lam=%o)", N1, tv, T2, lamstr));
            if res eq "hit" then nhit +:= 1; end if;
        end for;
    end for;
    if nt mod 10 eq 0 then
        printf "PROGRESS t#%o lam=%o pre=%o cand=%o glue=%o hit=%o %os\n",
            nt, nlam, npre, ncand, nglue, nhit, Cputime()-t0c;
    end if;
end for;
printf "SEARCH_DONE c33 N1=%o T2=%o TH=%o SH=%o part %o/%o t=%o lam=%o pre=%o cand=%o glue=%o hits=%o %.1o s\n",
    N1, T2, TH, SH, Part, NParts, nt, nlam, npre, ncand, nglue, nhit, Cputime()-t0c;
quit;
