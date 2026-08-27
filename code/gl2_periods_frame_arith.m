// 37-hunt v15 stage 1: recover the arithmetic natively in the Periods frame.
// The Periods basis is related to Basis/Lattice(Af) by an UNKNOWN GL4
// transform (all equivariance tests fail in both assumed frames), so we
// solve for the operators directly:
//   Shat: conj(PM) = Shat * PM          (integer 4x4, involution)
//   That_l: That_l * PM = PM * diag(a_l^(1), a_l^(2))  (integer 4x4)
// via row-wise integer-relation LLL at Prec digits, with exact a_l from the
// modular symbols Hecke charpolys.  Then: Eisenstein kernel mod 37, star
// eigenvectors, star action on the NS span <E1,E2> -- all in ONE frame.
// Prints everything needed by stage 2 (the PP enumeration).
//
// Run: magma -b Prec:=100 code/gl2_periods_frame_arith.m > results/gl2_periods_frame_arith.log

SetColumns(0);
SetSeed(1);
if not assigned Prec then Prec := 100; elif Type(Prec) eq MonStgElt then Prec := StringToInteger(Prec); end if;
SetMemoryLimit(24*10^9);

SetDefaultRealField(RealField(Prec + 10));
CC := ComplexField(Prec + 10);
Q := Rationals();

Lv := 2190;
PP := 37;
trtargets := [<7,2>, <11,-6>, <13,-2>, <17,6>, <19,10>];
badtargets := [<2,1>, <3,1>, <5,1>, <73,73>];

M := ModularSymbols(Lv, 2, 0);
S := CuspidalSubspace(M);
NS_ := NewSubspace(S);
D := NewformDecomposition(NS_);
target := 0;
for i in [1..#D] do
    if Dimension(D[i]) ne 4 then continue; end if;
    ok := true;
    for tt in trtargets do
        if Trace(HeckeOperator(D[i], tt[1])) ne 2*tt[2] then ok := false; break; end if;
    end for;
    if ok then target := i; break; end if;
end for;
error if target eq 0, "piece not found";
Af := D[target];
printf "PIECE %o\n", target;

nterms := Ceiling(Lv * Prec * Log(10)/(2*Pi(RealField(20)))) + 200;
printf "PERIOD_TERMS %o\n", nterms;
Pv := Periods(Af, nterms);
PMb := Matrix(CC, 4, 2, &cat[ [CC!Pv[i][1], CC!Pv[i][2]] : i in [1..4] ]);
printf "PERIODS ok\n";

// ---- integer-relation solver: find integer row v with
//      sum_j v[j]*targ[j] = rhs  (targ: list of C-numbers), i.e.
//      v satisfies  sum v[j] targ[j] - rhs = 0.  Solve by LLL on
//      real+imag columns scaled by 10^prec.
function IntRelRow(targ, rhs, prec, bound)
    n := #targ;
    sc := 10^prec;
    rows := [];
    for j in [1..n] do
        Append(~rows, [ j eq k select 1 else 0 : k in [1..n] ]
            cat [ Round(sc*Re(targ[j])), Round(sc*Im(targ[j])) ]);
    end for;
    Append(~rows, [ 0 : k in [1..n] ] cat [ Round(sc*Re(rhs)), Round(sc*Im(rhs)) ]);
    Lm := Matrix(Integers(), n+1, n+2, &cat rows);
    Lr := LLL(Lm);
    // look for a row with last-column-block ~0 and coefficient -1 or 1 on rhs...
    // simpler: we encoded rhs as an extra vector; a relation row has small norm
    for i in [1..Nrows(Lr)] do
        r := Lr[i];
        if Abs(r[n+1]) lt sc div 10^6 and Abs(r[n+2]) lt sc div 10^6 then
            // candidate relation: v . targ + c * rhs ~ 0 -- need coefficient of rhs
            // recover: recompute c = coefficient on the rhs generator; we lost it.
            continue;
        end if;
    end for;
    return false, [0 : j in [1..n]];
end function;

// cleaner: augment identity to track ALL coefficients incl. rhs
function IntRel(vals, prec)
    // find integer relation sum c_i vals_i = 0, return c (or fail)
    n := #vals;
    sc := 10^prec;
    rows := [];
    for j in [1..n] do
        Append(~rows, [ j eq k select 1 else 0 : k in [1..n] ]
            cat [ Round(sc*Re(vals[j])), Round(sc*Im(vals[j])) ]);
    end for;
    Lm := Matrix(Integers(), n, n+2, &cat rows);
    Lr := LLL(Lm);
    r := Lr[1];
    resid := Max(Abs(r[n+1]), Abs(r[n+2]));
    coefs := [ r[j] : j in [1..n] ];
    maxc := Max([Abs(c) : c in coefs]);
    ok := resid lt sc div 10^20 and maxc gt 0 and maxc lt 10^12;
    return ok, coefs;
end function;

// ---- star matrix: conj(PM[i]) = sum_j Shat[i][j] PM[j]
Shat := ZeroMatrix(Integers(), 4, 4);
okS := true;
for i in [1..4] do
    // relation on 5 values: PM[1][k]..PM[4][k], conj(PM[i][k]) jointly for k=1,2
    vals := [ CC | PMb[j][1] : j in [1..4] ] cat [ Conjugate(PMb[i][1]) ];
    vals2 := [ CC | PMb[j][2] : j in [1..4] ] cat [ Conjugate(PMb[i][2]) ];
    // need SAME coefficients for both coordinates: stack as one relation over
    // a random C-linear combination to avoid coordinate degeneracy
    mu := CC!(31/17) + CC.1*CC!(11/23);
    valsC := [ vals[j] + mu*vals2[j] : j in [1..5] ];
    ok, cf := IntRel(valsC, Prec - 15);
    if not ok or cf[5] eq 0 then okS := false; break; end if;
    // normalize: coefficient of conj term must be -1 (up to sign)
    c5 := cf[5];
    if Abs(c5) ne 1 then okS := false; break; end if;
    for j in [1..4] do Shat[i][j] := -c5*cf[j]; end for;
end for;
error if not okS, "star LLL failed";
printf "SHAT %o\n", Eltseq(Shat);
printf "SHAT_INVOLUTION %o\n", Shat^2 eq IdentityMatrix(Integers(),4);
// verify numerically
errS := Max([ Abs(x) : x in Eltseq(ChangeRing(ChangeRing(Shat,Q),CC)*PMb
       - Matrix(CC,4,2,[Conjugate(x) : x in Eltseq(PMb)])) ]);
printf "SHAT_RESID %o\n", RealField(6)!errS;

// ---- endomorphism lattice of the torus via LLL projection ----
// Endos: integer T with T*[X|Xbar] = [X|Xbar]*diag(A,Abar), A in M2(C).
// V := real span of the 8 generators (A = E_mn, i*E_mn) is 8-dim; the
// integer points of V = End(torus) = order <1, omega> in Q(sqrt3).
Y := HorizontalJoin(PMb, Matrix(CC,4,2,[Conjugate(x) : x in Eltseq(PMb)]));
Yi := Y^-1;
Gs := [];
for m in [1..2] do for n in [1..2] do
    for im in [0,1] do
        Amn := ZeroMatrix(CC,2,2);
        Amn[m][n] := im eq 0 select CC!1 else CC.1;
        Abar := Matrix(CC,2,2,[Conjugate(x) : x in Eltseq(Amn)]);
        Dg := DiagonalJoin(Amn, Abar);
        Tm := Y * Dg * Yi;
        imax := Max([Abs(Im(x)) : x in Eltseq(Tm)]);
        error if imax gt 10^(-Prec div 2), "endo generator not real";
        Append(~Gs, Matrix(RealField(Prec), 4,4, [Re(x) : x in Eltseq(Tm)]));
    end for;
end for; end for;
// Gram-Schmidt orthonormal basis of V in R^16
RR := RealField(Prec);
vecs := [ Vector(RR, Eltseq(g)) : g in Gs ];
onb := [];
for v in vecs do
    w := v;
    for u in onb do w := w - InnerProduct(w,u)*u; end for;
    nw := Sqrt(InnerProduct(w,w));
    if nw gt 10^(-Prec div 2) then Append(~onb, w/nw); end if;
end for;
printf "ENDO_V_DIM %o\n", #onb;
// lattice: e_i augmented with scaled orthogonal-complement projection
sc := 10^(Prec - 20);
rows := [];
for i in [1..16] do
    ei := Vector(RR, [ RR | j eq i select 1 else 0 : j in [1..16] ]);
    pr := ei;
    for u in onb do pr := pr - InnerProduct(ei,u)*u; end for;   // e_i - proj_V
    Append(~rows, [ Integers() | j eq i select 1 else 0 : j in [1..16] ]
        cat [ Round(sc*pr[j]) : j in [1..16] ]);
end for;
Lm := Matrix(Integers(), 16, 32, &cat rows);
Lr := LLL(Lm);
endos := [];
for i in [1..16] do
    r := Lr[i];
    resid := Max([Abs(r[16+j]) : j in [1..16]]);
    if resid lt sc div 10^(Prec div 2) then
        Tm := Matrix(Integers(),4,4,[ r[j] : j in [1..16] ]);
        Append(~endos, Tm);
    end if;
end for;
printf "N_ENDOS_FOUND %o\n", #endos;
error if #endos lt 2, "endomorphism lattice rank < 2";
// normalize: find Ghat not a multiple of identity
Ghat := 0;
for Tm in endos do
    if Tm ne Tm[1][1]*IdentityMatrix(Integers(),4) then Ghat := Tm; break; end if;
end for;
error if Type(Ghat) eq RngIntElt, "no non-scalar endo";
printf "GHAT %o\n", Eltseq(Ghat);
cpG := CharacteristicPolynomial(ChangeRing(Ghat,Q));
printf "GHAT_CHARPOLY %o\n", cpG;
facG := Factorization(cpG);
error if #facG ne 1 or facG[1][2] ne 2 or Degree(facG[1][1]) ne 2, "Ghat charpoly not (quadratic)^2";
mpol := facG[1][1];
K3<omg> := NumberField(mpol);
error if not IsIsomorphic(K3, QuadraticField(3)), "endo field is not Q(sqrt3)!";
printf "ENDO_FIELD Q(sqrt3) confirmed; omega has minpoly %o\n", mpol;
// commutation with star: star o phi = sigma(phi) o star for RM; check either
comm := ChangeRing(Shat,Q)*ChangeRing(Ghat,Q) eq ChangeRing(Ghat,Q)*ChangeRing(Shat,Q);
printf "GHAT_SHAT_COMMUTE %o\n", comm;

// ---- Hecke matrices exactly: a_l = c + d*omega in K3 ----
function HeckeCD(l)
    Tl := HeckeOperator(Af, l);
    cp := CharacteristicPolynomial(Tl);
    fac := Factorization(ChangeRing(cp, K3));
    rts := [];
    for t in fac do
        if Degree(t[1]) eq 1 then Append(~rts, -Coefficient(t[1],0)/Coefficient(t[1],1)); end if;
    end for;
    rts := Setseq(Seqset(rts));
    error if #rts eq 0, Sprintf("no roots in K3 for l=%o", l);
    out := [];
    for r in rts do
        cf := Eltseq(r);
        Append(~out, [cf[1], cf[2]]);   // r = cf1 + cf2*omega
    end for;
    return out;
end function;

F := GF(PP);
V4 := VectorSpace(F, 4);
GhF := ChangeRing(Ghat, F);
branches := [ V4 ];
ells := [ l : l in [7,11,13,17,19,23,29,31,41,43,47,53] | Lv mod l ne 0 ];
alltargets := [ <l, 1+l> : l in ells ] cat badtargets;
for tgt in alltargets do
    l := tgt[1]; av := tgt[2];
    cds := HeckeCD(l);
    newbr := [];
    for Wb in branches do
        for cd in cds do
            c := cd[1]; d := cd[2];
            if Denominator(c) mod PP eq 0 or Denominator(d) mod PP eq 0 then continue; end if;
            MlF := (F!Numerator(c)/F!Denominator(c))*IdentityMatrix(F,4)
                 + (F!Numerator(d)/F!Denominator(d))*GhF;
            Wn := Wb meet Kernel(MlF - ScalarMatrix(F, 4, F!av));
            if Dimension(Wn) ge 1 then Append(~newbr, Wn); end if;
        end for;
    end for;
    // dedupe branches
    seenb := {};
    ub := [];
    for Wb in newbr do
        kb := { Eltseq(b) : b in Basis(Wb) };
        if kb in seenb then continue; end if;
        Include(~seenb, kb);
        Append(~ub, Wb);
    end for;
    branches := ub;
    printf "ELL %o branches %o dims %o\n", l, #branches, [Dimension(w) : w in branches];
    if #branches eq 0 then break; end if;
end for;
error if #branches eq 0, "Eisenstein kernel empty in all branches";
W := branches[1];
printf "EIS_KERNEL_DIM %o\n", Dimension(W);

// star eigenvectors in W
SI37 := ChangeRing(Shat, F);
Bw := Basis(W);
wplus := 0; wminus := 0;
for co in CartesianPower([0..PP-1], #Bw) do
    wv := &+[ co[j]*Bw[j] : j in [1..#Bw] ];
    if wv eq V4!0 then continue; end if;
    im := wv*SI37;
    if im eq wv and Type(wplus) eq RngIntElt then wplus := wv; end if;
    if im eq -wv and Type(wminus) eq RngIntElt then wminus := wv; end if;
end for;
printf "STAR_EIGVECS plus=%o minus=%o\n", wplus, wminus;

// star action on the NS span (E-side)
prs := [<1,2>,<1,3>,<1,4>,<2,3>,<2,4>,<3,4>];
f1c := [ 0, 1, -1, -1, 2, 0 ];
f2c := [ -3, 4, -2, 2, 0, 3 ];
function FMat(f)
    FM := ZeroMatrix(Q, 4, 4);
    for k in [1..6] do
        FM[prs[k][1]][prs[k][2]] := f[k];
        FM[prs[k][2]][prs[k][1]] := -f[k];
    end for;
    return FM;
end function;
F1 := FMat(f1c); F2 := FMat(f2c);
function Pf(FM) return FM[1][2]*FM[3][4] - FM[1][3]*FM[2][4] + FM[1][4]*FM[2][3]; end function;
E1 := Pf(F1)*F1^-1; E2 := Pf(F2)*F2^-1;
ShQ := ChangeRing(Shat, Q);
SE1 := ShQ*ChangeRing(E1,Q)*Transpose(ShQ);
SE2 := ShQ*ChangeRing(E2,Q)*Transpose(ShQ);
v1 := Vector(Q, Eltseq(ChangeRing(E1,Q))); v2 := Vector(Q, Eltseq(ChangeRing(E2,Q)));
Mv := Matrix(Q, 2, 16, Eltseq(v1) cat Eltseq(v2));
r1 := Solution(Mv, Vector(Q, Eltseq(SE1)));
r2 := Solution(Mv, Vector(Q, Eltseq(SE2)));
A2 := Matrix(Q, 2, 2, [r1[1], r1[2], r2[1], r2[2]]);
printf "STAR_NS_ACTION %o\n", Eltseq(A2);
Km := KernelMatrix(A2 + IdentityMatrix(Q,2));
Kp := KernelMatrix(A2 - IdentityMatrix(Q,2));
printf "MINUS_EIG dim %o gen %o\n", Nrows(Km), Nrows(Km) gt 0 select Eltseq(Km[1]) else [];
printf "PLUS_EIG dim %o gen %o\n", Nrows(Kp), Nrows(Kp) gt 0 select Eltseq(Kp[1]) else [];
printf "GL2_FRAME_ARITH_DONE\n";
quit;
