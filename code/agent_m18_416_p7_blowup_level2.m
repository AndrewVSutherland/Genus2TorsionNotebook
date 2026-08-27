
//////////////////////////////////////////////////////////////////////
//  p=7 blowup, LEVEL 2 (mod 49 -> mod 343) for the [4,16] second halving.
//
//  Continues code/agent_m18_416_p7_blowup.m.  Level 1 left three flat
//  strata w0=0 (<0,0>,<1,0>,<5,0>) fully surviving (49/49) because the
//  first-order expansion is degenerate there.  This script resolves them
//  by explicit Hensel lifting of the exact E416 system through mod 343.
//
//  For a fixed (R,w) mod 343 we lift each mod-7 aux solution aux0=(a,b,c,d,e)
//  one 7-adic digit at a time:  aux -> aux + 7^k * sigma with
//     J_aux(aux) * sigma = -carry,   carry = E416(R,w,aux)/7^k mod 7,
//  enumerating the affine solution coset (particular + ker J_aux).  A
//  (R,w) mod 343 SURVIVES iff some aux0 lifts all the way to mod 343.
//
//  Strata with a large aux nullspace (e.g. <1,0>, where J_aux == 0) are
//  reported as UNRESOLVED beyond the kernel-dim cap rather than guessed.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
p := 7;
KerDimCap := 3;   // enumerate aux cosets only up to this nullspace dim
if assigned kercap then KerDimCap := StringToInteger(kercap); end if;

// residues to resolve (default: the three flat w0=0 strata)
resid := [[0,0],[1,0],[5,0]];
if assigned only then
    parts := Split(only, ";");
    resid := [];
    for pr in parts do xy := Split(pr, ","); Append(~resid, [StringToInteger(xy[1]), StringToInteger(xy[2])]); end for;
end if;

// ---- build cleared integer E416 ----
Rng<R,w,a,b,c,d,e> := PolynomialRing(Q, 7, "grevlex");
KF := FieldOfFractions(Rng); PX<x> := PolynomialRing(KF);
t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
A := x^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*x + R^4;
B := (R + 2 + 4*t)*x^2 + (R^2 + 4*R + 1 + 8*t)*x + (2*R^2 + R + 4*t);
f := x*A*B; c4 := R + 2 + 4*t;
q416 := x^2 + a*x + b; ell416 := c*x^2 + d*x + e;
F416 := f - ell416^2 - c4*(x+R)*q416^2;
E416 := [];
for i in [0..4] do
    ci := Numerator(Coefficient(F416, i));
    den := LCM([Denominator(cc) : cc in Coefficients(ci)]); ci := den*ci;
    g := GCD([Integers()!cc : cc in Coefficients(ci)]); Append(~E416, Rng!(ci/g));
end for;
auxidx := [3,4,5,6,7];
DerAux := [[Derivative(E416[i], auxidx[j]) : j in [1..5]] : i in [1..5]];
Fp := GF(p); Z := Integers();

function Ev(pol, v) return Z!Evaluate(pol, v); end function;

// mod-7 aux solutions at (R0,w0)
function AuxMod7(R0,w0)
    sols := [];
    for aa in [0..6] do for bb in [0..6] do for cc in [0..6] do for dd in [0..6] do for ee in [0..6] do
        v := [R0,w0,aa,bb,cc,dd,ee]; ok:=true;
        for i in [1..5] do if Ev(E416[i],v) mod p ne 0 then ok:=false; break; end if; end for;
        if ok then Append(~sols,[aa,bb,cc,dd,ee]); end if;
    end for; end for; end for; end for; end for;
    return sols;
end function;

// given full integer (R,w) and an aux vector solving mod 7^k, return the list
// of aux vectors solving mod 7^(k+1) (as integer 5-vectors), or the string
// "CAP" if the nullspace is too large to enumerate.
function LiftStep(Rw, aux, k)
    pk := p^k;
    base := [Rw[1], Rw[2]] cat aux;
    carry := [ Fp!((Ev(E416[i], base)) div pk) : i in [1..5] ];  // E416==0 mod pk assumed
    J := Matrix(Fp,5,5,[ Fp!(Ev(DerAux[i][j], base)) : i in [1..5], j in [1..5] ]);
    rhs := Vector(Fp,[ -carry[i] : i in [1..5] ]);
    consistent, part := IsConsistent(Transpose(J), rhs);   // solves x*J^T = rhs, x=sigma
    if not consistent then return []; end if;
    ker := KernelMatrix(Transpose(J));   // rows sigma with J*sigma=0
    dk := Nrows(ker);
    if dk gt KerDimCap then return "CAP"; end if;
    out := [];
    for cnt in [0..p^dk - 1] do
        coeffs := []; x0 := cnt;
        for _ in [1..dk] do Append(~coeffs, x0 mod p); x0 := x0 div p; end for;
        sig := part;
        for r in [1..dk] do sig := sig + (Fp!coeffs[r])*ker[r]; end for;
        Append(~out, [ aux[j] + pk*(Z!sig[j]) : j in [1..5] ]);
    end for;
    return out;
end function;

// does (R,w) mod 343 admit an aux solution mod 343?  returns true/false/"CAP"
function Survives343(Rw, aux7list)
    capped := false;
    for aux0 in aux7list do
        cur := [aux0];
        okchain := true;
        for k in [1,2] do
            nxt := [];
            for aux in cur do
                r := LiftStep(Rw, aux, k);
                if Type(r) eq MonStgElt then capped := true; okchain := false; break; end if;
                nxt cat:= r;
            end for;
            if not okchain then break; end if;
            cur := nxt;
            if #cur eq 0 then break; end if;
        end for;
        if okchain and #cur gt 0 then return true; end if;
    end for;
    if capped then return "CAP"; end if;
    return false;
end function;

printf "M18_416_P7_BLOWUP_LEVEL2 p=%o kercap=%o\n", p, KerDimCap;
for rw0 in resid do
    R0 := rw0[1]; w0 := rw0[2];
    aux7 := AuxMod7(R0,w0);
    nsurv := 0; ncap := 0; total := 0;
    // (R,w) mod 343 reducing to (R0,w0) mod 7: 49 lifts each
    for ri in [0..48] do for wi in [0..48] do
        Rv := R0 + p*ri; Wv := w0 + p*wi;   // ranges over all residues mod 343 == (R0,w0) mod 7
        total +:= 1;
        s := Survives343([Rv,Wv], aux7);
        if Type(s) eq MonStgElt then ncap +:= 1;
        elif s then nsurv +:= 1; end if;
    end for; end for;
    printf "residue <%o,%o>: #aux7=%o  (R,w) mod 343 total=%o survive=%o capped=%o killed=%o\n",
        R0, w0, #aux7, total, nsurv, ncap, total - nsurv - ncap;
end for;
print "DONE";
quit;
