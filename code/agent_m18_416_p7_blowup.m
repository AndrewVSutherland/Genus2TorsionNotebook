
//////////////////////////////////////////////////////////////////////
//  p = 7 blowup / Hensel tower for the [4,16] second-halving cover in
//  the M_1(8,4) family.
//
//  Background (see notes/agent_m18_416_search_notes.md,
//  code/m18_m14_second_halving_equations.m):
//    C : y^2 = f = x*A*B,   parameters (R,w),
//    t = (2R^2+(1-w^2)R-2w^2)/(4(w^2-1)),
//    A = x^2+(R^3+4R^2 t+R-8R t+4t)x+R^4,
//    B = (R+2+4t)x^2+(R^2+4R+1+8t)x+(2R^2+R+4t),
//    c4 = R+2+4t = lc(A*B).
//  The order-8 point P_R=(-R,Y_R) is divisible by 2 (=> torsion contains
//  Z/16) iff  f - ell^2 = c4*(x+R)*q^2  has a rational solution
//    q = x^2+a*x+b,  ell = c*x^2+d*x+e.
//  These give the 5 cleared equations E416(R,w,a,b,c,d,e) (x^0..x^4).
//
//  Finite diagnostic (Filip): the good-open [4,16] locus is empty over
//  F_7, so any rational [4,16] must reduce onto the p=7 boundary.  This
//  script performs the p-adic BLOWUP that the boundary CRT/height searches
//  do not: for each mod-7 solution it linearises the E416 system and finds
//  which next 7-adic digits of (R,w) survive, i.e. the mod-49 (and then
//  mod-343) closure.  Shrinking-to-empty ==> a genuine 7-adic obstruction;
//  a stable nonempty closure ==> concrete higher-precision residue targets
//  for a refined search.
//
//  Method.  Work with cleared INTEGER polynomials E416_i.  A point
//  base=(R,w,a,b,c,d,e) in (Z/7^k)^7 with E416(base)=0 mod 7^k lifts to
//  base + 7^k*delta (delta in F_7^7) with E416=0 mod 7^{k+1} iff
//     carry_i + sum_j J_ij delta_j = 0  (mod 7),
//  where J_ij = dE416_i/dvar_j (base) mod 7 and carry_i = E416_i(base)/7^k
//  mod 7.  This is a 5x7 linear system over F_7; its solution set is an
//  affine subspace.  We track the (R,w) coordinates of the lifts.
//
//  Level:=k controls how many blowup levels (mod 7^{k+1}) to enumerate.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();

Level := 2;
if assigned level then Level := StringToInteger(level); end if;
p := 7;
if assigned prime then p := StringToInteger(prime); end if;

// ---- build cleared integer equations E416 in (R,w,a,b,c,d,e) ----
Rng<R,w,a,b,c,d,e> := PolynomialRing(Q, 7, "grevlex");
KF := FieldOfFractions(Rng);
PX<x> := PolynomialRing(KF);
t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
A := x^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*x + R^4;
B := (R + 2 + 4*t)*x^2 + (R^2 + 4*R + 1 + 8*t)*x + (2*R^2 + R + 4*t);
f := x*A*B;
c4 := R + 2 + 4*t;
q416 := x^2 + a*x + b;
ell416 := c*x^2 + d*x + e;
F416 := f - ell416^2 - c4*(x+R)*q416^2;
E416 := [];
for i in [0..4] do
    ci := Numerator(Coefficient(F416, i));
    den := LCM([Denominator(cc) : cc in Coefficients(ci)]);
    ci := den*ci;                                  // integer coefficients
    g := GCD([Integers()!cc : cc in Coefficients(ci)]);
    ci := Rng!(ci / g);                            // primitive integer poly
    Append(~E416, ci);
end for;

// precompute formal derivatives dE416_i/dvar_j (var order R,w,a,b,c,d,e)
nv := 7;
Der := [[Derivative(E416[i], j) : j in [1..nv]] : i in [1..5]];

// integer evaluation of a poly at an integer vector
function Ev(pol, v)  // v is a sequence of 7 integers
    return Evaluate(pol, v);
end function;

Fp := GF(p);

// enumerate F_p solutions (a,b,c,d,e) of E416 at a fixed (R0,w0) in F_p
function SolAux(R0, w0)
    sols := [];
    // substitute R,w = R0,w0 (as integers), reduce mod p, brute (a..e)
    for aa in [0..p-1] do for bb in [0..p-1] do for cc in [0..p-1] do
      for dd in [0..p-1] do for ee in [0..p-1] do
        v := [Integers()!R0, Integers()!w0, aa, bb, cc, dd, ee];
        ok := true;
        for i in [1..5] do
            if (Integers()!Ev(E416[i], v)) mod p ne 0 then ok := false; break; end if;
        end for;
        if ok then Append(~sols, [aa,bb,cc,dd,ee]); end if;
      end for; end for; end for;
    end for; end for;
    return sols;
end function;

// ---- diagnostic: good-open vs boundary E416 solvability over F_p ----
boundaryFactors := [
  R, w, w-1, w+1, R-1, R+1, R-w, R+w,
  R*w-3*R+3*w-1, R*w+3*R+3*w+1,
  2*R^2-R*w^2+R-2*w^2, R^4-2*R^3+R^2*w^2-R^2+2*R*w^2-w^2 ];
function IsBoundary(R0,w0)
    v := [Integers()!R0,Integers()!w0,0,0,0,0,0];
    for bf in boundaryFactors do
        if (Integers()!Ev(bf,v)) mod p eq 0 then return true; end if;
    end for;
    return false;
end function;

print "M18_416_P7_BLOWUP p", p, "Level", Level;
goodopen_solvable := 0; boundary_solvable := 0;
targetResidues := [];   // (R0,w0) with E416 solvable
for R0 in [0..p-1] do for w0 in [0..p-1] do
    s := SolAux(R0,w0);
    if #s eq 0 then continue; end if;
    if IsBoundary(R0,w0) then
        boundary_solvable +:= 1;
        Append(~targetResidues, [R0,w0]);
    else
        goodopen_solvable +:= 1;
        Append(~targetResidues, [R0,w0]);
    end if;
end for; end for;
printf "E416 solvable residues: good_open=%o boundary=%o (good_open should be 0)\n",
    goodopen_solvable, boundary_solvable;

// ---- first blowup (mod p -> mod p^2): which (R,w) mod p^2 survive ----
// At a boundary residue (R0,w0) with a mod-p aux solution s=(a,b,c,d,e),
// lift R=R0+p*rho, w=w0+p*omega.  The aux next-digit delta in F_p^5 must
// solve  J_aux * delta = -carry(rho,omega)  where
//   J_aux = dE416/d(a,b,c,d,e)(base)  (5x5 over F_p),
//   carry_i(rho,omega) = c0_i + gR_i*rho + gw_i*omega   (mod p),
//   c0_i = E416_i(R0,w0,s)/p mod p,  gR_i=dE416_i/dR, gw_i=dE416_i/dw at base.
// Solvable  <=>  P*carry = 0, P = left-null (cokernel) basis of J_aux.
// That is a linear condition on (rho,omega); its solutions are the live
// tangent directions.  Union over all mod-p aux solutions s.
//
// var index in Rng: R=1,w=2,a=3,b=4,c=5,d=6,e=7.
auxidx := [3,4,5,6,7];

survivors := AssociativeArray();   // key "R0,w0" -> set of <rho,omega> live
smoothAt := AssociativeArray();     // key -> count of Hensel-smooth aux solutions
allrhoom := {<r,o> : r in [0..p-1], o in [0..p-1]};
for rw in targetResidues do
    R0 := rw[1]; w0 := rw[2];
    live := {};
    nsmooth := 0;
    for s in SolAux(R0, w0) do
        base := [R0, w0] cat s;
        // J_aux (5x5)
        Jaux := Matrix(Fp, 5, 5,
            [ Fp!(Integers()!Ev(Der[i][auxidx[j]], base)) : i in [1..5], j in [1..5] ]);
        // cokernel P = left kernel of Jaux (rows v with v*Jaux=0)
        P := KernelMatrix(Jaux);
        if Nrows(P) eq 0 then
            // J_aux invertible => s Hensel-lifts to a SMOOTH Q_7 point of E416
            // for EVERY (rho,omega): a 7-adic [4,16] with bad curve reduction
            // at 7.  These are the prime candidates for a rational point.
            nsmooth +:= 1;
            live := live join allrhoom;
            continue;
        end if;
        // c0, gR, gw as F_p^5 vectors
        c0 := [ Fp!((Integers()!Ev(E416[i], base)) div p) : i in [1..5] ];
        gR := [ Fp!(Integers()!Ev(Der[i][1], base)) : i in [1..5] ];
        gw := [ Fp!(Integers()!Ev(Der[i][2], base)) : i in [1..5] ];
        // conditions rows of P applied: (P*c0) + (P*gR)*rho + (P*gw)*omega = 0
        nr := Nrows(P);
        Pc  := [ &+[P[r][j]*c0[j]  : j in [1..5]] : r in [1..nr] ];
        PgR := [ &+[P[r][j]*gR[j]  : j in [1..5]] : r in [1..nr] ];
        Pgw := [ &+[P[r][j]*gw[j]  : j in [1..5]] : r in [1..nr] ];
        for rho in [0..p-1] do for om in [0..p-1] do
            ok := true;
            for r in [1..nr] do
                if Pc[r] + PgR[r]*(Fp!rho) + Pgw[r]*(Fp!om) ne 0 then ok := false; break; end if;
            end for;
            if ok then Include(~live, <rho,om>); end if;
        end for; end for;
    end for;
    survivors[Sprintf("%o,%o", R0, w0)] := live;
    smoothAt[Sprintf("%o,%o", R0, w0)] := nsmooth;
end for;

// report
total49 := 0; nSmoothResidues := 0;
print "FIRST_BLOWUP mod", p^2, ": surviving (R,w) mod", p^2, "per boundary residue";
print "(smooth = # mod-7 aux solutions with invertible aux-Jacobian => Hensel-liftable to Q_7)";
for rw in targetResidues do
    R0 := rw[1]; w0 := rw[2];
    live := survivors[Sprintf("%o,%o", R0, w0)];
    ns := smoothAt[Sprintf("%o,%o", R0, w0)];
    if ns gt 0 then nSmoothResidues +:= 1; end if;
    printf "  residue <%o,%o>: surviving = %o / %o   Hensel_smooth_aux = %o%o\n",
        R0, w0, #live, p^2, ns, (ns gt 0 select "   <== 7-adic [4,16] candidate" else "");
    total49 +:= #live;
end for;
printf "Hensel-smooth boundary residues (prime candidates): %o\n", nSmoothResidues;
printf "TOTAL surviving (R,w) mod %o = %o  (was %o residues mod %o, i.e. %o lifts before blowup)\n",
    p^2, total49, #targetResidues, p, #targetResidues*p^2;
if total49 eq 0 then
    print "OBSTRUCTION: first blowup kills every boundary lift => no rational [4,16] in M_1(8,4) reduces mod 7 at all.";
else
    print "Live tangent strata exist; these mod-49 residues are the refined search targets.";
end if;
print "DONE";
quit;
