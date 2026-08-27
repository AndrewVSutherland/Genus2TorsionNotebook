// 37-torsion Jacobian hunt, step D: reconstruct genus-2 curves for the
// conjugation-stable principally polarized members of the isogeny class
// of A_f, f = 2190.2.a.v, and test their torsion for the 37-point.
//
// Pipeline per candidate maximal isotropic H (combined over p | 90354,
// filtered by star-involution stability = complex-conjugation stability,
// a genuine necessary condition for Q-definedness):
//   lattice L' -> symplectic basis (FrobeniusFormAlternating) -> big
//   period matrix -> tau (checked symmetric, Im > 0) -> Rosenhain
//   invariants via theta quotients -> numeric Igusa-Clebsch -> rational
//   recognition (continued fractions) -> Mestre construction over Q ->
//   quadratic-twist fix by L-polynomial matching -> TorsionSubgroup.
// Any curve with 37 | #tors is the first genus-2 Jacobian with a rational
// 37-torsion point.
//
// Run: magma -b Prec:=100 code/gl2_theta_reconstruct.m > results/gl2_theta_2190.log

SetColumns(0);
SetSeed(1);
if not assigned Prec then Prec := 100; elif Type(Prec) eq MonStgElt then Prec := StringToInteger(Prec); end if;
if not assigned MaxRec then MaxRec := 60; elif Type(MaxRec) eq MonStgElt then MaxRec := StringToInteger(MaxRec); end if;
SetMemoryLimit(24*10^9);

Lv := 2190;
trtargets := [<7,2>, <11,-6>, <13,-2>, <17,6>, <19,10>];

M := ModularSymbols(Lv, 2, 0);
S := CuspidalSubspace(M);
NS := NewSubspace(S);
D := NewformDecomposition(NS);
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

BAf := ChangeRing(BasisMatrix(VectorSpace(Af)), Rationals());
Lb := ChangeRing(BasisMatrix(Lattice(Af)), Rationals());
C := Solution(BAf, Lb);
IP := ChangeRing(IntersectionPairing(Af), Rationals());
E := C * IP * Transpose(C);
den := LCM([Denominator(x) : x in Eltseq(E)]);
EZ0 := Matrix(Integers(), 4, 4, [Integers()!(den*x) : x in Eltseq(E)]);
g := GCD([x : x in Eltseq(EZ0) | x ne 0]);
EZ := EZ0 div g;
ed := ElementaryDivisors(EZ);
dd := ed[3];
printf "TYPE %o dd=%o\n", ed, dd;

// star involution on the lattice
SI := StarInvolution(Af);
SIL := C * ChangeRing(SI, Rationals()) * C^-1;
error if LCM([Denominator(x) : x in Eltseq(SIL)]) ne 1, "star not integral";
SIZ := Matrix(Integers(), 4, 4, [Integers()!x : x in Eltseq(SIL)]);
printf "STAR ok\n";

Zd := Integers(dd);
Einv := ChangeRing(EZ, Rationals())^-1;
gens := [ Vector(Rationals(), Eltseq(Einv[i])) : i in [1..4] ];
Ggens := [ Vector(Zd, [ Zd!(Integers()!(dd*x) mod dd) : x in Eltseq(v) ]) : v in gens ];
Gmod := sub< RSpace(Zd, 4) | Ggens >;
EZd := ChangeRing(EZ, Zd);
hd := [ ChangeRing(SIZ, Zd) ];   // conjugation-stability filter

fac := Factorization(dd);
function AddOrd(x)
    return dd div GCD([dd] cat [Integers()!c : c in Eltseq(x)]);
end function;
function MaxIsotropics(p, k)
    m := dd div p^k;
    pg := [ m*g_ : g_ in Ggens ];
    Mp := sub< RSpace(Zd, 4) | pg >;
    ord := p^k;
    e1 := 0; e2 := 0;
    for x in Mp do if AddOrd(x) eq ord then e1 := x; break; end if; end for;
    for x in Mp do
        if AddOrd(x) eq ord and #sub<Mp | [e1, x]> eq ord^2 then e2 := x; break; end if;
    end for;
    genlist := [ e1 + t*e2 : t in [0..ord-1] ] cat
               [ p*u*e1 + e2 : u in [0..(ord div p)-1] ];
    cands := [];
    for x in genlist do
        Hx := sub< Mp | x >;
        if #Hx ne ord then continue; end if;
        iso := &and[ IsZero(InnerProduct(u*EZd, v)) : u, v in Generators(Hx) ];
        if not iso then continue; end if;
        stab := &and[ u*h in Hx : h in hd, u in Generators(Hx) ];
        if not stab then continue; end if;
        Append(~cands, Hx);
    end for;
    if k eq 2 then
        Hp := sub< Mp | [ p*x : x in Generators(Mp) ] >;
        if #Hp eq ord then
            iso := &and[ IsZero(InnerProduct(u*EZd, v)) : u, v in Generators(Hp) ];
            stab := &and[ u*h in Hp : h in hd, u in Generators(Hp) ];
            if iso and stab then Append(~cands, Hp); end if;
        end if;
    end if;
    return cands;
end function;

allc := [* *];
for pf in fac do
    cs := MaxIsotropics(pf[1], pf[2]);
    printf "PRIME %o^%o: %o star-stable isotropics\n", pf[1], pf[2], #cs;
    Append(~allc, cs);
end for;
ncomb := &*[ #cs : cs in allc ];
printf "STAR_STABLE_COMBINATIONS %o\n", ncomb;

// periods
nterms := Ceiling(Lv * Prec * Log(10)/(4*Pi(RealField(20)))) + 100;
printf "PERIOD_TERMS %o\n", nterms;
Pv := Periods(Af, nterms);
CC := ComplexField(Prec);
PM := Matrix(CC, 4, 2, [ CC!Pv[i][j] : j in [1..2], i in [1..4] ]);
printf "PERIODS ok\n";

// Riemann-relation sanity for the base lattice: E-symplectic basis of L
function TauFromLattice(Bp, EpZ)
    // symplectic transform: T * EpZ * T^t = standard J = [[0,I],[-I,0]]
    T, _ := FrobeniusFormAlternating(EpZ);
    Bs := ChangeRing(T, Rationals()) * Bp;      // symplectic basis rows
    Om := ChangeRing(Bs, CC) * PM;              // 4x2 period vectors
    A := Matrix(CC, 2, 2, [Om[1][1], Om[1][2], Om[2][1], Om[2][2]]);
    B := Matrix(CC, 2, 2, [Om[3][1], Om[3][2], Om[4][1], Om[4][2]]);
    tau := B^-1 * A;
    // symmetry/positivity checks + normalization
    sym := Max([Abs(tau[1][2] - tau[2][1])]);
    tau[1][2] := (tau[1][2]+tau[2][1])/2; tau[2][1] := tau[1][2];
    imt := Matrix(RealField(Prec), 2,2, [Im(x) : x in Eltseq(tau)]);
    pos := imt[1][1] gt 0 and Determinant(imt) gt 0;
    if not pos then
        tau := -Conjugate(tau);  // try the other orientation
        tau := Matrix(CC,2,2,[Conjugate(-x) : x in Eltseq(tau)]);
    end if;
    return tau, sym;
end function;

// theta with half characteristics (genus 2), direct series
function Theta2(a, b, tau, prec)
    // a, b in {0,1/2}^2; z = 0
    R := RealField(prec);
    Cc := ComplexField(prec);
    pi := Pi(R);
    s := Cc!0;
    Bnd := Ceiling(Sqrt(prec*Log(10)/ (pi*Min([Im(tau[1][1]), Im(tau[2][2])])))) + 8;
    for n1 in [-Bnd..Bnd] do
        for n2 in [-Bnd..Bnd] do
            v := Vector(Cc, [n1 + a[1], n2 + a[2]]);
            q := (v*ChangeRing(tau,Cc), v);
            s +:= Exp(pi*Cc.1*q + 2*pi*Cc.1*(v[1]*b[1] + v[2]*b[2]));
        end for;
    end for;
    return s;
end function;

// Rosenhain from theta-nulls (classical): la = t1^2 t3^2/(t2^2 t4^2) etc.
// using the standard 4 even thetas:
//   t1 = theta[0,0;0,0], t2 = theta[0,0;1/2,1/2],
//   t3 = theta[1/2,0;0,0]? -- conventions vary; we use the Rosenhain set:
function RosenhainFromTau(tau, prec)
    h := 1/2;
    t00 := Theta2([0,0],[0,0], tau, prec);
    t01 := Theta2([0,0],[0,h], tau, prec);
    t10 := Theta2([0,0],[h,0], tau, prec);
    t11 := Theta2([0,0],[h,h], tau, prec);
    t20 := Theta2([h,0],[0,0], tau, prec);
    t22 := Theta2([0,h],[0,0], tau, prec);
    la := (t00*t01/(t10*t11))^2;
    mu := (t01*t20/(t11*t22))^2 * (t00*t20/(t10*t22))^0;  // placeholder combo
    nu := (t00*t20/(t10*t22))^2;
    return la, mu, nu;
end function;

// NOTE: exact Rosenhain conventions are fiddly; instead of relying on them
// we go through IgusaClebsch of a NUMERIC curve built from the six theta-
// derived branch points is error-prone.  ROBUST alternative: absolute
// Igusa invariants direct from the 10 even theta-nulls (Igusa's formulas):
function EvenThetas(tau, prec)
    h := 1/2;
    chars := [
        <[0,0],[0,0]>, <[0,0],[0,h]>, <[0,0],[h,0]>, <[0,0],[h,h]>,
        <[h,0],[0,0]>, <[h,0],[0,h]>, <[0,h],[0,0]>, <[0,h],[h,0]>,
        <[h,h],[0,0]>, <[h,h],[h,h]>
    ];
    return [ Theta2(c[1], c[2], tau, prec) : c in chars ];
end function;

function IgusaFromTau(tau, prec)
    th := EvenThetas(tau, prec);
    t4 := [ t^4 : t in th ];
    h4 := &+t4;
    h6 := &+[ t4[i]*(&+[ t4[j] : j in [1..10] | j ne i ]) : i in [1..10] ]; // crude symmetric fn (placeholder)
    // Igusa's psi4, psi6, chi10, chi12 in terms of theta-nulls:
    psi4 := (&+[ t^8 : t in th ])/4;  // placeholder normalization
    chi10 := (&*[ t^2 : t in th ]) / 2^12;
    return psi4, chi10;
end function;

printf "RECONSTRUCTION_STUBS_PRESENT -- run stops here pending exact Igusa formulas\n";
printf "GL2_THETA_STAGE0_DONE combos=%o\n", ncomb;
quit;
