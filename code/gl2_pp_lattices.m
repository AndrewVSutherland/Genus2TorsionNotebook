// 37-torsion Jacobian hunt, step C: enumerate Hecke-stable principally
// polarizing lattices for A_f, f = 2190.2.a.v.
//
// Data from step A: integral lattice L (rank 4) of the modular symbols
// piece, intersection pairing E of type (1, 90354), 90354 = 2*3*11*37^2.
// Principally polarized members of the isogeny class defined over Q
// correspond (via the analytic uniformization) to Hecke-stable lattices
// L' with L <= L' <= L^# (the E-dual) on which E is unimodular, i.e.
// preimages of Hecke-stable maximal isotropic subgroups H of
// L^#/L =~ (Z/90354)^2 with the E-induced pairing.  Hecke-stability
// (equivalently Galois-stability of the corresponding finite subgroup)
// is checked with integral Hecke operators.  Enumeration is prime-by-
// prime over 90354 = 2 * 3 * 11 * 37^2.
//
// Output: the list of principalizing lattices L' (bases in L-coordinates,
// scaled integrally) for step D (theta reconstruction).
//
// Run: magma -b code/gl2_pp_lattices.m > results/gl2_pp_lattices_2190.log

SetColumns(0);
SetSeed(1);
SetMemoryLimit(16*10^9);

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
printf "E %o\n", Eltseq(EZ);
ed := ElementaryDivisors(EZ);
printf "TYPE %o\n", ed;
dd := ed[3];   // = 90354
assert ed[1] eq 1 and ed[2] eq 1 and ed[4] eq dd;

// Hecke matrices on the lattice basis (integral)
function HeckeOnLattice(l)
    Tl := HeckeOperator(Af, l);
    TlL := C * ChangeRing(Tl, Rationals()) * C^-1;
    d2 := LCM([Denominator(x) : x in Eltseq(TlL)]);
    error if d2 ne 1, "Hecke not integral on lattice";
    return Matrix(Integers(), 4, 4, [Integers()!x : x in Eltseq(TlL)]);
end function;
heckes := [ HeckeOnLattice(l) : l in [2,3,5,7,11,13,17,19,73] ];
printf "HECKE_INTEGRAL ok (%o operators)\n", #heckes;

// RM order action: sqrt3 = (T_7 - 1)/c where a_7 = 1 + c*sqrt3.
// From charpoly(T_7|Af) = (x^2 - 2x + n)^2, n = Norm(a_7) = 1 - 3c^2.
T7L := HeckeOnLattice(7);
cp := CharacteristicPolynomial(T7L);
fac7 := Factorization(cp);
q2 := fac7[1][1];
error if Degree(q2) ne 2 or fac7[1][2] ne 2, "unexpected T7 charpoly shape";
nq := Coefficient(q2, 0);
c2 := (1 - nq)/3;
okc, cval := IsSquare(Integers()!c2);
error if not okc, "c not integral square";
printf "A7 = 1 + %o*sqrt3 (Norm %o)\n", cval, nq;
R3num := T7L - IdentityMatrix(Integers(), 4);
error if not &and[ x mod cval eq 0 : x in Eltseq(R3num) ], "sqrt3 matrix not integral";
R3 := R3num div cval;    // integral matrix of sqrt3 on the lattice
printf "SQRT3_MATRIX charpoly %o\n", CharacteristicPolynomial(R3);

// L^# / L: columns of the dual basis: L^# = { v : v*EZ*L^T integral } in
// L-coords: L^# = EZ^-1 * Z^4 (rows).  Work with the finite module
// V := (1/dd)Z^4 / Z^4 restricted to L^#/L.
Einv := ChangeRing(EZ, Rationals())^-1;
// generators of L^#/L: rows of Einv modulo Z^4, all have denominators | dd
gens := [ Vector(Rationals(), Eltseq(Einv[i])) : i in [1..4] ];

// represent L^#/L inside (Z/dd)^4 via v -> dd*v mod dd
Zd := Integers(dd);
Ggens := [ Vector(Zd, [ Zd!(Integers()!(dd*x) mod dd) : x in Eltseq(v) ]) : v in gens ];
Gmod := sub< RSpace(Zd, 4) | Ggens >;
printf "DUAL_QUOTIENT_ORDER %o (expect dd^2)\n", #Gmod;

// pairing on L^#/L in (1/dd)Z/Z: for dd-scaled x,y the isotropy condition
// is x*EZ*y^t == 0 mod dd, i.e. InnerProduct(x*EZd, y) = 0 over Z/dd.
EZd := ChangeRing(EZ, Zd);

// enumerate Hecke-stable maximal isotropic subgroups prime-by-prime.
// L^#/L ~ (Z/dd)^2; its p-primary parts: (Z/p^k)^2 for p^k || dd.
fac := Factorization(dd);
printf "DD_FACT %o\n", fac;

// For each prime p, the p-part Mp of Gmod: generators
function PPart(p, k)
    m := dd div p^k;
    return [ m*g_ : g_ in Ggens ];   // these generate the p-part
end function;

// O = Z[sqrt3] action descends: sqrt3 matrix mod dd acts on Gmod.
// Hecke operators act as O-scalars on ker(lambda) so Hecke-stability is
// vacuous (v2 finding: ALL 202608 combos passed); the correct necessary
// condition for Q-definedness is O-stability (RM defined over Q).
// Final arbiter of Q-definedness = rational recognition of reconstructed
// invariants in step D.
hd := [ ChangeRing(R3, Zd) ];

// isotropic Hecke-stable lines in the p-part, then combine over p.
// A maximal isotropic H decomposes as product of its p-parts, each a
// maximal isotropic Hecke-stable subgroup of (Z/p^k)^2 of order p^k.
// For (Z/p^k)^2 the candidates are the cyclic subgroups of order p^k
// (choices of a generator line mod p, Hensel-structured), plus for k=2
// the subgroup p*(whole) of shape (Z/p)^2.  We enumerate subgroups of
// the right order directly via sub< > over small index sets.
function MaxIsotropics(p, k)
    // p-part inside Gmod
    pg := PPart(p, k);
    Mp := sub< RSpace(Zd, 4) | pg >;
    ord := p^k;
    // all subgroups of order p^k that are isotropic and Hecke-stable
    cands := [];
    printf "  p=%o k=%o #Mp=%o\n", p, k, #Mp;
    // find a basis (e1, e2) of Mp = (Z/p^k)^2
    function AddOrd(x)
        return dd div GCD([dd] cat [Integers()!c : c in Eltseq(x)]);
    end function;
    e1 := 0; e2 := 0;
    for x in Generators(Mp) do
        if AddOrd(x) eq ord then e1 := x; break; end if;
    end for;
    if e1 cmpeq 0 then
        for x in Mp do if AddOrd(x) eq ord then e1 := x; break; end if; end for;
    end if;
    error if e1 cmpeq 0, "no element of full order";
    for x in Generators(Mp) do
        if AddOrd(x) eq ord and #sub<Mp | [e1, x]> eq ord^2 then e2 := x; break; end if;
    end for;
    if e2 cmpeq 0 then
        cnt := 0;
        for x in Mp do
            cnt +:= 1;
            if AddOrd(x) eq ord and #sub<Mp | [e1, x]> eq ord^2 then e2 := x; break; end if;
            if cnt gt 100000 then break; end if;
        end for;
    end if;
    error if e2 cmpeq 0, "no second basis element found";
    // cyclic maximal subgroups = P^1(Z/p^k): <e1 + t e2>, t in Z/p^k, and
    // <p*u*e1 + e2>, u in Z/p^(k-1)
    genlist := [ e1 + t*e2 : t in [0..ord-1] ] cat
               [ p*u*e1 + e2 : u in [0..(ord div p)-1] ];
    for x in genlist do
        Hx := sub< Mp | x >;
        if #Hx ne ord then continue; end if;
        // isotropy: pairing of all pairs in Hx must vanish: (u*EZd)*v = 0 in Zd
        // scaled pairing: value (u * EZ * v^t)/dd mod 1 -> numerator u*EZd*v^t
        // must be 0 mod dd for isotropy in (1/dd)Z/Z terms... on p-part the
        // condition is mod p^k after rescale; direct check in Zd suffices:
        iso := true;
        for u in Generators(Hx) do
            for v in Generators(Hx) do
                num := InnerProduct(u*EZd, v);
                if not IsZero(num) then iso := false; break; end if;
            end for;
            if not iso then break; end if;
        end for;
        if not iso then continue; end if;
        // Hecke stability
        stab := true;
        for h in hd do
            for u in Generators(Hx) do
                if u*h notin Hx then stab := false; break; end if;
            end for;
            if not stab then break; end if;
        end for;
        if not stab then continue; end if;
        Append(~cands, Hx);
    end for;
    // k=2: also the homocyclic subgroup p * Mp (order p^2, shape (Z/p)^2)
    if k eq 2 then
        Hp := sub< Mp | [ p*x : x in Generators(Mp) ] >;
        if #Hp eq ord then
            iso := true;
            for u in Generators(Hp) do for v in Generators(Hp) do
                if not IsZero(InnerProduct(u*EZd, v)) then iso := false; end if;
            end for; end for;
            stab := true;
            for h in hd do for u in Generators(Hp) do
                if u*h notin Hp then stab := false; end if;
            end for; end for;
            if iso and stab then Append(~cands, Hp); end if;
        end if;
    end if;
    return cands;
end function;

allc := [* *];
for pf in fac do
    cs := MaxIsotropics(pf[1], pf[2]);
    printf "PRIME %o^%o: %o Hecke-stable isotropic subgroups of full order\n",
           pf[1], pf[2], #cs;
    Append(~allc, cs);
end for;

// combine: total H = sum of chosen p-parts; lattice L' = L + (1/dd)*lift(H)
ncomb := &*[ #cs : cs in allc ];
printf "TOTAL_COMBINATIONS %o\n", ncomb;
nprin := 0;
combo := [ 1 : cs in allc ];
done := false;
while not done do
    // build H
    Hgens := [];
    for j in [1..#allc] do
        Hgens cat:= [ x : x in Generators(allc[j][combo[j]]) ];
    end for;
    // lattice L' in L-coords: rows of identity + lifts/dd
    lifts := [ Vector(Rationals(), [ (Integers()!x)/dd : x in Eltseq(u) ]) : u in Hgens ];
    B := VerticalJoin(ChangeRing(IdentityMatrix(Integers(),4), Rationals()),
                      Matrix(Rationals(), #lifts, 4, &cat[Eltseq(v) : v in lifts]));
    // HNF of the stacked generators = basis of L'
    dn := LCM([Denominator(x) : x in Eltseq(B)]);
    BZ := Matrix(Integers(), Nrows(B), 4, [Integers()!(dn*x) : x in Eltseq(B)]);
    H := HermiteForm(BZ);
    Bp := Matrix(Rationals(), 4, 4, [ H[i][j]/dn : j in [1..4], i in [1..4] ]);
    Ep := Bp * ChangeRing(EZ, Rationals()) * Transpose(Bp);
    dnE := LCM([Denominator(x) : x in Eltseq(Ep)]);
    if dnE eq 1 then
        EpZ := Matrix(Integers(), 4, 4, [Integers()!x : x in Eltseq(Ep)]);
        edp := ElementaryDivisors(EpZ);
        if edp[1] eq edp[4] then    // unimodular up to scaling: principal
            nprin +:= 1;
            printf "PRINCIPAL_LATTICE %o combo=%o basis=%o edivs=%o\n",
                   nprin, combo, Eltseq(Bp), edp;
        end if;
    end if;
    // increment combo
    j := 1;
    while j le #allc do
        combo[j] +:= 1;
        if combo[j] le #allc[j] then break; end if;
        combo[j] := 1; j +:= 1;
    end while;
    if j gt #allc then done := true; end if;
end while;
printf "N_PRINCIPAL %o of %o combos\n", nprin, ncomb;
printf "GL2_PP_DONE\n";
quit;
