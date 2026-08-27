// 37-torsion Jacobian hunt, step A: period lattice + polarization type of
// the modular abelian surface A_f, f = 2190.2.a.v.
// Computes: the 4 period vectors in C^2 (Periods on the modular symbols
// piece), and the intersection pairing on the integral lattice with its
// elementary divisors = the type of the modular polarization.  If the type
// is (d,d) the polarization is principal up to scaling; type (d1,d2) with
// d1 != d2 means genuine (d1/d,d2/d)-polarization and step C (symplectic
// modification / 37-quotient) follows.
//
// Run: magma -b Lv:=2190 Prec:=60 code/gl2_af_periods.m > results/gl2_periods_2190.log

SetColumns(0);
SetSeed(1);
if not assigned Lv then Lv := 2190; elif Type(Lv) eq MonStgElt then Lv := StringToInteger(Lv); end if;
if not assigned Prec then Prec := 60; elif Type(Prec) eq MonStgElt then Prec := StringToInteger(Prec); end if;
if not assigned MemGB then MemGB := 24; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

if Lv eq 2190 then
    trtargets := [<7,2>, <11,-6>, <13,-2>, <17,6>, <19,10>];
elif Lv eq 1830 then
    trtargets := [<7,0>, <11,8>, <13,-2>, <17,6>, <19,-6>];
else
    error "unknown level";
end if;

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

// intersection pairing on the integral lattice of Af
IP := IntersectionPairing(Af);
printf "IP_COMPUTED size=%o\n", Nrows(IP);
BAf := BasisMatrix(VectorSpace(Af));
Lb := ChangeRing(BasisMatrix(Lattice(Af)), Rationals());
C := Solution(ChangeRing(BAf, Rationals()), Lb);
E := C * ChangeRing(IP, Rationals()) * Transpose(C);
// scale to primitive integral alternating matrix
den := LCM([Denominator(x) : x in Eltseq(E)]);
EZ := Matrix(Integers(), 4, 4, [Integers()!(den*x) : x in Eltseq(E)]);
g := GCD([x : x in Eltseq(EZ) | x ne 0]);
EZ := EZ div g;
printf "PAIRING_MATRIX %o\n", Eltseq(EZ);
ediv := ElementaryDivisors(EZ);
printf "ELEMENTARY_DIVISORS %o\n", ediv;
printf "POLARIZATION_TYPE (%o, %o) up to scaling\n", ediv[1], ediv[3];

// number of q-expansion terms for the target precision
nterms := Ceiling(Lv * Prec * Log(10)/(4*Pi(RealField(20)))) + 100;
printf "PERIOD_TERMS %o\n", nterms;
P := Periods(Af, nterms);
printf "N_PERIOD_VECTORS %o\n", #P;
CC := ComplexField(Prec);
for i in [1..#P] do
    printf "PERIOD %o: %o\n", i, [CC!P[i][j] : j in [1..Degree(Universe(P))]];
end for;
printf "GL2_PERIODS_DONE level=%o\n", Lv;
quit;
