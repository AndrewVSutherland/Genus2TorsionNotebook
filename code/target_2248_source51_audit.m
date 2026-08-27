//////////////////////////////////////////////////////////////////////
// Independent audit of the two [2,4,8] row-51 Richelot sources.
//
// This deliberately does not use the v_2 reduction gate to decide any
// torsion group.  Every Jacobian neighbor is sent to TorsionSubgroup.
// It also tests geometric simplicity on the common, small row-51 base:
// simplicity/nonsimplicity is invariant under isogeny, so this is both
// cleaner and stronger than trying to find a witness on the huge dual
// models.
//
// Run from torsion_jac:
//   magma -b code/target_2248_source51_audit.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

if not assigned bank_file then
    bank_file := "data/tor2244_bank.txt";
end if;
if not assigned write_log then
    write_log := true;
elif Type(write_log) eq MonStgElt then
    write_log := write_log in {"true","True","1","yes"};
end if;
if not assigned log_file then
    log_file := "results/target_2248_source51_audit.log";
end if;
if write_log then SetLogFile(log_file : Overwrite := true); end if;

PrimeList := [
    3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,
    73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,
    151,157,163,167,173,179,181,191,193,197,199
];

function BoolString(b)
    return b select "true" else "false";
end function;

function ReadTupleFile(filename)
    rows := Split(Read(filename), "\n");
    tuples := [];
    for raw in rows do
        n := #raw;
        if n gt 0 and raw[n] eq "\r" then n -:= 1; end if;
        if n ge 2 and raw[1] eq "[" and raw[n] eq "]" then
            ss := Split(raw[2..n-1], ",");
            if #ss eq 4 then
                Append(~tuples,[StringToInteger(s) : s in ss]);
            end if;
        end if;
    end for;
    return tuples;
end function;

function SquareIntegralModel(f)
    d := LCM([Denominator(Coefficient(f,i)) : i in [0..Degree(f)]]);
    return P!(d^2*f);
end function;

function CompletedPolynomial(C)
    f,h := HyperellipticPolynomials(C);
    return h eq 0 select P!f else P!(h^2+4*f);
end function;

function NormalizeCurve(C)
    D := SimplifiedModel(HyperellipticCurve(SquareIntegralModel(
         CompletedPolynomial(C))));
    fD,hD := HyperellipticPolynomials(D);
    return D,P!fD,P!hD;
end function;

function FactorDegrees(C)
    F := CompletedPolynomial(C);
    return Sort(&cat[[Degree(z[1]) : j in [1..z[2]]]
                     : z in Factorization(F)]);
end function;

function FullWeierstrass(C)
    ds := FactorDegrees(C);
    return (Degree(CompletedPolynomial(C)) eq 5 and
            ds eq [1,1,1,1,1]) or
           (Degree(CompletedPolynomial(C)) eq 6 and
            ds eq [1,1,1,1,1,1]);
end function;

function SameQCurve(C,D)
    try
        ok,mp := IsIsomorphic(C,D);
        return ok;
    catch e
        return false;
    end try;
end function;

function AbsoluteInvariants(C)
    return [Q!z : z in G2Invariants(C)];
end function;

// This is the workspace's strict root-power standard: irreducible
// Frobenius quartic and no degree drop for pi^n, 2 <= n <= 12.
function RootPowerAtPrime(f,pp)
    disc := Discriminant(f);
    if Z!Denominator(LeadingCoefficient(f)) mod pp eq 0 or
       Z!Numerator(LeadingCoefficient(f)) mod pp eq 0 or
       Z!Denominator(disc) mod pp eq 0 or
       Z!Numerator(disc) mod pp eq 0 then
        return false,P!0,[];
    end if;
    Fp := GF(pp); Pp := PolynomialRing(Fp);
    fp := Pp![Fp!c : c in Coefficients(f)];
    if Degree(fp) ne Degree(f) or not IsSquarefree(fp) then
        return false,P!0,[];
    end if;
    Lp := LPolynomial(HyperellipticCurve(fp));
    cs := [Z!Coefficient(Lp,i) : i in [0..4]];
    R<T> := PolynomialRing(Q);
    chi := &+[Q!cs[5-i]*T^i : i in [0..4]];
    if not IsIrreducible(chi) then return false,P!chi,[]; end if;
    K<pi> := NumberField(chi);
    degrees := [Degree(MinimalPolynomial(pi^n)) : n in [2..12]];
    return &and[d eq 4 : d in degrees],P!chi,degrees;
end function;

function RootPowerCertificates(f)
    certs := [];
    tested := 0;
    for pp in PrimeList do
        ok,chi,degrees := RootPowerAtPrime(f,pp);
        if chi ne 0 then tested +:= 1; end if;
        if ok then Append(~certs,<pp,chi,degrees>); end if;
    end for;
    return certs,tested;
end function;

function Bracket(A,B)
    return Derivative(A)*B-A*Derivative(B);
end function;

function CoefficientMatrixDeterminant(A,B,C)
    return Determinant(Matrix(Q,3,3,
        [Coefficient(A,i) : i in [0..2]] cat
        [Coefficient(B,i) : i in [0..2]] cat
        [Coefficient(C,i) : i in [0..2]]));
end function;

function Pairings4()
    return [
        [<1,2>,<3,4>],
        [<1,3>,<2,4>],
        [<1,4>,<2,3>]
    ];
end function;

// For a sextic with one quadratic and four linear factors, return all three
// rational maximal isotropic kernels.  The leading coefficient is absorbed
// into A.  This correction is essential: using monic factors without it can
// silently replace the Richelot codomain by a quadratic twist.
function OneSplitKernelTriples(C)
    F := CompletedPolynomial(C);
    fac := Factorization(F);
    lins := [z[1] : z in fac | Degree(z[1]) eq 1];
    quads := [z[1] : z in fac | Degree(z[1]) eq 2];
    assert #lins eq 4 and #quads eq 1 and Degree(F) eq 6;
    out := [];
    for pairing in Pairings4() do
        A := LeadingCoefficient(F)*quads[1];
        B := lins[pairing[1][1]]*lins[pairing[1][2]];
        D := lins[pairing[2][1]]*lins[pairing[2][2]];
        label := Sprintf("Q | L%o*L%o | L%o*L%o",
            pairing[1][1],pairing[1][2],
            pairing[2][1],pairing[2][2]);
        assert A*B*D eq F;
        Append(~out,<A,B,D,label>);
    end for;
    return out;
end function;

function ManualRichelotCurve(A,B,D)
    delta := CoefficientMatrixDeterminant(A,B,D);
    if delta eq 0 then return false,_,_,delta,P!0; end if;
    g := delta*Bracket(B,D)*Bracket(D,A)*Bracket(A,B);
    if Degree(g) lt 5 or Discriminant(g) eq 0 then
        return false,_,_,delta,g;
    end if;
    C,f,h := NormalizeCurve(HyperellipticCurve(g));
    return true,C,f,delta,g;
end function;

print "TARGET_2248_SOURCE51_INDEPENDENT_AUDIT";
print "SCOPE row 51; reverse edges 13 and 14; every rational forward Richelot kernel; exact torsion on every Jacobian codomain";

rows := ReadTupleFile(bank_file);
assert #rows ge 51;
t := rows[51];
assert t eq [55216,56550,62234,64090];
a,b,c,d := Explode(t);
fbase := x*(x+a^2)*(x+b^2)*(x+c^2)*(x+d^2);
Cbase := HyperellipticCurve(fbase);
Gbase,mpbase := TorsionSubgroup(Jacobian(Cbase));
print "BASE_ROW",51,"tuple",t;
print "BASE_POLYNOMIAL",fbase;
print "BASE_FACTORIZATION",Factorization(fbase);
print "BASE_EXACT_TORSION",Invariants(Gbase),"order",#Gbase;

certs,good_tested := RootPowerCertificates(fbase);
print "BASE_STRICT_ROOT_POWER_GOOD_PRIMES_TESTED",good_tested;
print "BASE_STRICT_ROOT_POWER_WITNESSES",#certs,certs;

Aut := AutomorphismGroup(Cbase);
subcovers := Degree2Subcovers(Cbase);
print "BASE_AUTOMORPHISM_GROUP_ORDER",#Aut;
print "BASE_DEGREE2_SUBCOVER_COUNT",#subcovers;
for i in [1..#subcovers] do
    E := subcovers[i][1];
    phi := subcovers[i][2];
    Emin := MinimalModel(E);
    print "BASE_DEGREE2_SUBCOVER",i;
    print "ELLIPTIC_RAW",E;
    print "ELLIPTIC_MINIMAL",Emin;
    print "ELLIPTIC_MINIMAL_AINVARIANTS",aInvariants(Emin);
    print "DEGREE2_MAP",phi;
end for;
assert #subcovers eq 2;
print "BASE_GEOMETRIC_SIMPLICITY_REFUTED true reason two rational degree-2 elliptic subcovers";

Rbase := RichelotIsogenousSurfaces(Jacobian(Cbase));
assert #Rbase eq 15;
sources := [];
source_fs := [];
source_hs := [];
source_edges := [13,14];
for e in source_edges do
    assert Type(Rbase[e]) eq JacHyp;
    C,f,h := NormalizeCurve(Curve(Rbase[e]));
    G,mp := TorsionSubgroup(Jacobian(C));
    print "SOURCE","reverse_edge",e;
    print "SOURCE_POLYNOMIALS",f,h;
    print "SOURCE_FACTORIZATION",Factorization(CompletedPolynomial(C));
    print "SOURCE_FACTOR_DEGREES",FactorDegrees(C);
    print "SOURCE_EXACT_TORSION",Invariants(G),"order",#G;
    print "SOURCE_ABSOLUTE_G2",AbsoluteInvariants(C);
    print "SOURCE_GEOMETRICALLY_NONSIMPLE_BY_ISOGENY true common_base_row 51";
    assert Invariants(G) eq [2,4,8];
    Append(~sources,C); Append(~source_fs,f); Append(~source_hs,h);
end for;

source_qiso := SameQCurve(sources[1],sources[2]);
source_same_moduli := AbsoluteInvariants(sources[1]) eq
                      AbsoluteInvariants(sources[2]);
print "SOURCE_PAIR_Q_ISOMORPHIC",BoolString(source_qiso);
print "SOURCE_PAIR_SAME_ABSOLUTE_MODULI",BoolString(source_same_moduli);
print "SOURCE_PAIR_QUADRATIC_TWISTS",BoolString(source_same_moduli and not source_qiso);
print "SOURCE_PAIR_Q_ISOGENOUS true certificate each is (2,2)-isogenous to row-51 base; composite degree 16";
assert not source_qiso and not source_same_moduli;

total_kernels := 0;
total_product := 0;
total_jacobians := 0;
total_full := 0;
total_exact := 0;
target_hits := 0;
all_forward_curves := [];
all_forward_source := [];
all_forward_edge := [];

for si in [1..2] do
    C := sources[si];
    triples := OneSplitKernelTriples(C);
    assert #triples eq 3;
    Rs := RichelotIsogenousSurfaces(Jacobian(C));
    assert #Rs eq 3;
    print "FORWARD_SOURCE",si,"reverse_edge",source_edges[si],
          "rational_kernel_count",#triples,"builtin_count",#Rs;

    // First enumerate the three kernels intrinsically and identify the
    // two nondegenerate codomains among Magma's built-in results.
    for k in [1..#triples] do
        A := triples[k][1]; B := triples[k][2]; D := triples[k][3];
        label := triples[k][4];
        total_kernels +:= 1;
        delta := CoefficientMatrixDeterminant(A,B,D);
        print "RATIONAL_KERNEL","source",si,"kernel",k,"label",label,
              "factors",[A,B,D],"delta",delta;
        if delta eq 0 then
            total_product +:= 1;
            print "RATIONAL_KERNEL_CODOMAIN","source",si,"kernel",k,
                  "kind product_of_elliptic_curves";
            continue;
        end if;
        ok,Cm,fm,delta2,g0 := ManualRichelotCurve(A,B,D);
        assert ok and delta2 eq delta;
        matches := [];
        for j in [1..#Rs] do
            if Type(Rs[j]) ne JacHyp then continue; end if;
            Cj,fj,hj := NormalizeCurve(Curve(Rs[j]));
            if SameQCurve(Cm,Cj) then Append(~matches,j); end if;
        end for;
        print "RATIONAL_KERNEL_CODOMAIN","source",si,"kernel",k,
              "kind Jacobian","builtin_matches",matches,
              "manual_polynomial",fm;
        assert #matches eq 1;
    end for;

    // Exact torsion on every Jacobian codomain, including all full-J[2]
    // codomains.  No v_2 filter controls whether this call is made.
    for j in [1..#Rs] do
        if Type(Rs[j]) ne JacHyp then
            print "FORWARD_BUILTIN","source",si,"edge",j,
                  "kind product_of_elliptic_curves";
            continue;
        end if;
        total_jacobians +:= 1;
        D,fD,hD := NormalizeCurve(Curve(Rs[j]));
        G,mp := TorsionSubgroup(Jacobian(D));
        inv := Invariants(G);
        full := FullWeierstrass(D);
        if full then total_full +:= 1; end if;
        total_exact +:= 1;
        print "FORWARD_EXACT","source",si,"edge",j,
              "factor_degrees",FactorDegrees(D),
              "full_weierstrass",BoolString(full),
              "torsion",inv,"order",#G;
        print "FORWARD_POLYNOMIALS",fD,hD;
        Append(~all_forward_curves,D);
        Append(~all_forward_source,si);
        Append(~all_forward_edge,j);
        if inv eq [2,2,4,8] then
            target_hits +:= 1;
            print "TARGET_2248_HIT","source",si,"edge",j,
                  "f",fD,"h",hD;
        end if;
    end for;
end for;

// Q-isomorphism classes among the four Jacobian outputs.
classes := [];
representatives := [];
for i in [1..#all_forward_curves] do
    old := 0;
    for k in [1..#representatives] do
        if SameQCurve(all_forward_curves[i],representatives[k]) then
            old := k; break;
        end if;
    end for;
    if old eq 0 then
        Append(~representatives,all_forward_curves[i]);
        old := #representatives;
    end if;
    Append(~classes,old);
    print "FORWARD_Q_ISOMORPHISM_CLASS","source",all_forward_source[i],
          "edge",all_forward_edge[i],"class",old,
          "is_base",BoolString(SameQCurve(all_forward_curves[i],Cbase));
end for;

print "AUDIT_SUMMARY",
      "raw_sources",#sources,
      "unique_source_Q_curves",2,
      "source_quadratic_twist_pairs",0,
      "geometrically_simple_sources",0,
      "rational_forward_kernels",total_kernels,
      "product_codomains",total_product,
      "Jacobian_codomains",total_jacobians,
      "full_Weierstrass_Jacobians",total_full,
      "exact_torsion_tests",total_exact,
      "unique_forward_Q_curves",#representatives,
      "target_2248_hits",target_hits;

assert total_kernels eq 6;
assert total_product eq 2;
assert total_jacobians eq 4;
assert total_full eq 2;
assert total_exact eq 4;
assert target_hits eq 0;
print "TARGET_2248_SOURCE51_INDEPENDENT_AUDIT_DONE";
quit;
