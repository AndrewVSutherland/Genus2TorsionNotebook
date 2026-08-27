//////////////////////////////////////////////////////////////////////
// Exact verifier for a nondegenerate point on the A,B,C square surface.
//
// Default: the third orbit found by search_22212_abc.cpp.  Its smallest
// direct-search chart is (4165,19661,5364); the default symmetry chart
// (16660,78644,78793) gives a smaller primitive (s,m,n).
//
// Run:
//   magma -b code/verify_22212_abc_hit.m
// or:
//   magma -b AA:=4165 BB:=19661 CC:=5364 code/verify_22212_abc_hit.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);

if not assigned AA then AA := 16660;
elif Type(AA) eq MonStgElt then AA := StringToInteger(AA);
end if;
if not assigned BB then BB := 78644;
elif Type(BB) eq MonStgElt then BB := StringToInteger(BB);
end if;
if not assigned CC then CC := 78793;
elif Type(CC) eq MonStgElt then CC := StringToInteger(CC);
end if;

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);
A := Z!AA; B := Z!BB; Cc := Z!CC;
X := A^2; Y := B^2; Zc := Cc^2;
D := X-Y+Zc;
L := (X-Y)*(Y-Zc);
F := D*L;
G := D*(L+D^2);
assert D ne 0 and IsSquare(F) and IsSquare(G);
yroot := Isqrt(F); zroot := Isqrt(G);
assert yroot^2 eq F and zroot^2 eq G and G-F eq D^3;

delta_factors := [
    Cc, B, B-Cc, B+Cc, -A+B, A, A-Cc, A+Cc, A+B,
    -X+Y-2*Zc,
    -X+Y-Zc,
    -X+Y-A*Cc-Zc,
    -X+Y+A*Cc-Zc,
    -X*Y+Y^2+2*X*Zc-2*Y*Zc+Zc^2,
    -X^2+X*Y-X*Zc+Y*Zc-Zc^2,
    X^2-2*X*Y+Y^2+2*X*Zc-Y*Zc
];
assert &and[v ne 0 : v in delta_factors];

s := Q!(A*(X-Y+2*Zc))/(2*B*D);
m := Q!(-2*X*Y+2*Y^2+2*X*Zc-6*Y*Zc+4*Zc^2)/(2*A*B*D);
n := Q!A/B;
den := LCM([Denominator(v) : v in [s,m,n]]);
si := Z!(den*s); mi := Z!(den*m); ni := Z!(den*n);
gg := GCD([Abs(si),Abs(mi),Abs(ni)]);
si div:= gg; mi div:= gg; ni div:= gg;
if ni lt 0 then si *:= -1; mi *:= -1; ni *:= -1; end if;

function CurveFromSMN(ss,mm,nn)
    avec := [1,1,1,2,2];
    bvec := [
        2*ss^2-ss*nn,
        2*ss^2+ss*mm-2*ss*nn-mm*nn,
        2*ss^2+ss*mm-ss*nn-mm*nn,
        -mm*nn,
        4*ss^2-4*ss*nn-mm*nn
    ];
    ff := &*[avec[i]+bvec[i]*x : i in [1..5]];
    assert Degree(ff) eq 5 and Discriminant(ff) ne 0;
    return HyperellipticCurve(ff), ff, bvec;
end function;

C0, f0, bvec := CurveFromSMN(si,mi,ni);
printf "ABC = (%o,%o,%o)\n", A,B,Cc;
printf "D = %o, y = %o, z = %o\n", D,yroot,zroot;
printf "(s,m,n) = (%o,%o,%o)\n", si,mi,ni;
printf "B-vector = %o\n", bvec;

// Simplifying before the arithmetic avoids carrying the discovery model's
// very large branch coefficients through every subsequent computation.
Cs := SimplifiedModel(C0);
fs, hs := HyperellipticPolynomials(Cs);
assert Discriminant(fs) ne 0;
printf "SIMPLIFIED_CURVE = %o\n", Cs;
printf "SIMPLIFIED_POLYNOMIALS = %o, %o\n", fs,hs;

Cmin := ReducedMinimalWeierstrassModel(C0);
fmin, hmin := HyperellipticPolynomials(Cmin);
assert Discriminant(Cmin) ne 0;
printf "REDUCED_MINIMAL_CURVE = %o\n", Cmin;
printf "REDUCED_MINIMAL_POLYNOMIALS = %o, %o\n", fmin,hmin;
mindisc := Z!Discriminant(Cmin);
printf "REDUCED_MINIMAL_DISCRIMINANT = %o\n", mindisc;
printf "REDUCED_MINIMAL_DISCRIMINANT_FACTORIZATION = %o\n", Factorization(mindisc);

// Magma's rational torsion routine currently requires an integral y^2=f
// model, so use Cs here; Cmin is retained for the small displayed equation.
J := Jacobian(Cs);
T := TorsionSubgroup(J);
invs := Invariants(T);
printf "TORSION_INVARIANTS = %o\n", invs;
printf "TORSION_ORDER = %o\n", #T;
assert invs eq [2,2,2,12] and #T eq 96;

// Strict Frobenius root-power witnesses for geometric simplicity.  Requiring
// degree four for pi^e, 2<=e<=12, rules out the usual power-induced splitting
// pathologies, not merely rational decomposability.
witnesses := [];
disc := Z!Discriminant(Cmin);
for pp in PrimesInInterval(17,500) do
    if disc mod pp eq 0 then continue; end if;
    Cp := ChangeRing(Cmin,GF(pp));
    lp := LPolynomial(Cp);
    chi := P!Reverse(Coefficients(lp));
    if not IsIrreducible(chi) or Degree(chi) ne 4 then continue; end if;
    K<pi> := NumberField(chi);
    degrees := [Degree(MinimalPolynomial(pi^e)) : e in [2..12]];
    if &and[d eq 4 : d in degrees] then
        Append(~witnesses,<pp,chi>);
        printf "SIMPLICITY_WITNESS p=%o chi=%o powers_2_to_12=%o\n",
               pp,chi,degrees;
        if #witnesses ge 4 then break; end if;
    end if;
end for;
assert #witnesses ge 2;

g2 := G2Invariants(Cmin);
orbit_smn := [
    [1500518600,253638081,138777800],
    [1591920768,162235913,230179968],
    [-2493301632,5240700751,460359936],
    [-2584703800,5606309423,277555600],
    [3254675281,-507276162,6208336849],
    [3346077449,-324471826,6208336849],
    [254097487,-10481401502,555111200],
    [436901823,-11212618846,920719872],
    [944177985,22342964000,920719872],
    [12115659985,-22342964000,12416673698],
    [578569313,22708572672,555111200],
    [11932855649,-22708572672,12416673698]
];
same_orbit := true;
for v in orbit_smn do
    Corbit, _, _ := CurveFromSMN(v[1],v[2],v[3]);
    same_orbit and:= G2Invariants(Corbit) eq g2;
end for;
printf "ALL_12_CHARTS_SAME_G2 = %o\n", same_orbit;
assert same_orbit;

Cknown1, _, _ := CurveFromSMN(336396,-689185,-166464);
Cknown2, _, _ := CurveFromSMN(2208,-8303,-7200);
printf "G2_INVARIANTS = %o\n", g2;
printf "DISTINCT_FROM_CURVE1 = %o\n", g2 ne G2Invariants(Cknown1);
printf "DISTINCT_FROM_CURVE2 = %o\n", g2 ne G2Invariants(Cknown2);
assert g2 ne G2Invariants(Cknown1) and g2 ne G2Invariants(Cknown2);

printf "VERIFIED_NEW_22212_ORBIT\n";
quit;
