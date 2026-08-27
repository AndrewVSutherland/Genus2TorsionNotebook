// claude_ov_b4_bfs.m -- Lane 4 (route B4): Galois-stable Richelot BFS from one
// seed of the infinite (6,2) stream of generic [22] curves (Flynn order-11).
//
// One process per seed (sharded); pass IDX (1-based index into b4seeds).
//   magma -b IDX:=7 DEPTH:=3 MemGB:=2 claude_ov_b4_bfs.m
//
// At each node: 2-rank via TwoTorsionSubgroup; if the 2-rank is >= 2 the node
// has torsion containing [2,2] and (by odd-injectivity of a (2,2)-isogeny) the
// rational 11-torsion of the seed, i.e. it is a [2,22] CANDIDATE -> run exact
// TorsionSubgroup and print HIT.  A (2,2)-isogeny also preserves End^0, so
// genericity of the seed transfers to every node.
//
// Markers: SEEDINFO, NODE, HIT, SEEDDONE.

SetColumns(0);
if not assigned IDX   then IDX := 1;   elif Type(IDX) eq MonStgElt then IDX := StringToInteger(IDX); end if;
if not assigned DEPTH then DEPTH := 3; elif Type(DEPTH) eq MonStgElt then DEPTH := StringToInteger(DEPTH); end if;
if not assigned MemGB then MemGB := 2; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
if not assigned SEEDFILE then SEEDFILE := "../data/claude_ov_b4_seeds.m"; end if;
SetMemoryLimit(MemGB*10^9);

load "../data/claude_ov_b4_seeds.m";

Q := Rationals();
P<x> := PolynomialRing(Q);

// integral model: y^2 = f, f = g/D with g in Z[x] -> (D y)^2 = D*g
function IntSextic(f)
    cs := Coefficients(f);
    D := LCM([Denominator(c) : c in cs]);
    g := P![ Numerator(c)*(D div Denominator(c)) : c in cs ];
    h := D*g;
    // strip square content
    ct := GCD([Numerator(c) : c in Coefficients(h)]);
    if ct ne 0 then
        fac := Factorization(ct);
        sq := 1;
        for tt in fac do if tt[2] ge 2 then sq *:= tt[1]^(2*(tt[2] div 2)); end if; end for;
        h := h div sq;
    end if;
    return h;
end function;

function TwoRankOf(J)
    return #Invariants(TwoTorsionSubgroup(J));
end function;

// normalise any genus-2 curve to a small y^2 = f(x) model
function Norm2(C)
    try
        C := ReducedMinimalWeierstrassModel(C);
    catch e
        ;
    end try;
    C := SimplifiedModel(C);
    try
        C := ReducedModel(C);
        C := SimplifiedModel(C);
    catch e
        ;
    end try;
    return C;
end function;

sd := b4seeds[IDX];
n := sd[1]; sg := sd[2]; uu := sd[3]; vv := sd[4]; tt := sd[5];
F := x^6 + 2*x^5 + (2*tt+3)*x^4 + 2*x^3 + (tt^2+1)*x^2 + 2*tt*(1-tt)*x + tt^2;
assert IsDivisibleBy(F, x^2 + uu*x + vv);
f0 := IntSextic(F);
printf "SEEDINFO idx=%o n=%o sg=%o deg=%o coefheight=%o\n", IDX, n, sg, Degree(f0),
        Maximum([#Sprint(Abs(Numerator(c))) : c in Coefficients(f0) | c ne 0]);

C0 := HyperellipticCurve(f0);
J0 := Jacobian(C0);
r0 := TwoRankOf(J0);
printf "SEEDINFO idx=%o seed2rank=%o factortype=%o\n", IDX, r0,
        Sort([Degree(t[1]) : t in Factorization(f0)]);

// quick 11-gate: 11 | #J(F_p) at a few good primes
gate := [];
p := 5;
while #gate lt 4 and p lt 200 do
    if Valuation(Integers()!Discriminant(C0), p) eq 0 then
        Append(~gate, #BaseChange(J0, GF(p)) mod 11);
    end if;
    p := NextPrime(p);
end while;
printf "SEEDINFO idx=%o 11gate=%o\n", IDX, gate;

visited := {};
queue := [* <C0, J0, 0> *];
nnodes := 0; nraise := 0; nchild := 0;

while #queue gt 0 do
    item := queue[1];
    Remove(~queue, 1);
    C := item[1]; J := item[2]; d := item[3];
    inv := G2Invariants(C);
    if inv in visited then continue; end if;
    Include(~visited, inv);
    nnodes +:= 1;
    r := TwoRankOf(J);
    ff := HyperellipticPolynomials(C);
    ft := Sort([Degree(t[1]) : t in Factorization(ff)]);
    printf "NODE idx=%o depth=%o 2rank=%o ftype=%o\n", IDX, d, r, ft;
    if r ge 2 and d ge 1 then
        nraise +:= 1;
        printf "RAISE idx=%o depth=%o 2rank=%o f=%o\n", IDX, d, r, ff;
        tsub := TorsionSubgroup(J);
        printf "HIT idx=%o depth=%o TORSION=%o f=%o\n", IDX, d, Invariants(tsub),
               HyperellipticPolynomials(C);
    end if;
    if d ge DEPTH then continue; end if;
    ok := true;
    RS := [];
    try
        RS := RichelotIsogenousSurfaces(J);
    catch e
        printf "RICHFAIL idx=%o depth=%o\n", IDX, d;
        ok := false;
    end try;
    if not ok then continue; end if;
    printf "RSCOUNT idx=%o depth=%o nRS=%o nJac=%o\n", IDX, d, #RS, #[1 : S in RS | Type(S) eq JacHyp];
    for S in RS do
        if Type(S) ne JacHyp then
            printf "SPLITCODOMAIN idx=%o depth=%o\n", IDX, d;
            continue;
        end if;
        nchild +:= 1;
        CS := Norm2(Curve(S));
        Append(~queue, <CS, Jacobian(CS), d+1>);
    end for;
end while;

printf "SEEDDONE idx=%o n=%o nodes=%o children=%o raises=%o\n", IDX, n, nnodes, nchild, nraise;
quit;
