//////////////////////////////////////////////////////////////////////
// fable_2248_richelot_bfs.m  (2026-07-18, Fable session)
//
// Richelot isogeny-graph BFS hunting [2,2,4,8] (order 128).
//
// Rationale: the all-squares model is closed for d <= 65535 by the
// prod-07 sign-reduction theorem + tier-1 audit; Richelot walks from
// known geometrically simple seeds jump to HIGH-height curves
// structurally, outside the audited box.  Every vertex reached is
// isogenous to a simple seed, hence geometrically simple; any vertex
// with exact torsion [2,2,4,8] is an instant record (order 128).
//
// Seeds (from certify_frontier_existing_banks_20260718 / frontier note):
//   S1: [2,2,4,4]  y^2 = x(x+1296)(x+3249)(x+4096)(x+17424)
//   S2: [2,2,2,8]  y^2 = x(x+1)(x+3025)(x+9801)(x+15625)
// (extend seed list via seeds_extra below as recon supplies more)
//
// Also runs a divisibility jackpot check at every vertex: for each
// order-4 torsion class T4, test IsDivisibleBy(T4, 2) — a rational
// half of a 4-class on a (2,2,4,*) vertex is the 8-chain appearing.
//
// Run:  magma -b code/fable_2248_richelot_bfs.m
// Opts: max_depth:=2  max_vertices:=80
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(2*10^9);

if not assigned max_depth then max_depth := 2;
elif Type(max_depth) eq MonStgElt then max_depth := StringToInteger(max_depth); end if;
if not assigned max_vertices then max_vertices := 80;
elif Type(max_vertices) eq MonStgElt then max_vertices := StringToInteger(max_vertices); end if;

Q := Rationals();
P<x> := PolynomialRing(Q);

function IntegralSquareScale(f)
    den := LCM([Denominator(Coefficient(f,i)) : i in [0..Degree(f)]]);
    return P!(den^2*f), den;
end function;

function NormalizePoly(D)   // curve -> integral y^2 = F model polynomial
    fD, hD := HyperellipticPolynomials(D);
    FD := hD eq 0 select P!fD else P!(hD^2 + 4*fD);
    FI, _ := IntegralSquareScale(FD);
    return FI;
end function;

function G2Key(f)
    C := HyperellipticCurve(f);
    g2 := G2Invariants(C);
    return Sprintf("%o", [Q!a : a in g2]);
end function;

// Torsion + jackpot analysis of one vertex.
procedure AnalyzeVertex(f, tag)
    J := Jacobian(HyperellipticCurve(f));
    T, mp := TorsionSubgroup(J);
    inv := Invariants(T);
    printf "VERTEX %o inv=%o\n", tag, inv;
    if inv eq [2,2,4,8] then
        printf "JACKPOT %o exact [2,2,4,8]! f=%o\n", tag, f;
    end if;
    // divisibility jackpot: try to halve each order-4 class
    n4 := 0; halved := 0;
    for g in T do
        if Order(g) eq 4 then
            n4 +:= 1;
            D := mp(g);
            ok := false;
            try
                ok, half := IsDivisibleBy(D, 2);
            catch ee
                ok := false;
            end try;
            if ok then
                halved +:= 1;
                printf "HALF4 %o a 4-class halves rationally: inv=%o\n", tag, inv;
            end if;
        end if;
    end for;
    if n4 gt 0 then
        printf "DIV4 %o order4classes=%o halved=%o\n", tag, n4, halved;
    end if;
end procedure;

// seeds
seeds := [
    x*(x+1296)*(x+3249)*(x+4096)*(x+17424),
    x*(x+1)*(x+3025)*(x+9801)*(x+15625)
];
seed_tags := ["S1_2244", "S2_2228"];
// extra seeds may be appended by editing here:
seeds_extra := [];
tags_extra := [];
seeds cat:= seeds_extra;
seed_tags cat:= tags_extra;

visited := {};          // G2 keys
queue := [];            // <poly, depth, tag>

for i in [1..#seeds] do
    fI, _ := IntegralSquareScale(seeds[i]);
    k := G2Key(fI);
    if k notin visited then
        Include(~visited, k);
        Append(~queue, <fI, 0, seed_tags[i]>);
    end if;
end for;

cursor := 1;
while cursor le #queue and #queue le max_vertices do
    node := queue[cursor];
    f := node[1]; depth := node[2]; tag := node[3];
    AnalyzeVertex(f, tag);
    if depth lt max_depth then
        J := Jacobian(HyperellipticCurve(f));
        Rs := [];
        try
            Rs := RichelotIsogenousSurfaces(J);
        catch ee
            printf "RICHELOT_FAIL %o\n", tag;
        end try;
        for i in [1..#Rs] do
            S := Rs[i];
            tS := Type(S);
            if tS eq JacHyp or tS eq CrvHyp then
                D := tS eq JacHyp select Curve(S) else S;
                fN := NormalizePoly(D);
                k := G2Key(fN);
                if k notin visited then
                    Include(~visited, k);
                    Append(~queue, <fN, depth+1, Sprintf("%o.%o", tag, i)>);
                end if;
            else
                printf "NONJAC %o.%o type=%o\n", tag, i, tS;
            end if;
        end for;
    end if;
    cursor +:= 1;
end while;

printf "BFS_DONE vertices=%o (max_depth=%o)\n", #queue, max_depth;
quit;
