// [2,22] hunt, route 3: Richelot/2-isogeny neighborhood of the known
// NON-RM cyclic-[22] curve LMFDB 1192.a.19072.1
//   y^2 + (x^3+x)y = x^3 - 2x^2 - x + 1   (End(J_Qbar) = Z, geom. simple).
// (2,2)-isogenies preserve the 11-part of torsion and End tensor Q, but can
// change the rational 2-torsion structure: any vertex of the graph with
// 2-rank 2 is a NON-RM [2,22] Jacobian.  BFS via RichelotIsogenousSurfaces
// up to depth 3, exact TorsionSubgroup at every genus-2 vertex.
// Also runs the same BFS from the Flynn-family [22] members t=-36, t=-static
// second seed for coverage.
//
// Run: magma -b code/z22_richelot_bfs.m > results/z22_richelot_bfs.log

SetColumns(0);
SetSeed(1);
SetMemoryLimit(3*10^9);

Q := Rationals();
P<x> := PolynomialRing(Q);

function TorsOf(C)
    Cc := C;
    try
        Cm := ReducedMinimalWeierstrassModel(Cc);
        Cc := SimplifiedModel(Cm);
    catch e ; end try;
    return Invariants(TorsionSubgroup(Jacobian(Cc))), Cc;
end function;

procedure BFS(C0, name, depth)
    printf "SEED %o\n", name;
    inv0, Cs0 := TorsOf(C0);
    printf "SEED_TORSION %o = %o\n", name, inv0;
    seen := { G2Invariants(Cs0) };
    frontier := [ Cs0 ];
    for level in [1..depth] do
        newfrontier := [];
        for C in frontier do
            surfs := [];
            try
                surfs := RichelotIsogenousSurfaces(Jacobian(C));
            catch e
                printf "LEVEL %o: RichelotIsogenousSurfaces error on a vertex\n", level;
                continue;
            end try;
            for S in surfs do
                // only genus-2 Jacobian vertices matter (products = split)
                if Type(S) eq SetCart then
                    printf "LEVEL %o SPLIT_VERTEX (product of elliptic curves)\n", level;
                    continue;
                end if;
                Cn := Curve(S);
                g2 := G2Invariants(Cn);
                if g2 in seen then continue; end if;
                Include(~seen, g2);
                invn, Csn := TorsOf(Cn);
                printf "LEVEL %o VERTEX torsion=%o\n", level, invn;
                if #invn ge 2 and invn[#invn] mod 11 eq 0 and invn[#invn-1] mod 2 eq 0 then
                    printf "HIT22 %o level=%o torsion=%o curve=%o\n",
                           name, level, invn, HyperellipticPolynomials(Csn);
                end if;
                Append(~newfrontier, Csn);
            end for;
        end for;
        printf "LEVEL %o DONE new_vertices=%o\n", level, #newfrontier;
        frontier := newfrontier;
        if #frontier eq 0 then break; end if;
    end for;
    printf "BFS_DONE %o total_vertices=%o\n", name, #seen;
end procedure;

// Seed 1: LMFDB 1192.a.19072.1 (End = Z, torsion [22])
C1 := HyperellipticCurve(x^3 - 2*x^2 - x + 1, x^3 + x);
BFS(C1, "1192a", 3);

// Seed 2: Flynn-family [22] member t = -36 (order-22 class; End generic)
F36 := x^6 + 2*x^5 + (2*(-36)+3)*x^4 + 2*x^3 + ((-36)^2+1)*x^2 + 2*(-36)*(1-(-36))*x + (-36)^2;
C2 := HyperellipticCurve(F36);
BFS(C2, "flynn_t-36", 3);

printf "Z22_RICHELOT_BFS_DONE\n";
quit;
