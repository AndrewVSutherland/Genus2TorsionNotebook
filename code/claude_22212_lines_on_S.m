// Do any lines through pairs of nodes lie on S: x1^2+x2^2+x3^2 = x4^2+x5^2,
// x1^4+...= x4^4+x5^4?  A line in S with generic point in S^o would give
// INFINITELY many [2,2,2,12] fibers.
SetColumns(0);
Q := Rationals();
// nodes: orbit of (1,1,0,0,0) (12) and (1,1,0,1,1) (24), per the paper.
nodes := [];
// orbit 1: two of x,y,z are +-1, rest 0 with the quadric/quartic satisfied:
// (1,\pm1,0,0,0)-type: check membership generically below; enumerate all
// vectors with entries in {-1,0,1} and test the two equations + singularity-free test skipped.
cands := CartesianPower([-1,0,1], 5);
for c in cands do
    v := [c[i] : i in [1..5]];
    if &and[x eq 0 : x in v] then continue; end if;
    q1 := v[1]^2+v[2]^2+v[3]^2 - v[4]^2 - v[5]^2;
    q2 := v[1]^4+v[2]^4+v[3]^4 - v[4]^4 - v[5]^4;
    if q1 eq 0 and q2 eq 0 then Append(~nodes, v); end if;
end for;
// dedupe projectively (v ~ -v)
uniq := [];
seen := {};
for v in nodes do
    w1 := v; w2 := [-x : x in v];
    if w1 in seen or w2 in seen then continue; end if;
    Include(~seen, w1);
    Append(~uniq, v);
end for;
printf "candidate node-like points (proj): %o\n", #uniq;
P2<u,w> := PolynomialRing(Q, 2);
nlines := 0;
for i in [1..#uniq] do
    for j in [i+1..#uniq] do
        a := uniq[i]; b := uniq[j];
        L := [ u*a[k] + w*b[k] : k in [1..5] ];
        Q1 := L[1]^2+L[2]^2+L[3]^2-L[4]^2-L[5]^2;
        Q2 := L[1]^4+L[2]^4+L[3]^4-L[4]^4-L[5]^4;
        if Q1 eq 0 and Q2 eq 0 then
            nlines +:= 1;
            // generic point nondegenerate? need all coords nonzero and squares distinct
            gen := [ 2*a[k] + 3*b[k] : k in [1..5] ];   // u=2,w=3
            allnz := &and[ x ne 0 : x in gen ];
            sq := [ x^2 : x in gen ];
            dist := #Set(sq) eq 5;
            printf "LINE through %o and %o : generic nonzero=%o distinct-squares=%o gen=%o\n",
                a, b, allnz, dist, gen;
        end if;
    end for;
end for;
printf "lines found: %o\n", nlines;
printf "LINES_DONE\n";
quit;
