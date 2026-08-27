// The genus-1 nondegenerate quartic component of S cap {x2 = x3 + x5}:
// ideal, rational points, elliptic model, rank.
SetColumns(0);
SetMemoryLimit(8*10^9);
Q := Rationals();
P3<y1,y2,y3,y4> := ProjectiveSpace(Q, 3);   // (x1, x3, x4, x5)
X := [ y1, y2+y4, y2, y3, y4 ];             // x2 = x3+x5
Q2 := X[1]^2+X[2]^2+X[3]^2-X[4]^2-X[5]^2;
F4 := X[1]^4+X[2]^4+X[3]^4-X[4]^4-X[5]^4;
Se := Scheme(P3, [Q2, F4]);
comps := IrreducibleComponents(Se);
for co in comps do
    cor := ReducedSubscheme(co);
    if Dimension(cor) ne 1 or Degree(cor) ne 4 then continue; end if;
    printf "quartic component ideal:\n";
    B := MinimalBasis(Ideal(cor));
    for b in B do printf "  %o\n", b; end for;
    Cu := Curve(cor);
    printf "genus %o\n", Genus(Cu);
    // rational points by small search
    pts := PointSearch(Cu, 500);
    printf "small points (h<=500): %o\n", #pts;
    ndeg := [];
    for pt in pts do
        xx := [ Evaluate(f, Eltseq(pt)) : f in [y1, y2+y4, y2, y3, y4] ];
        sq := [ t^2 : t in xx ];
        if &and[ t ne 0 : t in xx ] and #Set(sq) eq 5 then Append(~ndeg, xx); end if;
        printf "  pt %o -> x = %o%o\n", Eltseq(pt), xx,
            (&and[ t ne 0 : t in xx ] and #Set(sq) eq 5) select "  NONDEGENERATE" else "";
    end for;
    printf "nondegenerate points found: %o\n", #ndeg;
    if #pts gt 0 then
        E, mp := EllipticCurve(Cu, pts[1]);
        Em, _ := MinimalModel(E);
        printf "elliptic Jacobian: %o\n", aInvariants(Em);
        printf "torsion: %o\n", Invariants(TorsionSubgroup(Em));
        ok, r := RankBounds(Em);
        printf "rank bounds: %o %o\n", ok, r;
    end if;
end for;
printf "ELLSEC_DONE\n";
quit;
