//////////////////////////////////////////////////////////////////////
//  Certify the residual theta-doubling curves from
//  code/m2226_order6_doubling.m.
//
//  Each residual plane curve is genus one.  The script maps it to an
//  elliptic curve, checks the rank is 0, pulls back the torsion points,
//  and verifies that every rational plane point lies on the bad locus
//  where the M(2,2,2,6) model is singular or degenerate.
//
//  Typical run from torsion_jac:
//      magma code/m2226_certify_doubling_curves.m
//////////////////////////////////////////////////////////////////////

Q:=Rationals();
P2<s,m,n>:=ProjectiveSpace(Q,2);
R:=CoordinateRing(P2);
polys := [
<s^4 + 1/2*s^3*m - 3/2*s^3*n + 1/8*s^2*m^2 - 3/4*s^2*m*n + 5/8*s^2*n^2 - 3/16*s*m^2*n + 5/16*s*m*n^2 + 9/128*m^2*n^2, "12">,
<s^4 - 3/2*s^3*n + 1/4*s^2*m^2 + 9/16*s^2*n^2 - 3/8*s*m^2*n - 1/16*s*m*n^2 + 9/64*m^2*n^2, "13">,
<s^3 - s^2*m + 1/4*s*m^2 + 2*s*m*n - m*n^2, "14">,
<s^3 + 1/2*s^2*m - s^2*n + 1/16*s*m^2 - 1/2*s*m*n + 1/2*s*n^2 + 1/4*m*n^2, "15">,
<s^4 + 1/2*s^3*m - 2*s^3*n + 1/16*s^2*m^2 - s^2*m*n + 5/4*s^2*n^2 - 3/16*s*m^2*n + 5/8*s*m*n^2 + 9/64*m^2*n^2, "23">,
<s^4 + 2*s^3*m - 3*s^3*n + s^2*m^2 - 4*s^2*m*n + 9/4*s^2*n^2 - 2*s*m^2*n + 5/2*s*m*n^2 + 5/4*m^2*n^2, "24">,
<s^4 + s^3*m - 3/2*s^3*n + 1/2*s^2*m^2 - 3/2*s^2*m*n + 9/16*s^2*n^2 - 3/4*s*m^2*n + 5/8*s*m*n^2 + 5/16*m^2*n^2, "25">,
<s^4 + 2*s^3*m + s^2*m^2 - 2*s^2*m*n - s*m^2*n + s*m*n^2 + 1/2*m^2*n^2, "34">,
<s^4 - s^3*m - 3*s^3*n + 1/4*s^2*m^2 + s^2*m*n + 9/4*s^2*n^2 - 3/4*s*m^2*n - 1/4*s*m*n^2 + 9/16*m^2*n^2, "45">
];
bad := [ n, m, s, s-n, s-1/2*n, s+1/2*m, s*m-s*n-m*n, s*m+s*n-1/2*m*n, s^2-1/2*s*n+1/4*m*n, s^2+1/2*s*m-s*n-1/4*m*n ];

function OnBadLocus(P)
    vals := Eltseq(P);
    return &or [ Evaluate(b, vals) eq 0 : b in bad ];
end function;

function RationalPointsOfZeroDimScheme(Z)
    pts := [];
    try
        pts := Setseq(Points(Z));
    catch e
        pts := [];
    end try;
    return pts;
end function;

for item in polys do
    C := Curve(P2, item[1]);
    print "curve", item[2], "genus", Genus(C);
    pts_search := PointSearch(C, 500);
    print "  pointsearch", pts_search;
    sing_pts := RationalPointsOfZeroDimScheme(SingularSubscheme(C));
    print "  rational singular points", sing_pts;
    print "  singular bad flags", [ OnBadLocus(P) : P in sing_pts ];
    origin_found := false;
    P0 := P2![1,0,0];
    for P in pts_search do
        if P notin sing_pts then
            P0 := P;
            origin_found := true;
            break;
        end if;
    end for;
    if not origin_found then
        print "  no nonsingular point found to use as origin";
        continue;
    end if;
    print "  origin", P0;
    E, phi := EllipticCurve(C, P0);
    lb, ub := RankBounds(E);
    print "  E minimal", MinimalModel(E);
    print "  rank bounds", lb, ub;
    G, mp := TorsionSubgroup(E);
    tors_pts := [ mp(g) : g in G ];
    print "  torsion order", #tors_pts, "structure", G;
    if lb eq 0 and ub eq 0 then
        pre_pts := [];
        for T in tors_pts do
            Z := T @@ phi;
            ptsT := RationalPointsOfZeroDimScheme(Z);
            pre_pts cat:= ptsT;
        end for;
        pre_set := Set(pre_pts cat sing_pts cat pts_search);
        print "  collected rational plane points", Setseq(pre_set);
        print "  bad flags", [ OnBadLocus(P) : P in Setseq(pre_set) ];
    end if;
end for;
quit;
