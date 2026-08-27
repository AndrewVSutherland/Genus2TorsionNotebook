//////////////////////////////////////////////////////////////////////
// Complete rational Richelot components of the two known simple Z/24
// Jacobians.  The 3-primary part is preserved along every edge; we look
// for a vertex where the 2-primary part becomes [2,2,2,8].
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetLogFile("results/target_22224_z24_richelot_bfs.log" : Overwrite := true);
Q := Rationals();
P<x> := PolynomialRing(Q);

function A8f(rv,pv,tv)
    e := tv^2-2*pv*tv/rv;
    d := e+2*pv-rv^2;
    lambda := rv/tv;
    aa := rv^2-lambda;
    bb := 2*rv*pv-2*lambda*(pv+rv*tv)+2*rv*lambda;
    cc := pv^2+2*pv*rv^2-rv^4-rv^3*tv-rv*pv^2/tv
          -lambda*(rv^2+e)
          +2*lambda*(rv*pv+rv^2*tv-3*pv*tv+rv*tv^2);
    q := aa*x^2+bb*x+cc;
    Qq := x^2+d;
    return q*(Qq^2+q);
end function;

function Normalize(C)
    f,h := HyperellipticPolynomials(C);
    F := h eq 0 select P!f else P!(h^2+4*f);
    den := LCM([Denominator(Coefficient(F,i)) : i in [0..Degree(F)]]);
    return SimplifiedModel(HyperellipticCurve(P!(den^2*F)));
end function;

procedure Explore(name,f)
    C0 := Normalize(HyperellipticCurve(f));
    curves := [C0];
    invs := [[Q!z : z in G2Invariants(C0)]];
    depth := [0];
    cursor := 1;
    edges := 0;
    print "Z24_COMPONENT_START",name;
    while cursor le #curves do
        C := curves[cursor];
        T := TorsionSubgroup(Jacobian(C));
        print "Z24_NODE",name,cursor,"depth",depth[cursor],
              "torsion",Invariants(T);
        if Invariants(T) eq [2,2,2,24] then
            ff,hh := HyperellipticPolynomials(C);
            print "TARGET_22224_HIT",name,cursor,ff,hh;
        end if;
        if depth[cursor] ge 8 then
            cursor +:= 1;
            continue;
        end if;
        Rs := RichelotIsogenousSurfaces(Jacobian(C));
        print "Z24_EDGES",name,cursor,#Rs;
        for j in [1..#Rs] do
            edges +:= 1;
            if Type(Rs[j]) ne JacHyp then
                print "Z24_NONJAC",name,cursor,j,Type(Rs[j]);
                continue;
            end if;
            D := Normalize(Curve(Rs[j]));
            gi := [Q!z : z in G2Invariants(D)];
            old := Index(invs,gi);
            if old eq 0 then
                Append(~curves,D); Append(~invs,gi);
                Append(~depth,depth[cursor]+1);
                print "Z24_NEW",name,#curves,"parent",cursor,"edge",j;
            else
                print "Z24_OLD",name,cursor,j,"node",old;
            end if;
        end for;
        cursor +:= 1;
    end while;
    print "Z24_COMPONENT_DONE",name,"nodes",#curves,"edges",edges;
end procedure;

print "TARGET_22224_Z24_RICHELOT_BFS_START";
Explore("A",A8f(5,-5/2,-9/2));
Explore("B",A8f(1/3,-1/9,-1));
print "TARGET_22224_Z24_RICHELOT_BFS_DONE";
quit;
