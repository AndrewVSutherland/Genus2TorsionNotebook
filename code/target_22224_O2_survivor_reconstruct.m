//////////////////////////////////////////////////////////////////////
// Exact reconstruction of the two modular survivors of the O2
// million-bound rank-one Mordell--Weil sieve.
//////////////////////////////////////////////////////////////////////
Q:=Rationals(); R<T>:=PolynomialRing(Q);
SetLogFile("results/target_22224_O2_rank1_survivor_reconstruct.log":Overwrite:=true);
fixed:=[Q!-1,Q!6,Q!9]; d0:=Q!-507/98;
Rx:=(fixed[1]+d0*T^2)/(fixed[1]+d0);
Ry:=(fixed[3]+d0*T^2)/(fixed[3]+d0);
C:=HyperellipticCurve(Rx*Ry); P0:=C![Q!1,Q!1,Q!1];
Eraw,phi:=EllipticCurve(C,P0); Einv:=Inverse(phi);
E,minmap:=MinimalModel(Eraw); minmapinv:=Inverse(minmap);
gens:=Generators(E); G:=[P:P in gens|Order(P) eq 0][1];
TG,tmap:=TorsionSubgroup(E); tors:=[tmap(g):g in TG];

print "O2_SURVIVOR_RECONSTRUCT_START",E,"G",G,"torsion_cosets",tors;
for m in [-3,0] do
    P:=m*G+tors[1]; cp:=Einv(minmapinv(P));
    if cp[3] eq 0 then
        print "O2_SURVIVOR",m,1,"t", "infinity";
    else
        t:=Q!(cp[1]/cp[3]);
        print "O2_SURVIVOR",m,1,"t",t,
              "v13_den",Valuation(Denominator(t),13),
              "v11_den",Valuation(Denominator(t),11);
    end if;
end for;
UnsetLogFile();
quit;
