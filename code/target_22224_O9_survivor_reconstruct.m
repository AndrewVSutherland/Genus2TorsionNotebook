//////////////////////////////////////////////////////////////////////
// Exact reconstruction of the two modular survivors of the O9
// million-by-million rank-two Mordell--Weil sieve.
//////////////////////////////////////////////////////////////////////
SetSeed(8); Q:=Rationals(); R<T>:=PolynomialRing(Q);
SetLogFile("results/target_22224_O9_12_rank2_survivor_reconstruct.log":Overwrite:=true);
fixed:=[Q!1,Q!55,Q!99]; d0:=Q!5/9;
Rx:=(fixed[1]+d0*T^2)/(fixed[1]+d0);
Ry:=(fixed[2]+d0*T^2)/(fixed[2]+d0);
C:=HyperellipticCurve(Rx*Ry); P0:=C![Q!1,Q!1,Q!1];
Eraw,phi:=EllipticCurve(C,P0); Einv:=Inverse(phi);
E,minmap:=MinimalModel(Eraw); minmapinv:=Inverse(minmap);
gens:=Generators(E); free:=[P:P in gens|Order(P) eq 0];
TG,tmap:=TorsionSubgroup(E); tors:=[tmap(g):g in TG];

print "O9_SURVIVOR_RECONSTRUCT_START",E,"free",free,"torsion_cosets",tors;
for s in [<1,-1,1>,<0,0,1>] do
    m:=s[1]; n:=s[2]; ti:=s[3];
    P:=m*free[1]+n*free[2]+tors[ti]; cp:=Einv(minmapinv(P));
    if cp[3] eq 0 then
        print "O9_SURVIVOR",m,n,ti,"t","infinity";
    else
        t:=Q!(cp[1]/cp[3]); dd:=d0*t^2;
        squares:=[IsSquare((z+dd)/(z+d0)):z in fixed];
        print "O9_SURVIVOR",m,n,ti,"t",t,"full_square_flags",squares;
    end if;
end for;
UnsetLogFile(); quit;
