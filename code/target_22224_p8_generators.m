SetColumns(0); Q:=Rationals();
curves:=[
    EllipticCurve([Q|0,0,0,1520244,202972176]),
    EllipticCurve([Q|0,0,0,1837869,458936786]),
    EllipticCurve([Q|0,0,0,-13836108,19808836432])
];
for i in [1..3] do
    E:=curves[i]; print "P8_GENERATORS",i,Generators(E),"height_pairing",
        HeightPairingMatrix(Generators(E));
end for;
quit;
