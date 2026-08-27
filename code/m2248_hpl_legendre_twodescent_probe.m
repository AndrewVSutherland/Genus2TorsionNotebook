//////////////////////////////////////////////////////////////////////
//  Probe 2-descent covers after removing the HPL point on the
//  fixed-rho HPL Legendre model.
//////////////////////////////////////////////////////////////////////

Q := Rationals();

rho0 := Q!58466134224 / Q!53109477625;
tau0 := Q!307598400 / Q!352612321;
dB := Q!72946054224 / Q!53109477625;
Y0 := dB*(tau0^2-1);
V0 := tau0*Y0;

E0 := EllipticCurve([Q!0, -(Q!1+rho0^2), Q!0, rho0^2, Q!0]);
P0 := E0![tau0^2, V0, 1];

print "E0", E0;
print "P0", P0;
print "RankBounds", RankBounds(E0);
print "Starting TwoDescent with torsion and P0 removed";
SetVerbose("TwoDescent", 1);
covers, maps := TwoDescent(E0 : RemoveTorsion := true, RemoveGens := {P0});
print "num_covers", #covers;
for i in [1..#covers] do
    print "COVER", i, covers[i];
    print "MAP", i, maps[i];
    print "DOMAIN", Domain(maps[i]);
    print "CODOMAIN", Codomain(maps[i]);
    try
        pts := Points(covers[i] : Bound := 1000);
        print "POINTS_BOUND_1000", i, #pts;
        for j in [1..Min(#pts, 5)] do
            print "cover_point", i, j, pts[j];
            try
                print "image", maps[i](pts[j]);
            catch e
                print "image_failed", e`Object;
            end try;
        end for;
    catch e
        print "points_failed", i, e`Object;
    end try;
end for;
