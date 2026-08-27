//////////////////////////////////////////////////////////////////////
// Geometry probes for the R = -25/4 S_B pullback cover.
//
// Historical filename note: an earlier version missed the 609/256 scale
// factor in the reduced fiber coordinate and produced a spurious genus-5
// model.  This file now uses the corrected genus-3 model.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();

if not assigned PointBound then PointBound := 100; end if;
if Type(PointBound) eq MonStgElt then PointBound := StringToInteger(PointBound); end if;
if not assigned DoAuto then DoAuto := false; end if;
if Type(DoAuto) eq MonStgElt then DoAuto := DoAuto in {"true", "True", "1", "yes"}; end if;
if not assigned DoSubgroups then DoSubgroups := false; end if;
if Type(DoSubgroups) eq MonStgElt then DoSubgroups := DoSubgroups in {"true", "True", "1", "yes"}; end if;
if not assigned DoRank then DoRank := false; end if;
if Type(DoRank) eq MonStgElt then DoRank := DoRank in {"true", "True", "1", "yes"}; end if;

A3<m,g,Y> := AffineSpace(Q, 3);

// g is the reduced quartic coordinate on g^2 = Gamma(m).
eq1 := -1024*m^4 + 865600*m^2 + g^2 - 231800625;

// Y=m*y, with the corrected 609/256 scale factor included.
eq2 := 51200*Y^2 - 29696*m^4 - 928*m^2*g
       + 27051200*m^2 + 441525*g - 6722218125;

Caff := Curve(A3, [eq1, eq2]);
printf "affine dimension=%o\n", Dimension(Caff);
Cp := ProjectiveClosure(Caff);
printf "projective ambient=%o degree=%o dimension=%o\n", Ambient(Cp), Degree(Cp), Dimension(Cp);
printf "normalized genus=%o\n", Genus(Cp);

try
    printf "is nonsingular: %o\n", IsNonsingular(Cp);
catch e
    printf "IsNonsingular failed: %o\n", e`Object;
end try;

try
    sing := SingularPoints(Cp);
    printf "singular points count=%o\n", #sing;
    if #sing le 20 then printf "singular points=%o\n", sing; end if;
catch e
    printf "SingularPoints failed: %o\n", e`Object;
end try;

try
    pts := Points(Cp : Bound := PointBound);
    printf "projective rational points Bound %o count=%o\n", PointBound, #pts;
    if #pts le 50 then printf "points=%o\n", pts; end if;
catch e
    printf "Points(Cp) failed: %o\n", e`Object;
end try;

try
    ishyp, hypmap := IsHyperelliptic(Cp);
    printf "IsHyperelliptic=%o\n", ishyp;
    if ishyp then
        printf "hyperelliptic map=%o\n", hypmap;
    end if;
catch e
    printf "IsHyperelliptic failed: %o\n", e`Object;
end try;

try
    can := CanonicalMap(Cp);
    printf "CanonicalMap succeeded: domain genus=%o codomain=%o\n", Genus(Domain(can)), Codomain(can);
catch e
    printf "CanonicalMap failed: %o\n", e`Object;
end try;

if Genus(Cp) eq 5 then
    try
        gon, mp, X4, f4 := Genus5GonalMap(Cp);
        printf "Genus5GonalMap succeeded: gonality=%o\n", gon;
        printf "map=%o\n", mp;
        printf "X4=%o\n", X4;
        printf "f4=%o\n", f4;
    catch e
        printf "Genus5GonalMap failed: %o\n", e`Object;
    end try;
else
    printf "Genus5GonalMap skipped: normalized genus is %o\n", Genus(Cp);
end if;

// Explicit quotients by visible involutions.
// 1.  Y -> -Y gives the elliptic fiber, already known.
// 2.  m -> -m is represented by X=m^2 and is the best chance for a
//     hyperelliptic lower-genus auxiliary curve.
// 3.  both signs gives the singular/easy quotient checked in the
//     descent probe.
A3q<Xq,gq,Yq> := AffineSpace(Q, 3);
eq1_msign := -1024*Xq^2 + 865600*Xq + gq^2 - 231800625;
eq2_msign := 51200*Yq^2 - 29696*Xq^2 - 928*Xq*gq
             + 27051200*Xq + 441525*gq - 6722218125;
Qm := ProjectiveClosure(Curve(A3q, [eq1_msign, eq2_msign]));
printf "quotient by m -> -m: degree=%o genus=%o nonsingular=%o\n",
    Degree(Qm), Genus(Qm), IsNonsingular(Qm);
try
    ishypQm, hypQm := IsHyperelliptic(Qm);
    printf "quotient m-sign IsHyperelliptic=%o\n", ishypQm;
    if ishypQm then printf "quotient m-sign hyperelliptic map=%o\n", hypQm; end if;
catch e
    printf "quotient m-sign IsHyperelliptic failed: %o\n", e`Object;
end try;
try
    ptsQm := Points(Qm : Bound := PointBound);
    printf "quotient m-sign points Bound %o count=%o\n", PointBound, #ptsQm;
    if #ptsQm le 50 then printf "quotient m-sign points=%o\n", ptsQm; end if;
catch e
    printf "quotient m-sign Points failed: %o\n", e`Object;
end try;
try
    G1m := GenusOneModel(Qm);
    printf "quotient m-sign GenusOneModel degree=%o\n", Degree(G1m);
    printf "quotient m-sign G1 equations=%o\n", Equations(Curve(G1m));
    Jm := MinimalModel(Jacobian(G1m));
    printf "quotient m-sign Jacobian=%o\n", Jm;
    if DoRank then
        printf "quotient m-sign Jacobian rank bounds=%o\n", RankBounds(Jm);
    end if;
    printf "quotient m-sign Jacobian torsion=%o\n", Invariants(TorsionSubgroup(Jm));
    for pp in [2,3,5,7,11,13,17,19,23,29,31,37,41] do
        try
            ok, pt := IsLocallySolvable(G1m, pp);
            printf "quotient m-sign local p=%o: %o", pp, ok;
            if ok then printf " pt=%o", pt; end if;
            printf "\n";
        catch e
            printf "quotient m-sign local p=%o failed: %o\n", pp, e`Object;
        end try;
    end for;
    try
        ptsG1m := Points(Curve(G1m) : Bound := PointBound);
        printf "quotient m-sign G1 points Bound %o count=%o\n", PointBound, #ptsG1m;
        if #ptsG1m le 50 then printf "quotient m-sign G1 points=%o\n", ptsG1m; end if;
    catch e
        printf "quotient m-sign G1 Points failed: %o\n", e`Object;
    end try;
catch e
    printf "quotient m-sign GenusOneModel/Jacobian failed: %o\n", e`Object;
end try;

if DoAuto then
    try
        G := AutomorphismGroup(Cp);
        printf "AutomorphismGroup order=%o\n", #G;
        printf "AutomorphismGroup=%o\n", G;
        try
            Qc, pi := CurveQuotient(G);
            printf "full quotient genus=%o curve=%o\n", Genus(Qc), Qc;
            printf "quotient map=%o\n", pi;
        catch e
            printf "CurveQuotient(full group) failed: %o\n", e`Object;
        end try;
        if DoSubgroups then
            try
                subs := Subgroups(G);
                printf "subgroup count=%o\n", #subs;
                for i in [1..#subs] do
                    H := subs[i]`subgroup;
                    try
                        QH, piH := CurveQuotient(H);
                        printf "subgroup %o order=%o quotient genus=%o curve=%o\n", i, #H, Genus(QH), QH;
                    catch e
                        printf "subgroup %o order=%o quotient failed: %o\n", i, #H, e`Object;
                    end try;
                end for;
            catch e
                printf "Subgroups/quotients failed: %o\n", e`Object;
            end try;
        end if;
    catch e
        printf "AutomorphismGroup failed: %o\n", e`Object;
    end try;
end if;

quit;
