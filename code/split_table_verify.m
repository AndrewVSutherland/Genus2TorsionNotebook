// Verify the displayed non-hyperlinked ("new"/literature) rows of
// split_torsion_table.tex: exact torsion of each equation.
SetColumns(0);
Q := Rationals();
P<x> := PolynomialRing(Q);
rows := [
 <"[45]",  P!0, 13981*x^6+29240200*x^4+49996210000*x^2+168300000000>,
 <"[63]",  P!0, 897*x^6-197570*x^4+79136353*x^2-146398496>,
 <"[70]",  2*x^3-3*x^2-41*x+110, x^3-51*x^2+425*x+179>,
 <"[6,12]", P!0, 132*x^6+396*x^5-6347*x^4-13354*x^3+75207*x^2+81950*x+88825>,
 <"[7,7]", P!0, x^6+3025*x^4+3232987*x^2+869675859>,
 <"[8,8]", P!0, 836*x^6+88596*x^5+88597*x^4+1800118*x^3-4045487*x^2-4535664*x+84285504>,
 <"[2,2,24]", P!0, 581449680*x^6-4134794160*x^5-574778279*x^4+10516435114*x^3-15881777387*x^2+11223443268*x-1729978236>,
 <"[2,2,4,4]", P!0, x*(x+3048806656)*(x+3197902500)*(x+3873070756)*(x+4107528100)>,
 <"[2,2,4,8]", P!0, -144295356865660*x^6+289009358554092*x^5+860505249465645*x^4-1006367755763986*x^3-1568053362370059*x^2+812225828599020*x+949701503960100>
];
for r in rows do
    C := HyperellipticCurve(r[3], r[2]);
    T := TorsionSubgroup(Jacobian(SimplifiedModel(C)));
    printf "ROW %o : computed torsion %o (order %o)\n", r[1], Invariants(T), #T;
end for;
printf "TABLE_VERIFY_DONE\n";
quit;
