//////////////////////////////////////////////////////////////////////
//  Genus check for the plane component in the Elkies N=32 reconstruction.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned output_file then
    output_file := "data/elkies32_component_genus.txt";
end if;

Q := Rationals();
P2<Z,R,W> := ProjectiveSpace(Q, 2);

F := 3*Z^4*R^4 + 9*Z^4*R^3*W - 16*Z^3*R^4*W
   + 10*Z^4*R^2*W^2 - 56*Z^3*R^3*W^2 + 32*Z^2*R^4*W^2
   + 5*Z^4*R*W^3 - 72*Z^3*R^2*W^3 + 144*Z^2*R^3*W^3
   - 64*Z*R^4*W^3 + Z^4*W^4 - 40*Z^3*R*W^4
   + 208*Z^2*R^2*W^4 - 224*Z*R^3*W^4 + 80*R^4*W^4
   - 8*Z^3*W^5 + 120*Z^2*R*W^5 - 288*Z*R^2*W^5
   + 160*R^3*W^5 + 24*Z^2*W^6 - 160*Z*R*W^6
   + 160*R^2*W^6 - 32*Z*W^7 + 80*R*W^7 + 16*W^8;

C := Curve(P2, F);
P0 := C![11/6, 1/15, 1];

out := Open(output_file, "w");
fprintf out, "# Elkies N=32 reconstructed component: genus check\n\n";
fprintf out, "Projective plane closure degree: %o\n", Degree(C);
fprintf out, "Irreducible: %o\n", IsIrreducible(C);
fprintf out, "Genus: %o\n", Genus(C);
fprintf out, "Printed point: %o\n", P0;
fprintf out, "Printed point multiplicity: %o\n", Multiplicity(C, P0);

S := SingularPoints(C);
fprintf out, "Rational singular points: %o\n", #S;
for P in S do
    fprintf out, "  %o multiplicity %o\n", P, Multiplicity(C, P);
end for;
fprintf out, "Singular subscheme degree: %o\n", Degree(SingularSubscheme(C));

phi := Parametrization(C, P0);
eqs := DefiningEquations(phi);
fprintf out, "Parametrization from printed point exists; coordinate degrees: %o\n", [Degree(e) : e in eqs];
fprintf out, "The full parametrization is intentionally not printed here; Magma's default one is large.\n";

delete out;
print "Wrote", output_file;
quit;
