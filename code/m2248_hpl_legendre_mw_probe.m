//////////////////////////////////////////////////////////////////////
//  Probe Mordell-Weil generators on the fixed-rho HPL Legendre model.
//////////////////////////////////////////////////////////////////////

Q := Rationals();

rho0 := Q!58466134224 / Q!53109477625;
sigma0 := Q!719363573659505664 / Q!749082246897952705;
tau0 := Q!307598400 / Q!352612321;
dB := Q!72946054224 / Q!53109477625;
Y0 := dB*(tau0^2-1);
V0 := tau0*Y0;

E0 := EllipticCurve([Q!0, -(Q!1+rho0^2), Q!0, rho0^2, Q!0]);
P0 := E0![tau0^2, V0, 1];

print "E0", E0;
print "P0", P0;
print "RankBounds", RankBounds(E0);
print "TorsionSubgroup", TorsionSubgroup(E0);

print "Trying Generators(E0)";
gens := Generators(E0);
print "Generators", gens;
print "number_generators", #gens;

for i in [1..#gens] do
    print "gen", i, gens[i];
end for;

print "Trying MordellWeilGroup(E0)";
G, phi := MordellWeilGroup(E0);
print "MW abstract", G;
print "MW rank", Rank(G);
print "MW invariants", Invariants(G);
for g in Generators(G) do
    print "MW gen", g, phi(g);
end for;
