SetColumns(0);
if not assigned factor_index then
    factor_index := 1;
elif Type(factor_index) eq MonStgElt then
    factor_index := StringToInteger(factor_index);
end if;
analysis_mode := "genus";
load "code/contact6_m36_266_r4_plane_normalization.m";
print "AUTOMORPHISM_START", factor_index;
G, autmap := AutomorphismGroup(AF);
print "AUTOMORPHISM_GROUP", G;
print "AUTOMORPHISM_ORDER", Order(G);
print "AUTOMORPHISM_GENERATORS",
      [autmap(G.i) : i in [1..Ngens(G)]];
print "AUTOMORPHISM_DONE";
