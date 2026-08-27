//////////////////////////////////////////////////////////////////////
// Exact different divisor over Q(e), after exact reconstruction.
//////////////////////////////////////////////////////////////////////
print "RELATIVE12_EXACT_DIFFERENT_START";
td := Cputime();
Drel := DifferentDivisor(F12);
Ps,ns := Support(Drel);
print "RELATIVE12_EXACT_DIFFERENT_DEGREE", Degree(Drel),
      "SECONDS",Cputime(td);
for i in [1..#Ps] do
    print "DIFFERENT_PLACE",i,"COEFF",ns[i],"PLACE_DEGREE",Degree(Ps[i]),
          "RAMIFICATION_INDEX",RamificationIndex(Ps[i]),
          "INERTIA_DEGREE",InertiaDegree(Ps[i]),"MINIMUM",Minimum(Ps[i]);
end for;
quit;
