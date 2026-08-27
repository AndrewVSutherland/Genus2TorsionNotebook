//////////////////////////////////////////////////////////////////////
// Exact divisor of M on the degree-12 support field over Q(e).
// The odd valuations over e=0 make the relative signed cover ramified.
// The separate P8 ramification script shows that the P8 map has index 2
// at all points above e=0, making their pullback valuations even.
//////////////////////////////////////////////////////////////////////
print "RELATIVE12_EXACT_M_DIVISOR_START";
tm := Cputime();
Dm := Divisor(melt);
Pm,nm := Support(Dm);
print "RELATIVE12_EXACT_M_DIVISOR_DEGREE",Degree(Dm),
      "SUPPORT_LENGTH",#Pm,"SECONDS",Cputime(tm);
for i in [1..#Pm] do
    print "M_DIVISOR_PLACE",i,"VALUATION",nm[i],
          "PLACE_DEGREE",Degree(Pm[i]),
          "RAMIFICATION_INDEX",RamificationIndex(Pm[i]),
          "INERTIA_DEGREE",InertiaDegree(Pm[i]),
          "MINIMUM",Minimum(Pm[i]);
end for;
ramdeg:=&+[Degree(Pm[i]):i in [1..#Pm]|IsOdd(nm[i])];
grel:=Genus(F12);
gsigned:=2*grel-1+ramdeg div 2;
print "RELATIVE12_M_DIVISOR_ALL_EVEN",&and[IsEven(n):n in nm];
print "RELATIVE12_SIGN_QUADRATIC_RAMIFICATION_DEGREE",ramdeg;
print "RELATIVE12_SUPPORT_GENUS",grel,
      "RELATIVE12_SIGNED_GENUS_BY_RH",gsigned;
assert not &and[IsEven(n):n in nm];
assert ramdeg eq 2 and grel eq 5 and gsigned eq 10;
quit;
