// probe_ht.m — calibrate Magma height conventions for the full-saturation
// certificate: (1) Height(P) (canonical) vs log H(x(P)); (2) the exact
// semantics of CPSHeightBounds / SilvermanBound (which difference, sign);
// (3) Points(E : Bound := n) semantics (what n bounds).
SetColumns(0);
E := EllipticCurve([0,0,1,-1,0]);   // 37a1
G := E![0,0];
lg := func<x | x eq 0 select -999.0 else Log(x)>;
for n in [1..6] do
    P := n*G;
    xn := Numerator(P[1]); xd := Denominator(P[1]);
    Hx := Max(Abs(xn), Abs(xd));
    printf "n=%o: hhat=%o  logHx=%o  ratio hhat/logHx=%o\n",
        n, Height(P), lg(Hx), Hx eq 1 select 0.0 else Height(P)/lg(Hx);
end for;
ok := true;
try
    lb, ub := CPSHeightBounds(E);
    printf "CPSHeightBounds: lb=%o ub=%o\n", lb, ub;
catch e
    printf "CPSHeightBounds unavailable: %o\n", e`Object;
    ok := false;
end try;
try
    sb := SilvermanBound(E);
    printf "SilvermanBound: %o\n", sb;
catch e
    printf "SilvermanBound unavailable\n";
end try;
// verify the difference direction empirically on multiples of G:
// which of  hhat - logHx  /  hhat - logHx/2  lies within [lb, ub]?
for n in [3..6] do
    P := n*G;
    Hx := Max(Abs(Numerator(P[1])), Abs(Denominator(P[1])));
    printf "n=%o: hhat-logHx=%o  hhat-logHx/2=%o\n",
        n, Height(P) - lg(Hx), Height(P) - lg(Hx)/2;
end for;
// Points bound semantics: count points with given bound
pts := Points(E : Bound := 10);
printf "Points(Bound:=10): %o points; x-coords %o\n", #pts, [ P[1] : P in pts | P[3] ne 0 ];
printf "PROBE_HT_DONE\n";
quit;
