// lane_richelot.m — capped Richelot-isogeny-class closure over the largest
// known split-torsion curves (LMFDB alpha minima + repo anchors), funneling
// every new node.  Targets: 2-group enhancements not in the DB disc box
// ([2,2,2,16]? [8,16]? [2,96]? [2,2,2,10] split? [2,2,2,12] split? ...).
// Usage: magma -b lane_richelot.m > ../logs/richelot.log
SetColumns(0);
SetMemoryLimit(4*10^9);
load "split_lab.m";  // run from product/code/

R<x> := PolynomialRing(Rationals());

// <name, f-coeffs (ascending), h-coeffs (ascending)>
RSEEDS := [*
  <"3.24_8100.e4",   [4,6,-4,-6,5,-3,1],        [0,1,1]>,
  <"2.48_18900.e2",  [1,-13,19,111,19,-13,1],   [0,1,1]>,
  <"4.16_50400.dj1", [252,0,-24,0,-4],          [0,1,0,1]>,
  <"48_1764.a1",     [3,0,9,0,-4],              [0,1,0,1]>,
  <"60_13500.r4",    [1,-7,25,-24,25,-7,1],     [0,1,1]>,
  <"2.2.16_69300.bu1",[9,-33,3,61,3,-33,-11],   [0,1,1]>,
  <"5.10_85500.cj1", [-135,135,11,-33,15,-3,9], [0,1,1]>,
  <"35_7436.a2",     [160,80,166,182,116,30,10],[0,1,1]>,
  <"36_4860.f1",     [9,-9,-1,-5,5,-1,1],       [0,1,1]>,
  <"2.6.6_132300.eo2",[18,-54,9,71,-40,-5,25],  [0,1,1]>,
  <"2.2.2.8_4410.d2",[-45,42,30,-23,-9,3,1],    [0,1,1]>,
  <"40_3168.i1",     [66,0,12,0,2],             [0,1,0,1]>,
  <"2.30_13068.f2",  [12,-36,18,23,-13,-5,1],   [0,1,1]>,
  <"2.24_1440.c2",   [15,0,-6,0,-1],            [0,1,0,1]>,
  <"6.6_196.a1",     [1,3,6,7,6,3,1],           [0,1,1]>,
  <"30_4356.j1",     [2,6,10,11,10,6,2],        [1,0,0,1]>,
  <"2.2.12_1350.d1", [1,-9,19,19,-55,-45],      [0,1,1]>,
  <"3.12_2700.c1",   [-8,12,6,-11,-3,3,1],      [1,1,1]>,
  <"24_300.a1",      [375,0,97,0,8],            [0,1,0,1]>,
  <"28_4732.c1",     [4,10,-4,-1,5,-3,1],       [0,1,1]>,
  <"27_7776.e2",     [0,4,20,17,9,0,1],         [1,1]>,
  <"21_324.a1",      [0,0,1,2,2,1],             [1,1,0,1]>,
  <"16_2700.d1",     [-125,0,97,0,-26,0,2],     [0,1,0,1]>,
  <"2.16_4608.h5",   [2,7,16,18,16,7,2],        [1,1,1,1]>,
  <"5.5_121.a1",     [-2,4,2,5,2,1],            [1,1,1,1]>,
  <"4.12_8100.j1",   [9,-27,53,-62,53,-27,9],   [0,1,1]>,
  <"4.8_225.a2",     [0,-15,-5,-5],             [0,1,0,1]>,
  <"2.4.8_573300.ir1",[7524,-1530,-10562,797,3630,140],[0,1,1]>,
  <"2.20_30492.q2",  [0,-18,42,-20,15,-3,1],    [0,1,1]>,
  <"2.2.10_30492.q1",[-20,-12,66,-2,-45,3,9],   [0,1,1]>,
  <"2.18_64980.bo1", [36,66,12,-5,-3,-3,1],     [0,1,1]>,
  <"2.14_86436.z3",  [10,-6,-24,-5,-9,-3,9],    [0,1,1]>,
  <"25_42336.dw1",   [1,8,13,9,14,4,2],         [0,0,1,1]>,
  <"19_169.a1",      [0,1,1],                   [1,0,1,1]>,
  <"2.10_256.a1",    [0,0,-1,-1,-1,-1],         [1,1,1,1]>,
  // repo anchors
  <"6.12_HLP",       Coefficients(183*(x^2+1)*(32*x^2+61*x+32)*(32*x^2-61*x+32)), [Integers()|]>,
  <"4.16_palin",     [0,2025,11484,9846,11484,2025], [Integers()|]>,
  <"63_lane8",       [-146398496,0,79136353,0,-197570,0,897], [Integers()|]>,
  <"7.7_lane8",      [869675859,0,3232987,0,3025,0,1], [Integers()|]>,
  <"45_lane8",       [168300000000,0,49996210000,0,29240200,0,13981], [Integers()|]>,
  <"70_lane8",       [5184,-23328,37620,-23220,697,3168], [Integers()|]>,
  <"2.24_lane8",     [1,0,46,0,409,0,840],      [Integers()|]>
*];

NODECAP := 24;

for s in RSEEDS do
    f := R!s[2]; h := R!s[3];
    C0 := HyperellipticCurve(f, h);
    seen := { Sprintf("%o", G2Invariants(C0)) };
    frontier := [C0]; all := [C0];
    t0 := Cputime();
    while #frontier gt 0 and #all lt NODECAP do
        newf := [];
        for c in frontier do
            L := [];
            try L := RichelotIsogenousSurfaces(c); catch e L := []; end try;
            for X in L do
                if Type(X) ne CrvHyp then continue; end if;
                key := "";
                try key := Sprintf("%o", G2Invariants(X)); catch e key := "bad"; end try;
                if key eq "bad" or key in seen then continue; end if;
                Include(~seen, key);
                Append(~newf, X); Append(~all, X);
                if #all ge NODECAP then break; end if;
            end for;
            if #all ge NODECAP then break; end if;
        end for;
        frontier := newf;
    end while;
    printf "RICHELOT %o class-size %o (%o s)\n", s[1], #all, Cputime()-t0;
    for i in [2..#all] do
        st := Funnel(all[i], Sprintf("rich|%o|%o", s[1], i));
    end for;
    printf "RICHDONE %o\n", s[1];
end for;
printf "SEARCH_DONE richelot\n";
quit;
