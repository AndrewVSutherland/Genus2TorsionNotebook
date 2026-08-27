// ===========================================================================
// Lane 3 : THE HUMBERT TEST.
//
// Question: is the CF-order-11 locus  Sigma'' = {ord(D_inf) = 11}  on the
// monic chart  y^2 = x(x-1)(x^2+p1x+p2)(x^2+p3x+p4)  CONTAINED in the
// Humbert surface H_5 (RM by sqrt5)?
//
// Sigma'' is a SURFACE in A^4 and H_5 pulls back to a HYPERSURFACE (3-fold),
// so containment is a codimension question that mod-p point counts settle:
//   * if Sigma'' subset H5~ then EVERY point of Sigma''(F_p) satisfies H5 = 0;
//   * if not, then Sigma'' cap H5~ is a curve, so only O(p) of the ~c*p^2
//     points of Sigma''(F_p) satisfy H5 = 0 -- a big, unmistakable signal.
//
// Input: a dump of Sigma''(F_p) produced by
//        code/claude_ov_b2p_scan.c   "surface P dumpfile"   (lines "a b c d").
//
// usage: magma -b P:=61 DUMP:=... NCH:=8 code/claude_ov_l3_humbtest.m
// ===========================================================================
SetColumns(0);
SetMemoryLimit(4*10^9);
if not assigned P then P := 61; elif Type(P) eq MonStgElt then P := StringToInteger(P); end if;
if not assigned NCH then NCH := 8; elif Type(NCH) eq MonStgElt then NCH := StringToInteger(NCH); end if;
if not assigned DUMP then DUMP := Sprintf("/tmp/claude-1000/-home-claude-torsion-jac/1c22f408-33bb-486e-870f-c5fca84bb522/scratchpad/sigma_p%o.txt", P); end if;
if not assigned OUT then OUT := Sprintf("results/claude_ov_l3_humbtest_p%o.log", P); end if;

Z := Integers();  Q := Rationals();
R<J2,J4,J6,J10> := PolynomialRing(Q, 4);
H5Q := R ! eval Read("results/claude_ov_l3_humbert5_eqn.txt");
printf "H5 loaded: #terms=%o\n", #Terms(H5Q);

Fp := GF(P);
Rp<K2,K4,K6,K10> := PolynomialRing(Fp, 4);
H5p := Rp ! H5Q;
Pp<x> := PolynomialRing(Fp);

function H5at(a,b,c,d)
    f := x*(x-1)*(x^2+a*x+b)*(x^2+c*x+d);
    if Degree(f) ne 6 or not IsSquarefree(f) then return -1, 0; end if;
    C := HyperellipticCurve(f);
    I := IgusaClebschInvariants(C);
    v := Evaluate(H5p, [I[1],I[2],I[3],I[4]]);
    return (v eq 0) select 1 else 0, I;
end function;

// ------------------------------------------------- controls over Q (once)
PQ<X> := PolynomialRing(Q);
function H5Q_of(f)
    C := HyperellipticCurve(f);
    I := IgusaClebschInvariants(C);
    return Evaluate(H5Q, [Q!I[1],Q!I[2],Q!I[3],Q!I[4]]) eq 0;
end function;

printf "\n-- control: the 18 BLP2009 order-11 curves (ambient BLP stratum) --\n";
blp := [ [-101/48,-61/48,1/4,-5/12], [473/147,-4013/343,6/7,207/49],
 [8/49,-134/49,3/7,47/49], [1159/81,-277/243,40/9,13/27],
 [-1/13,-191/2197,8/13,15/169], [-28/169,103/2197,3/13,-4/169],
 [594/1805,13348/34295,8/19,-64/361], [208/867,1338/4913,5/17,-39/289],
 [415/1089,-2207/1089,8/33,119/121], [4989/2500,-13599/12500,27/50,-81/250],
 [Q|-3,59,4,-7], [-163/1215,-367/3645,2/3,13/243], [-13/18,71/6,5/3,-13/3],
 [-2287/27,-1171/9,10/3,-323/3], [121/147,-141/343,2/7,15/49],
 [-1494/847,19480/9317,2/11,-256/121], [125/121,-223/1331,6/11,29/121],
 [187/361,-649/6859,6/19,23/361] ];
bn := ["C1","C2","C3","C4corr","C5=X023","C6","C7","C8","C9","C10",
       "Ct1","Ct2","Ct3","Ct4","Ct5","Ct7","Ct8","Ct9"];
nblp := 0;
for i in [1..#blp] do
    r := blp[i];
    f := (X^3-X^2+r[1]*X+r[2])^2 - 4*r[3]^2*(X^2+r[4])^2;
    z := H5Q_of(f);
    if z then nblp +:= 1; end if;
    printf "  BLP %o : onH5=%o\n", bn[i], z;
end for;
printf "CONTROL_BLP onH5=%o of %o\n", nblp, #blp;

printf "\n-- control: the 17 exact-[2,22] chart tuples of the H=60 sweep --\n";
hits := [ [4,27,-3,8], [-3,8,4,27], [-6,32,1,6], [-1/12,-2/9,-7/36,-5/324],
 [-23/12,25/36,-65/36,64/81], [4/15,-5/9,-61/15,128/45], [-15/4,27/4,-19/4,6],
 [11/4,9/4,7/4,4], [-19/4,6,-15/4,27/4], [Q|1,6,-6,32], [-34/15,32/45,31/15,-2/9],
 [-19/6,8/3,-25/6,9/2], [31/15,-2/9,-34/15,32/45], [7/4,4,11/4,9/4],
 [7/6,1/2,13/6,4/3], [13/6,4/3,7/6,1/2], [-25/6,9/2,-19/6,8/3] ];
nh := 0;
for r in hits do
    f := X*(X-1)*(X^2+r[1]*X+r[2])*(X^2+r[3]*X+r[4]);
    if H5Q_of(f) then nh +:= 1; end if;
end for;
printf "CONTROL_HITS onH5=%o of %o\n", nh, #hits;

// ---------------------------------------------------------------- controls
printf "\n-- control: RANDOM points of the chart (ambient stratum) --\n";
nz := 0; nt := 0;
for i in [1..4000] do
    a := Random(Fp); b := Random(Fp); c := Random(Fp); d := Random(Fp);
    r := H5at(a,b,c,d);
    if r lt 0 then continue; end if;
    nt +:= 1; nz +:= r;
end for;
printf "CONTROL_AMBIENT p=%o tested=%o onH5=%o  (expect ~ %o = tested/p)\n",
    P, nt, nz, RealField(4)!(nt/P);

// ------------------------------------------------------------- the test
lines := [l : l in Split(Read(DUMP), "\n") | #l gt 2];
printf "SIGMA_POINTS %o   p^2=%o   ratio=%o\n", #lines, P^2, RealField(5)!(#lines/P^2);

dir := "/tmp/claude_ov_l3_humb";
System("mkdir -p " cat dir cat " && rm -f " cat dir cat "/part*.txt");
for ch in [0..NCH-1] do
    pid := Fork();
    if pid eq 0 then
        F := Open(dir cat "/part" cat IntegerToString(ch) cat ".txt", "w");
        non := 0; on := 0; bad := 0; sample := [];
        for idx in [ch+1..#lines by NCH] do
            s := Split(lines[idx], " ");
            a := Fp!StringToInteger(s[1]); b := Fp!StringToInteger(s[2]);
            c := Fp!StringToInteger(s[3]); d := Fp!StringToInteger(s[4]);
            r := H5at(a,b,c,d);
            if r lt 0 then bad +:= 1;
            elif r eq 1 then on +:= 1;
            else
                non +:= 1;
                if #sample lt 8 then Append(~sample, [a,b,c,d]); end if;
            end if;
        end for;
        Puts(F, Sprintf("PART %o on=%o off=%o bad=%o", ch, on, non, bad));
        for v in sample do
            Puts(F, Sprintf("OFFH5 %o %o %o %o", v[1], v[2], v[3], v[4]));
        end for;
        Flush(F); delete F;
        quit;
    end if;
end for;
WaitForAllChildren();

G := Open(OUT, "w");
on := 0; off := 0; bad := 0;
Puts(G, Sprintf("HUMBTEST p=%o sigma_points=%o", P, #lines));
Puts(G, Sprintf("CONTROL_AMBIENT tested=%o onH5=%o", nt, nz));
for ch in [0..NCH-1] do
    for l in Split(Read(dir cat "/part" cat IntegerToString(ch) cat ".txt"), "\n") do
        if #l eq 0 then continue; end if;
        Puts(G, l);
        s := Split(l, " ");
        if s[1] eq "PART" then
            on  +:= StringToInteger(Split(s[3],"=")[2]);
            off +:= StringToInteger(Split(s[4],"=")[2]);
            bad +:= StringToInteger(Split(s[5],"=")[2]);
        end if;
    end for;
end for;
res := Sprintf("HUMBTEST_RESULT p=%o total=%o onH5=%o offH5=%o bad=%o  frac_on=%o",
    P, on+off+bad, on, off, bad, RealField(6)!(on/Maximum(1,on+off)));
Puts(G, res);
Flush(G); delete G;
printf "%o\n", res;
quit;
