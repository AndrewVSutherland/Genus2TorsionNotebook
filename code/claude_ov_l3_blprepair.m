// ===========================================================================
// Lane 3 : re-run the ambient BLP control with the REPAIRED RM screen.
//
// The claim "15 of the 18 BLP2009 order-11 curves have End = Z" (addendum 2 of
// the [2,22] plan note) was computed with the NAIVE real-subfield-disc census,
// which code/claude_ov_l3_h5diag.m shows is wrong for every curve whose
// geometric RM is not defined over Q (a1 = 0 at the inert primes).  So the
// control has to be recomputed with the repaired screen:
//     a prime is informative only if  a1 <> 0,  n_p <> 0,  n_p not a square.
//
// Also reports H5 (the Humbert-5 equation) for each curve, and, for the curves
// the repaired screen calls RM, the discriminant it finds.
//
// usage: code/claude_magma_slot.sh -b code/claude_ov_l3_blprepair.m
// ===========================================================================
SetColumns(0);
SetMemoryLimit(4*10^9);
if not assigned PB then PB := 600; elif Type(PB) eq MonStgElt then PB := StringToInteger(PB); end if;
if not assigned OUT then OUT := "results/claude_ov_l3_blprepair.log"; end if;

Z := Integers();  Q := Rationals();
X := PolynomialRing(Q).1;
Px<x> := PolynomialRing(Q);
R<J2,J4,J6,J10> := PolynomialRing(Q, 4);
H5 := R ! eval Read("results/claude_ov_l3_humbert5_eqn.txt");

function Integralise(f)
    m := LCM([Denominator(c) : c in Coefficients(f)] cat [Z|1]);
    return Px ! [ Z ! (m^2 * c) : c in Coefficients(m^4 * f) ];   // f -> m^6 f
end function;

function Census(f, B)
    g := Integralise(f);
    dsc := Z ! Discriminant(g);
    dn := {};  dr := {};  np := 0; na1z := 0; nsq := 0;
    for p in PrimesInInterval(3, B) do
        if dsc mod p eq 0 then continue; end if;
        fp := PolynomialRing(GF(p)) ! g;
        if Degree(fp) lt 5 or not IsSquarefree(fp) then continue; end if;
        L := LPolynomial(HyperellipticCurve(fp));
        a1 := Z!(-Coefficient(L,1));  a2 := Z!Coefficient(L,2);
        n := a1^2 - 4*(a2 - 2*p);
        np +:= 1;
        if n eq 0 then continue; end if;
        d := SquarefreeFactorization(n);
        if d eq 1 then nsq +:= 1; continue; end if;
        Include(~dn, d);
        if a1 eq 0 then na1z +:= 1; continue; end if;
        Include(~dr, d);
    end for;
    return dn, dr, np, na1z, nsq;
end function;

function H5Zero(f)
    C := HyperellipticCurve(Integralise(f));
    I := IgusaClebschInvariants(C);
    return Evaluate(H5, [Q!I[1],Q!I[2],Q!I[3],Q!I[4]]) eq 0;
end function;

G := Open(OUT, "w");
procedure Say(G, s)
    printf "%o\n", s;  Puts(G, s);
end procedure;

Say(G, Sprintf("BLP_REPAIR PB=%o", PB));

// ------------------------------------------- the 18 BLP2009 order-11 curves
// f = (X^3 - X^2 + a X + b)^2 - 4 c^2 (X^2 + d)^2
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

nRMnaive := 0; nRMrep := 0; nH5 := 0;
for i in [1..#blp] do
    r := blp[i];
    f := (X^3-X^2+r[1]*X+r[2])^2 - 4*r[3]^2*(X^2+r[4])^2;
    dn, dr, np, na1z, nsq := Census(f, PB);
    h5 := H5Zero(f);
    if #dn eq 1 then nRMnaive +:= 1; end if;
    if #dr eq 1 then nRMrep  +:= 1; end if;
    if h5 then nH5 +:= 1; end if;
    Say(G, Sprintf("BLP %-8o H5zero=%o | naive #=%o %o | REPAIRED #=%o %o | np=%o na1z=%o nsq=%o",
        bn[i], h5, #dn, (#dn le 3) select Sort(Setseq(dn)) else "scatter",
        #dr, (#dr le 3) select Sort(Setseq(dr)) else "scatter", np, na1z, nsq));
end for;
Say(G, Sprintf("BLP_SUMMARY n=%o  naive-RM=%o  REPAIRED-RM=%o  onH5=%o",
    #blp, nRMnaive, nRMrep, nH5));

// ---------------------------------- the 17 exact-[2,22] tuples of the sweep
hits := [ [4,27,-3,8], [-3,8,4,27], [-6,32,1,6], [-1/12,-2/9,-7/36,-5/324],
 [-23/12,25/36,-65/36,64/81], [4/15,-5/9,-61/15,128/45], [-15/4,27/4,-19/4,6],
 [11/4,9/4,7/4,4], [-19/4,6,-15/4,27/4], [Q|1,6,-6,32], [-34/15,32/45,31/15,-2/9],
 [-19/6,8/3,-25/6,9/2], [31/15,-2/9,-34/15,32/45], [7/4,4,11/4,9/4],
 [7/6,1/2,13/6,4/3], [13/6,4/3,7/6,1/2], [-25/6,9/2,-19/6,8/3] ];
nh := 0; nhr := 0;
for r in hits do
    f := X*(X-1)*(X^2+r[1]*X+r[2])*(X^2+r[3]*X+r[4]);
    dn, dr, np, na1z, nsq := Census(f, PB);
    h5 := H5Zero(f);
    if h5 then nh +:= 1; end if;
    if dr eq {5} then nhr +:= 1; end if;
    Say(G, Sprintf("HIT %o H5zero=%o naive=%o REPAIRED=%o np=%o na1z=%o",
        r, h5, Sort(Setseq(dn)), Sort(Setseq(dr)), np, na1z));
end for;
Say(G, Sprintf("HIT_SUMMARY n=%o onH5=%o repaired-core5=%o", #hits, nh, nhr));

Flush(G); delete G;
printf "BLPREPAIR_DONE\n";
quit;
