// claude_ov_b4_sporadic.m -- Lane 4 (route B4): Galois-stable Richelot from the
// SPORADIC order-11/22 seed bank (BLP2009 Table 1 as corrected in this repo,
// plus the literature order-22 curves), using Magma's RichelotIsogenousSurfaces
// which enumerates ALL Galois-stable (2,2)-kernels, not just rational ones.
//
// For every codomain: exact TorsionSubgroup, 2-rank, and the real-Weil-subfield
// disc signature (RM screen).  Any codomain with torsion containing [2,22] --
// or [44] or [2,2,11], which are absent from the extended DB -- is a hit.
//
// Parallelised with Fork (one child per seed) per code/claude_fork_template.m.
//
// Run: magma -b OUTDIR:=/path claude_ov_b4_sporadic.m

SetColumns(0);
if not assigned MemGB then MemGB := 6; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
if not assigned OUTDIR then OUTDIR := "/tmp/b4sporadic"; end if;
SetMemoryLimit(MemGB*10^9);
Q := Rationals(); P<x> := PolynomialRing(Q);

// ---- seed bank -----------------------------------------------------------
blp := [[-101/48,-61/48,1/4,-5/12], [473/147,-4013/343,6/7,207/49], [8/49,-134/49,3/7,47/49],
 [1159/81,-277/243,40/9,13/27], [-1/13,-191/2197,8/13,15/169], [-28/169,103/2197,3/13,-4/169],
 [594/1805,13348/34295,8/19,-64/361], [208/867,1338/4913,5/17,-39/289], [415/1089,-2207/1089,8/33,119/121],
 [4989/2500,-13599/12500,27/50,-81/250], [-3,59,4,-7], [-163/1215,-367/3645,2/3,13/243],
 [-13/18,71/6,5/3,-13/3], [-2287/27,-1171/9,10/3,-323/3], [121/147,-141/343,2/7,15/49],
 [-1494/847,19480/9317,2/11,-256/121], [125/121,-223/1331,6/11,29/121], [187/361,-649/6859,6/19,23/361]];
blpnames := ["C1","C2","C3","C4corr","C5=X023","C6","C7","C8","C9","C10","Ct1","Ct2","Ct3","Ct4","Ct5","Ct7","Ct8","Ct9"];

seeds := [];   // <name, f>
for i in [1..#blp] do
    a := blp[i][1]; b := blp[i][2]; c := blp[i][3]; d := blp[i][4];
    f := (x^3-x^2+a*x+b)^2 - 4*c^2*(x^2+d)^2;
    Append(~seeds, <"BLP_" cat blpnames[i], f>);
end for;
Append(~seeds, <"DB_1192a", 4*(x^3-2*x^2-x+1) + (x^3+x)^2>);          // generic [22]
Append(~seeds, <"Nicholls_22", (x^2+1)*(4*x^4-12*x^3+13*x^2-4*x+1)>); // ord(D_inf)=11, [2,4]
Append(~seeds, <"Leprevost_f22", (2*x^2-2*x+1)*(2*x^4-2*x^3+x^2-4*x+4)>);
Append(~seeds, <"Petrunin_q22", (16*x^4-16*x^3+x+1)*(x+1)>);
Append(~seeds, <"Flynn62_n1", 16*(x^6 + 2*x^5 + 7/2*x^4 + 2*x^3 + 17/16*x^2 + 3/8*x + 1/16)>);
// positive control: a 2-rank-1 curve with a verified non-rational Galois-stable kernel
Append(~seeds, <"CONTROL_rank1", (x^2+x+1)*(x^4-4*x^2+2)>);

function IntModel(f)
    cs := Coefficients(f);
    d := LCM([Denominator(c) : c in cs]);
    g := f*d;
    n := GCD([Numerator(c) : c in Coefficients(g)]);
    g := g/n;
    c := d/n;                  // f = (1/c)*g
    return (Numerator(1/c)*Denominator(1/c))*g;
end function;

function DiscSig(f, np)
    S := {};
    p := 11; n := 0;
    D := Integers()!Discriminant(HyperellipticCurve(f));
    while n lt np and p lt 400 do
        if D mod p ne 0 then
            try
                chi := FrobeniusPolynomial(BaseChange(Jacobian(HyperellipticCurve(f)), GF(p)));
                if Degree(chi) eq 4 and IsIrreducible(chi) then
                    c3 := Coefficient(chi,3); c2 := Coefficient(chi,2);
                    dd := c3^2 - 4*(c2 - 2*p);
                    if dd ne 0 then S := S join {Squarefree(Numerator(dd))}; n +:= 1; end if;
                end if;
            catch e ;
            end try;
        end if;
        p := NextPrime(p);
    end while;
    return S;
end function;

System("mkdir -p " cat OUTDIR);
NCH := #seeds;
for c in [1..NCH] do
  pid := Fork();
  if pid eq 0 then
    F := Open(OUTDIR cat "/seed" cat IntegerToString(c) cat ".txt", "w");
    nm := seeds[c][1];
    try
        f := IntModel(seeds[c][2]);
        C := HyperellipticCurve(f);
        J := Jacobian(C);
        r0 := #Invariants(TwoTorsionSubgroup(J));
        T0 := Invariants(TorsionSubgroup(J));
        sg := DiscSig(f, 20);
        Puts(F, "SEED " cat nm cat " type=" cat Sprint(Sort([Degree(t[1]) : t in Factorization(f)]))
                cat " 2rank=" cat IntegerToString(r0) cat " torsion=" cat Sprint(T0)
                cat " discsig=" cat Sprint(sg) cat " f=" cat Sprint(f));
        RS := RichelotIsogenousSurfaces(J);
        nb := 0;
        for S in RS do
            if Type(S) ne JacHyp then
                Puts(F, "NB " cat nm cat " NONJAC " cat Sprint(Type(S)));
                continue;
            end if;
            nb +:= 1;
            CS := Curve(S);
            fs := HyperellipticPolynomials(CS);
            r1 := #Invariants(TwoTorsionSubgroup(S));
            T1 := Invariants(TorsionSubgroup(S));
            sg1 := DiscSig(fs, 20);
            Puts(F, "NB " cat nm cat " idx=" cat IntegerToString(nb)
                    cat " type=" cat Sprint(Sort([Degree(t[1]) : t in Factorization(fs)]))
                    cat " 2rank=" cat IntegerToString(r1) cat " torsion=" cat Sprint(T1)
                    cat " discsig=" cat Sprint(sg1)
                    cat (r1 gt r0 select "  *** 2-RANK RAISED ***" else "")
                    cat " f=" cat Sprint(fs));
        end for;
        Puts(F, "DONE " cat nm cat " neighbours=" cat IntegerToString(nb));
    catch e
        Puts(F, "ERROR " cat nm cat " " cat Sprint(e`Object));
    end try;
    Flush(F); delete F;
    quit;
  end if;
end for;
WaitForAllChildren();
for c in [1..NCH] do
    fn := OUTDIR cat "/seed" cat IntegerToString(c) cat ".txt";
    ok, txt := ReadTest(fn);
    if ok then
        for l in Split(txt, "\n") do if #l gt 0 then printf "%o\n", l; end if; end for;
    else
        printf "MISSING %o\n", seeds[c][1];
    end if;
end for;
printf "SPORADIC_DONE seeds=%o\n", NCH;
quit;
