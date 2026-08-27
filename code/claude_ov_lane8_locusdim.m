// Lane 8 (2026-07-25, resumed session): DIMENSION AND DEGREE of the HLP split
// torsion locus, and the structural dichotomy that controls it.
//
// Question (lane item 1): every recovered HLP model is split.  Put the split
// point in a chart, identify the split divisor through it, and say whether the
// torsion locus has a component transverse to the split locus.
//
// What this script measures:
//  (A) For each universal family E_N^t, the factorization type of the
//      2-division cubic over Q(t).  If it is REDUCIBLE (i.e. N even, so E has a
//      rational 2-torsion point) then "E_N^s[2] iso E_{N'}^u[2]" is EQUIVALENT
//      to "Delta_N(s)*Delta_{N'}(u) is a square" -- ONE condition.  If it is
//      IRREDUCIBLE (N odd) the discriminant match is only necessary; the two
//      cubic FIELDS must agree, which is one more condition.
//  (B) The measured ratio  #(disc-matching pairs) : #(same-cubic-field pairs)
//      : #(rational Prop-4 sextics) in a height box, for each target.
//  (C) For a fixed s, the fiber curve  w^2 = core(Delta_N(s)) * Delta_{N'}(u):
//      genus, and (genus 1) rank.  Positive rank => infinitely many u for that
//      s => the disc-matching locus is a surface with Zariski-dense rational
//      points => the HLP locus is dense in the whole 2-dimensional split locus.
SetColumns(0);
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
if not assigned HB then HB := 8; elif Type(HB) eq MonStgElt then HB := StringToInteger(HB); end if;
SetMemoryLimit(MemGB*10^9);

QQ := Rationals();

// ---- universal families (Tate normal form data), verbatim from prop4b.m ----
function BCgen(N, R, t)
    case N:
        when 5:  return t, t;
        when 6:  return t^2+t, t;
        when 7:  return t^3-t^2, t^2-t;
        when 9:  return t^5-2*t^4+2*t^3-t^2, t^3-t^2;
        when 10: return (2*t^5-3*t^4+t^3)/(t^2-3*t+1)^2, (-2*t^3+3*t^2-t)/(t^2-3*t+1);
        when 12: return (12*t^6-30*t^5+34*t^4-21*t^3+7*t^2-t)/(t-1)^4,
                        (-6*t^4+9*t^3-5*t^2+t)/(t-1)^3;
        when 26: return (-2*t^3+14*t^2-22*t+10)/((t+3)^2*(t-3)^2), (-2*t+10)/((t+3)*(t-3));
        when 28: return (16*t^3+16*t^2+6*t+1)/(8*t^2-1)^2,
                        (16*t^3+16*t^2+6*t+1)/(2*t*(4*t+1)*(8*t^2-1));
    end case;
end function;

function UnivCubicGen(N, R, t)
    P := PolynomialRing(R); x := P.1;
    b, c := BCgen(N, R, t);
    return x^3 - b*x^2 + ((1-c)*x - b)^2/4;
end function;

function SqCore(q)
    if q eq 0 then return 0; end if;
    n := Numerator(q)*Denominator(q);
    s := Squarefree(AbsoluteValue(n));
    return n gt 0 select s else -s;
end function;

// ================= (A) generic factorization of the 2-division cubic ========
printf "########## (A) 2-division cubic of E_N^t over Q(t) ##########\n";
FT<tt> := RationalFunctionField(QQ);
DeltaPoly := AssociativeArray();
for N in [5,6,7,9,10,12,26,28] do
    cub := UnivCubicGen(N, FT, tt);
    fac := Factorization(cub);
    D := Discriminant(cub);
    num := Numerator(D); den := Denominator(D);
    sf  := SquarefreeFactorization(num*den);
    core := &*[ f[1]^(f[2] mod 2) : f in sf ];
    lead := (num*den) div &*[ f[1]^f[2] : f in sf ];   // dropped unit: reconstruct
    DeltaPoly[N] := <core, lead>;
    printf "  N=%2o : cubic factor degrees %o    Delta_N(t) squarefree core = (%o) * (%o)  [deg %o]\n",
           N, [<Degree(f[1]), f[2]> : f in fac], lead, core, Degree(core);
end for;

// ================= (B) disc-match vs field-match vs rational sextic =========
function DeltaAt(N, t)
    R := PolynomialRing(QQ);
    cub := UnivCubicGen(N, QQ, t);
    return Discriminant(cub), cub;
end function;

HT := function(H)
    L := [];
    for d in [1..H] do for n in [-H*4..H*4] do
        if GCD(n,d) eq 1 then Append(~L, QQ!n/d); end if;
    end for; end for;
    return Sort(Setseq(Seqset(L)));
end function;

printf "\n########## (B) disc-match : field-match : rational-sextic, height box %o ##########\n", HB;
pairs := [ <"Z/5xZ/10", 10, 10>, <"Z/7xZ/7", 7, 7>, <"Z/63", 7, 9>, <"Z/45", 9, 5>, <"Z/35", 7, 5>, <"Z/3xZ/12", 12, 6> ];
box := HT(HB);
for pr in pairs do
    tag := pr[1]; M := pr[2]; N := pr[3];
    nsq := 0; nfield := 0; ndiffj := 0;
    // precompute
    dataM := []; dataN := [];
    for s in box do
        ok := true; d := 0; c := 0;
        try d, c := DeltaAt(M, s); catch e ok := false; end try;
        if ok and d ne 0 then Append(~dataM, <s, d, c>); end if;
    end for;
    for u in box do
        ok := true; d := 0; c := 0;
        try d, c := DeltaAt(N, u); catch e ok := false; end try;
        if ok and d ne 0 then Append(~dataN, <u, d, c>); end if;
    end for;
    for a in dataM do
        for b in dataN do
            if not IsSquare(a[2]*b[2]) then continue; end if;
            nsq +:= 1;
            EA := EllipticCurve(a[3]); EB := EllipticCurve(b[3]);
            if jInvariant(EA) eq jInvariant(EB) then continue; end if;
            ndiffj +:= 1;
            // same 2-torsion field?
            fa := a[3]; fb := b[3];
            faca := Factorization(fa); facb := Factorization(fb);
            dega := Sort([Degree(f[1]) : f in faca]); degb := Sort([Degree(f[1]) : f in facb]);
            same := false;
            if dega eq degb then
                if dega eq [1,1,1] then
                    same := true;
                elif dega eq [1,2] then
                    qa := [f[1] : f in faca | Degree(f[1]) eq 2][1];
                    qb := [f[1] : f in facb | Degree(f[1]) eq 2][1];
                    same := IsSquare(Discriminant(qa)*Discriminant(qb));
                else
                    same := IsIsomorphic(NumberField(fa/LeadingCoefficient(fa)),
                                         NumberField(fb/LeadingCoefficient(fb)));
                end if;
            end if;
            if same then nfield +:= 1; end if;
        end for;
    end for;
    printf "  %-9o E_%o x E_%o : #box=%o x %o  disc-matching pairs=%o  (distinct j)=%o  SAME 2-TORSION FIELD=%o   ratio %o\n",
           tag, M, N, #dataM, #dataN, nsq, ndiffj, nfield,
           ndiffj eq 0 select "n/a" else RealField(4)!(nfield/ndiffj);
end for;

// ================= (C) fiber curves: genus and rank =========================
printf "\n########## (C) fiber curve  w^2 = core(Delta_M(s)) * Delta_N(u)  ##########\n";
PU<u> := PolynomialRing(QQ);
anchors := [ <"Z/5xZ/10", 10, 10, [QQ|-1, 3, 4, -7/2]>,
             <"Z/7xZ/7",   7,  7, [QQ|7, -1/6, 6/7]>,
             <"Z/63",      7,  9, [QQ|-16/3]>,
             <"Z/45",      9,  5, [QQ|-5, 1/6, 6/5]>,
             <"Z/35",      7,  5, [QQ|-1, -10, 2]> ];
for an in anchors do
    tag := an[1]; M := an[2]; N := an[3];
    coreN := DeltaPoly[N][1];  leadN := DeltaPoly[N][2];
    // core polynomial in u with its unit; substitute u
    cp := Evaluate(coreN, PU.1);
    printf "  --- %o : E_%o glued to E_%o ; squarefree core of Delta_%o(u) has degree %o\n",
           tag, M, N, N, Degree(cp);
    for s in an[4] do
        ds := 0; ok := true;
        try ds := DeltaAt(M, s); catch e ok := false; end try;
        if not ok or ds eq 0 then continue; end if;
        c0 := SqCore(ds);
        // also fold in the (constant-in-u) unit leadN
        c1 := SqCore(QQ!c0 * QQ!leadN);
        rhs := c1 * cp;
        // strip square factors / content to get a clean model
        rhs := rhs / LeadingCoefficient(rhs) * SqCore(Numerator(LeadingCoefficient(rhs))*Denominator(LeadingCoefficient(rhs)));
        dg := Degree(rhs);
        g := (dg mod 2 eq 0) select (dg-2) div 2 else (dg-1) div 2;
        printf "      s=%-8o core(Delta_%o(s))=%-8o  fiber curve w^2 = %o   deg %o  genus %o\n",
               s, M, c0, rhs, dg, g;
        if g eq 1 then
            try
                Cf := HyperellipticCurve(rhs);
                if #Points(Cf : Bound := 2000) gt 0 then
                    Ef, mp := EllipticCurve(Cf, Points(Cf : Bound := 2000)[1]);
                    Ef := MinimalModel(Ef);
                    rk, prv := RankBounds(Ef);
                    printf "         genus-1 fiber: cond %o  rank in [%o,%o]  torsion %o\n",
                           Conductor(Ef), rk, prv, Invariants(TorsionSubgroup(Ef));
                else
                    printf "         genus-1 fiber: NO small rational point found (bound 2000)\n";
                end if;
            catch e printf "         genus-1 fiber: EllipticCurve failed (%o)\n", e`Object;
            end try;
        elif g eq 0 then
            Cf := Conic(ProjectiveSpace(QQ,2), 0);
            printf "         genus-0 fiber: rational points are dense iff there is one (deg %o)\n", dg;
        else
            printf "         genus %o fiber: FINITELY many u (Faltings)\n", g;
        end if;
    end for;
end for;

printf "SEARCH_DONE locusdim\n";
quit;
