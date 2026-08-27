// Lane 8 (overnight 2026-07-25): a general implementation of Howe-Leprevost-Poonen
// (Forum Math. 12 (2000) 315-364) Proposition 4 -- the (2,2)-gluing of two elliptic
// curves along a Galois-module isomorphism of their 2-torsion -- plus HLP's Table 3
// universal elliptic curves E_N^t / E_{2,2N}^t, used to RECOVER explicit genus-2
// models for the split-realized targets the repo has never written down.
//
// Positive controls are run FIRST:
//   (a) the universal curves must have the advertised torsion;
//   (b) the gluing must reproduce HLP's own published Z/63 sextic
//       y^2 = 897x^6-197570x^4+79136353x^2-146398496  from (E_7^{-16/3}, E_9^{4}).
//
// Stage selected by Stage:=1|2|3 (1 = controls, 2 = odd targets Z/35, Z/45,
// 3 = the 2-torsion targets [6,6] control, [2,24], [3,12], [5,10]).
SetColumns(0);
if not assigned Stage then Stage := 1;
elif Type(Stage) eq MonStgElt then Stage := StringToInteger(Stage); end if;
if not assigned MemGB then MemGB := 3;
elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

QQ := Rationals();
P<x> := PolynomialRing(QQ);
Rz<zz> := PolynomialRing(QQ);

// ---------------------------------------------------------------------------
// HLP Table 3: universal elliptic curve y^2 + (1-c)xy - b y = x^3 - b x^2,
// with (0,0) of maximal order.  Completing the square gives the monic cubic model
//     y^2 = f(x) = x^3 - b x^2 + ((1-c)x - b)^2/4 .
// Labels 24/26/28 stand for the X_1(2,2N) families E_{2,4}, E_{2,6}, E_{2,8}.
// ---------------------------------------------------------------------------
function BC(N, t)
    case N:
        when 4:  return t, QQ!0;
        when 5:  return t, t;
        when 6:  return t^2+t, t;
        when 7:  return t^3-t^2, t^2-t;
        when 8:  return 2*t^2-3*t+1, (2*t^2-3*t+1)/t;
        when 9:  return t^5-2*t^4+2*t^3-t^2, t^3-t^2;
        when 10: return (2*t^5-3*t^4+t^3)/(t^2-3*t+1)^2, (-2*t^3+3*t^2-t)/(t^2-3*t+1);
        when 12: return (12*t^6-30*t^5+34*t^4-21*t^3+7*t^2-t)/(t-1)^4,
                        (-6*t^4+9*t^3-5*t^2+t)/(t-1)^3;
        when 24: return t^2-1/16, QQ!0;                                             // E_{2,4}
        when 26: return (-2*t^3+14*t^2-22*t+10)/((t+3)^2*(t-3)^2),
                        (-2*t+10)/((t+3)*(t-3));                                    // E_{2,6}
        when 28: return (16*t^3+16*t^2+6*t+1)/(8*t^2-1)^2,
                        (16*t^3+16*t^2+6*t+1)/(2*t*(4*t+1)*(8*t^2-1));              // E_{2,8}
        else:    return QQ!0, QQ!0;
    end case;
end function;

function UnivCubic(N, t)
    b, c := BC(N, t);
    return x^3 - b*x^2 + ((1-c)*x - b)^2/4;
end function;

// HLP Table 4: discriminant of E_N^t modulo squares (odd N: the cubic field disc).
function DeltaN(N, t)
    case N:
        when 3:  return t*(1-27*t);
        when 4:  return 16*t+1;
        when 5:  return t*(t^2-11*t-1);
        when 6:  return (t+1)*(9*t+1);
        when 7:  return t*(t-1)*(t^3-8*t^2+5*t+1);
        when 8:  return 8*t^2-8*t+1;
        when 9:  return t*(t-1)*(t^2-t+1)*(t^3-6*t^2+3*t+1);
        when 10: return (2*t-1)*(4*t^2-2*t-1);
        when 12: return (2*t^2-2*t+1)*(6*t^2-6*t+1);
        else:    return QQ!0;
    end case;
end function;

function SqCore(q)
    if q eq 0 then return 0; end if;
    n := Numerator(q)*Denominator(q);
    s := Squarefree(AbsoluteValue(n));
    return (n lt 0) select -s else s;
end function;

// ---------------------------------------------------------------------------
// HLP Proposition 4.
// ---------------------------------------------------------------------------
function Prop4(f, g)
    L := SplittingField(f*g);
    PL<X> := PolynomialRing(L);
    ra := [r[1] : r in Roots(PL!f)];
    rb := [r[1] : r in Roots(PL!g)];
    if #ra ne 3 or #rb ne 3 then return [], "roots not separable"; end if;
    Df := QQ!Discriminant(f);  Dg := QQ!Discriminant(g);
    out := [];
    for perm in Permutations({1,2,3}) do
        al := ra;
        be := [rb[perm[1]], rb[perm[2]], rb[perm[3]]];
        a1 := (al[3]-al[2])^2/(be[3]-be[2]) + (al[2]-al[1])^2/(be[2]-be[1])
            + (al[1]-al[3])^2/(be[1]-be[3]);
        b1 := (be[3]-be[2])^2/(al[3]-al[2]) + (be[2]-be[1])^2/(al[2]-al[1])
            + (be[1]-be[3])^2/(al[1]-al[3]);
        a2 := al[1]*(be[3]-be[2]) + al[2]*(be[1]-be[3]) + al[3]*(be[2]-be[1]);
        b2 := be[1]*(al[3]-al[2]) + be[2]*(al[1]-al[3]) + be[3]*(al[2]-al[1]);
        if a1 eq 0 or a2 eq 0 or b1 eq 0 or b2 eq 0 then continue; end if;
        A := Dg*a1/a2;  B := Df*b1/b2;
        // NB the minus sign in HLP's displayed formula for h is GLOBAL: the three
        // big brackets are cyclic images of one another.  (Absorbing it into the
        // first bracket destroys the S3-symmetry and h then never lands in Q[x] --
        // that was the debug signature: only the x^0 and x^6 coefficients stayed
        // rational.)
        hh := -( A*(al[2]-al[1])*(al[1]-al[3])*X^2 + B*(be[2]-be[1])*(be[1]-be[3]))
             * ( A*(al[3]-al[2])*(al[2]-al[1])*X^2 + B*(be[3]-be[2])*(be[2]-be[1]))
             * ( A*(al[1]-al[3])*(al[3]-al[2])*X^2 + B*(be[1]-be[3])*(be[3]-be[2]));
        cs := Coefficients(hh);
        ok := true; csQ := [];
        for cc in cs do
            fl, q := IsCoercible(QQ, cc);
            if not fl then ok := false; break; end if;
            Append(~csQ, q);
        end for;
        if not ok then continue; end if;
        hQ := P!csQ;
        if Degree(hQ) ne 6 or not IsSquarefree(hQ) then continue; end if;
        Append(~out, <perm, hQ>);
    end for;
    return out, "ok";
end function;

// clear denominators (y -> y/d) and drop square factors of the content
function CleanSextic(h)
    d := LCM([Denominator(c) : c in Coefficients(h)]);
    h2 := P!(d^2*h);
    cont := GCD([Integers()!c : c in Coefficients(h2)]);
    sq := 1;
    for pr in Factorization(cont) do sq *:= pr[1]^(2*(pr[2] div 2)); end for;
    return P!(h2/sq);
end function;

procedure Analyse(tag, h)
    print "--------------------------------------------------";
    print tag;
    hc := CleanSextic(h);
    print "MODEL  y^2 =", hc;
    C := HyperellipticCurve(hc);
    print "genus", Genus(C);
    J := Jacobian(C);
    T, mp := TorsionSubgroup(J);
    print "TORSION invariants:", Invariants(T), "order", #T;
    print "2-torsion:", Invariants(TwoTorsionSubgroup(J));
    print "G2Invariants:", G2Invariants(C);
    dsc := Integers()!Discriminant(hc);
    badp := Set(PrimeDivisors(dsc)) join {2};
    nirr := 0; ntest := 0; cores := {};
    for p in PrimesInInterval(3,200) do
        if p in badp then continue; end if;
        Cp := ChangeRing(C, GF(p));
        if not IsNonsingular(Cp) then continue; end if;
        chi := Rz!Reverse(Coefficients(LPolynomial(Cp)));
        ntest +:= 1;
        if IsIrreducible(chi) then nirr +:= 1; end if;
        dd := Integers()!(Coefficient(chi,3)^2 - 4*(Coefficient(chi,2)-2*p));
        Include(~cores, dd eq 0 select 0 else Squarefree(dd));
    end for;
    print "good primes 3..200:", ntest, " irreducible chi:", nirr,
          " RM/split-screen cores:", Sort(Setseq(cores));
end procedure;

procedure RunPair(tag, N1, t1, N2, t2)
    f := UnivCubic(N1, t1);  g := UnivCubic(N2, t2);
    printf "\n############ %o ############\n", tag;
    if Discriminant(f) eq 0 or Discriminant(g) eq 0 then
        print "DEGENERATE parameter (zero discriminant) -- skipped"; return;
    end if;
    printf "E = E_%o^(%o)   f = %o\n", N1, t1, f;
    printf "F = E_%o^(%o)   g = %o\n", N2, t2, g;
    printf "disc f core %o ; disc g core %o\n",
        SqCore(Discriminant(f)), SqCore(Discriminant(g));
    E := EllipticCurve(f); F := EllipticCurve(g);
    printf "E: cond %o torsion %o rank %o ; F: cond %o torsion %o rank %o\n",
        Conductor(E), Invariants(TorsionSubgroup(E)), Rank(E),
        Conductor(F), Invariants(TorsionSubgroup(F)), Rank(F);
    if jInvariant(E) eq jInvariant(F) then print "WARNING: equal j-invariants"; end if;
    res, msg := Prop4(f, g);
    printf "Prop4 returned %o rational sextics (%o)\n", #res, msg;
    for r in res do
        Analyse(Sprintf("%o -- root labelling %o", tag, r[1]), r[2]);
    end for;
end procedure;

// ===========================================================================
if Stage eq 1 then
    print "########## CONTROL A: universal curves have the advertised torsion ##########";
    for spec in [<5,QQ!(1/26)>, <5,QQ!(93/10)>, <5,QQ!3>, <7,QQ!(-1)>, <7,QQ!(-16/3)>,
                 <7,QQ!7>, <7,QQ!(-14/13)>, <9,QQ!4>, <9,QQ!(-5)>, <6,QQ!2>,
                 <10,QQ!2>, <12,QQ!3>, <24,QQ!2>, <26,QQ!2>, <28,QQ!2>] do
        N := spec[1]; t := spec[2];
        f := UnivCubic(N,t);
        if Discriminant(f) eq 0 then printf "N=%o t=%o DEGENERATE\n", N, t; continue; end if;
        E := EllipticCurve(f);
        printf "E_%o^(%o): torsion %o   disc-core %o   Delta_N-core %o\n",
            N, t, Invariants(TorsionSubgroup(E)), SqCore(Discriminant(f)),
            SqCore(DeltaN(N lt 20 select N else 0, t));
    end for;

    print "";
    print "########## CONTROL B: rebuild HLP's published Z/63 curve ##########";
    Cpub := HyperellipticCurve(897*x^6 - 197570*x^4 + 79136353*x^2 - 146398496);
    print "published Z/63 G2Invariants:", G2Invariants(Cpub);
    RunPair("CONTROL Z/63 = E_7^{-16/3} glued to E_9^{4}", 7, QQ!(-16/3), 9, QQ!4);

    print "";
    print "########## CONTROL C: rebuild HLP's published Z/7 x Z/7 curve ##########";
    Cpub2 := HyperellipticCurve(x^6 + 3025*x^4 + 3232987*x^2 + 869675859);
    print "published [7,7] G2Invariants:", G2Invariants(Cpub2);
    RunPair("CONTROL [7,7] = E_7^{7} glued to E_7^{-14/13}", 7, QQ!7, 7, QQ!(-14/13));
end if;

if Stage eq 2 then
    RunPair("TARGET Z/35 = E_7^{-1} glued to E_5^{1/26}", 7, QQ!(-1), 5, QQ!(1/26));
    RunPair("TARGET Z/45 = E_9^{-5} glued to E_5^{93/10}", 9, QQ!(-5), 5, QQ!(93/10));
end if;

if Stage eq 3 then
    RunPair("CONTROL [6,6] = E_{2,6}^{2} glued to E_{2,6}^{4}", 26, QQ!2, 26, QQ!4);
    RunPair("TARGET [2,24] = E_{2,6}^{2} glued to E_{2,8}^{2}", 26, QQ!2, 28, QQ!2);
    RunPair("TARGET [2,24] = E_{2,6}^{4} glued to E_{2,8}^{3}", 26, QQ!4, 28, QQ!3);
    RunPair("TARGET [2,24] = E_{2,6}^{7} glued to E_{2,8}^{1}", 26, QQ!7, 28, QQ!1);

    print "\n########## solving Delta_12(t) = Delta_6(u) mod squares ##########";
    sols312 := [];
    for tn in [-6..6] do for td in [1..4] do
        if GCD(tn,td) ne 1 then continue; end if;
        t := QQ!tn/td;
        d12 := DeltaN(12,t); if d12 eq 0 then continue; end if;
        for un in [-40..40] do for ud in [1..8] do
            if GCD(un,ud) ne 1 then continue; end if;
            u := QQ!un/ud;
            d6 := DeltaN(6,u); if d6 eq 0 then continue; end if;
            if IsSquare(d12*d6) then Append(~sols312, <t,u>); end if;
        end for; end for;
    end for; end for;
    printf "found %o (t,u) with Delta_12(t)Delta_6(u) square; first few: %o\n",
        #sols312, sols312[1..Min(12,#sols312)];
    done := 0;
    for s in sols312 do
        if done ge 3 then break; end if;
        f := UnivCubic(12,s[1]); g := UnivCubic(6,s[2]);
        if Discriminant(f) eq 0 or Discriminant(g) eq 0 then continue; end if;
        if jInvariant(EllipticCurve(f)) eq jInvariant(EllipticCurve(g)) then continue; end if;
        RunPair(Sprintf("TARGET [3,12] = E_12^{%o} glued to E_6^{%o}", s[1], s[2]),
                12, s[1], 6, s[2]);
        done +:= 1;
    end for;

    print "\n########## solving Delta_10(t) = Delta_10(u) mod squares ##########";
    sols510 := [];
    for tn in [-30..30] do for td in [1..6] do
        if GCD(tn,td) ne 1 then continue; end if;
        t := QQ!tn/td; d1 := DeltaN(10,t); if d1 eq 0 then continue; end if;
        for un in [-30..30] do for ud in [1..6] do
            if GCD(un,ud) ne 1 then continue; end if;
            u := QQ!un/ud; if u eq t then continue; end if;
            d2 := DeltaN(10,u); if d2 eq 0 then continue; end if;
            if IsSquare(d1*d2) then Append(~sols510, <t,u>); end if;
        end for; end for;
    end for; end for;
    printf "found %o (t,u) with Delta_10(t)Delta_10(u) square, t<>u; first few: %o\n",
        #sols510, sols510[1..Min(12,#sols510)];
    done := 0;
    for s in sols510 do
        if done ge 3 then break; end if;
        f := UnivCubic(10,s[1]); g := UnivCubic(10,s[2]);
        if Discriminant(f) eq 0 or Discriminant(g) eq 0 then continue; end if;
        if jInvariant(EllipticCurve(f)) eq jInvariant(EllipticCurve(g)) then continue; end if;
        RunPair(Sprintf("TARGET [5,10] = E_10^{%o} glued to E_10^{%o}", s[1], s[2]),
                10, s[1], 10, s[2]);
        done +:= 1;
    end for;
end if;


if Stage eq 4 then
    // Banks of split anchors: sweep the HLP gluing conditions Delta_M(s) = Delta_N(t)
    // mod squares and record every (s,t) for which Prop 4 actually returns a rational
    // sextic (equal discriminant class is necessary but NOT sufficient -- the cubic
    // fields must be isomorphic; Prop4's rationality test is the exact criterion).
    procedure Bank(tag, M, N, srange, trange, cap)
        printf "\n########## BANK %o : Delta_%o(s) = Delta_%o(t) mod squares ##########\n", tag, M, N;
        nsq := 0; nhit := 0;
        for s in srange do
            if nhit ge cap then break; end if;
            ds := DeltaN(M, s); if ds eq 0 then continue; end if;
            fs := UnivCubic(M, s); if Discriminant(fs) eq 0 then continue; end if;
            for t in trange do
                if nhit ge cap then break; end if;
                dt := DeltaN(N, t); if dt eq 0 then continue; end if;
                if not IsSquare(ds*dt) then continue; end if;
                nsq +:= 1;
                ft := UnivCubic(N, t); if Discriminant(ft) eq 0 then continue; end if;
                if jInvariant(EllipticCurve(fs)) eq jInvariant(EllipticCurve(ft)) then continue; end if;
                res := Prop4(fs, ft);
                if #res eq 0 then continue; end if;
                for r in res do
                    hc := CleanSextic(r[2]);
                    C := HyperellipticCurve(hc);
                    J := Jacobian(C);
                    T := TorsionSubgroup(J);
                    inv := Invariants(T);
                    nhit +:= 1;
                    printf "ANCHOR %o  s=%o t=%o perm=%o  torsion=%o  y^2 = %o\n",
                        tag, s, t, r[1], inv, hc;
                end for;
            end for;
        end for;
        printf "BANK %o: %o discriminant-matching (s,t) pairs, %o rational sextics built\n", tag, nsq, nhit;
    end procedure;

    HT := function(H)
        L := [];
        for d in [1..H] do for n in [-H*4..H*4] do
            if GCD(n,d) eq 1 then Append(~L, Rationals()!n/d); end if;
        end for; end for;
        return Sort(Setseq(Seqset(L)));
    end function;

    // Z/35 = E_7^s glued to E_5^t  (HLP: s=-1, t=1/26 is one point of a rank-2 curve)
    Bank("Z35", 7, 5, HT(6), HT(30), 40);
    // Z/45 = E_9^s glued to E_5^t  (HLP: s=-5, t=93/10)
    Bank("Z45", 9, 5, HT(6), HT(30), 40);
    // [7,7] = E_7^s glued to E_7^t  (HLP: s=7, t=-14/13)
    Bank("77", 7, 7, HT(8), HT(16), 40);
    // Z/63 = E_7^s glued to E_9^t   (HLP: s=-16/3, t=4)  -- the P^0 entry
    Bank("Z63", 7, 9, HT(12), HT(12), 40);
end if;

print "SEARCH_DONE stage", Stage;
quit;
