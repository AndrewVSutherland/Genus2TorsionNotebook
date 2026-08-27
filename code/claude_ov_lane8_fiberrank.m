// Lane 8 (2026-07-25, resumed session): ranks of the HLP gluing-condition fibers.
//
// Companion to code/claude_ov_lane8_locusdim.m.  The HLP split locus for the pair
// (M,N) is carried by the surface   S_{M,N} : w^2 = Delta_M(s)*Delta_N(u).
// Fixing s gives the fiber curve  w^2 = core(Delta_M(s)) * core(Delta_N(u)),
// whose genus is (deg core(Delta_N) - 1 or -2)/2.  This script computes the RANK
// of the genus-1 fibers (so: are there infinitely many gluing partners u for a
// given s?), on an integral model.
SetColumns(0);
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
QQ := Rationals();
PU<u> := PolynomialRing(QQ);

function BCgen(N, t)
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

function UnivCubicIn(N, R, t)
    P := PolynomialRing(R); x := P.1;
    b, c := BCgen(N, t);
    return x^3 - b*x^2 + ((1-c)*x - b)^2/4;
end function;

function SqCore(q)
    if q eq 0 then return 0; end if;
    n := Numerator(q)*Denominator(q);
    s := Squarefree(AbsoluteValue(n));
    return n gt 0 select s else -s;
end function;

// squarefree core of Delta_N(u) as an integral polynomial times a rational unit
function CoreDelta(N)
    FT<tt> := RationalFunctionField(QQ);
    D := Discriminant(UnivCubicIn(N, FT, tt));
    num := Numerator(D); den := Denominator(D);
    sf  := SquarefreeFactorization(num*den);
    core := &*[ f[1]^(f[2] mod 2) : f in sf ];
    unit := (num*den) div &*[ f[1]^f[2] : f in sf ];
    return core, unit;
end function;

function DeltaAt(N, t)
    return Discriminant(UnivCubicIn(N, QQ, t));
end function;

targets := [ <"Z/5xZ/10", 10, 10>, <"Z/35", 7, 5>, <"Z/45", 9, 5>, <"Z/3xZ/12", 12, 6>,
             <"Z/2xZ/24", 26, 28>, <"Z/7xZ/7", 7, 7>, <"Z/63", 7, 9> ];

// s-values to test per (M,N): the HLP anchors plus a systematic sweep
svals := AssociativeArray();
svals[10] := [QQ| -1, 3, 4, -7/2, -3/4, -1/3, 2/3, 4/5, -2, 5, -5, 1/2, -1/2, 6, -6, 7 ];
svals[7]  := [QQ| 7, -1, -16/3, -1/6, 6/7, 2, 3, -2, -3, 1/2, -1/2, 4, -4 ];
svals[9]  := [QQ| -5, 4, 2, -2, 3, -3, 1/2, -1/2, 1/3, 5 ];
svals[12] := [QQ| 1/3, -4, 2, -2, 3, 1/2, -1/2, 4 ];
svals[26] := [QQ| 2, 4, -1, 1/2, 3 ];

for T in targets do
    tag := T[1]; M := T[2]; N := T[3];
    coreN, unitN := CoreDelta(N);
    cp := Evaluate(coreN, PU.1);
    dg := Degree(cp);
    printf "\n######## %o : E_%o glued to E_%o ; core(Delta_%o) degree %o\n", tag, M, N, N, dg;
    if dg eq 0 then
        printf "   Delta_%o(u) is ALWAYS a square: the gluing condition is VACUOUS in u.\n", N;
        printf "   => the split locus is the full 2-dimensional (s,u) plane.\n";
        continue;
    end if;
    g := (dg mod 2 eq 0) select (dg-2) div 2 else (dg-1) div 2;
    printf "   generic fiber genus = %o\n", g;
    if g ge 2 then
        printf "   => genus >= 2: by Faltings only FINITELY many gluing partners u for each s.\n";
        continue;
    end if;
    nrk := [];
    for s in svals[M] do
        ds := 0; ok := true;
        try ds := DeltaAt(M, s); catch e ok := false; end try;
        if not ok or ds eq 0 then continue; end if;
        c0 := SqCore(ds);
        c1 := SqCore(QQ!c0 * QQ!unitN);
        rhs := c1 * cp;
        // clear denominators: (d*w)^2 = d^2 * rhs
        d := LCM([Denominator(cc) : cc in Coefficients(rhs)]);
        rhs := d^2 * rhs;
        if g eq 0 then
            // conic w^2 = quadratic(u): rational points dense iff one exists
            Cf := HyperellipticCurve(rhs);
            pts := Points(Cf : Bound := 500);
            printf "   s=%-8o core=%-8o genus 0 : %o small points => %o\n",
                   s, c0, #pts, (#pts gt 0) select "INFINITELY MANY u" else "none found (bound 500)";
            continue;
        end if;
        Cf := HyperellipticCurve(rhs);
        pts := Points(Cf : Bound := 3000);
        if #pts eq 0 then
            printf "   s=%-8o core=%-8o genus 1 : no rational point with bound 3000 (may be a nontrivial torsor)\n", s, c0;
            continue;
        end if;
        Ef := MinimalModel(EllipticCurve(Cf, pts[1]));
        lo, hi := RankBounds(Ef);
        Append(~nrk, lo);
        printf "   s=%-8o core=%-8o genus 1 : E cond %o  rank in [%o,%o]  torsion %o   => %o\n",
               s, c0, Conductor(Ef), lo, hi, Invariants(TorsionSubgroup(Ef)),
               (lo ge 1) select "INFINITELY MANY u" else (hi eq 0 select "FINITELY many u" else "undecided");
    end for;
    if #nrk gt 0 then
        printf "   summary: %o fibers with a rational point, %o of them of rank >= 1\n",
               #nrk, #[r : r in nrk | r ge 1];
    end if;
end for;

printf "SEARCH_DONE fiberrank\n";
quit;
