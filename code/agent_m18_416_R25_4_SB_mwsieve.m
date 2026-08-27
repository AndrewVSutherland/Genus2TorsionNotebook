//////////////////////////////////////////////////////////////////////
// Mordell-Weil sieve for the R = -25/4 S_B square-condition cover.
//
// The cover is the double cover of the rank-3 elliptic fiber cut out by
//
//     Y^2 = 2*alpha_B(m)*m^2 + 2*h_B(m)*m^2*g,
//
// with g on the normalized quartic fiber.  A rational point on the
// genus-5 cover gives, for every good prime p, a point of E(F_p) where
// this value is a square.  This script sieves classes in E(Q)/N E(Q).
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();

if not assigned N then N := 8; end if;
if Type(N) eq MonStgElt then N := StringToInteger(N); end if;
if not assigned MaxPrime then MaxPrime := 199; end if;
if Type(MaxPrime) eq MonStgElt then MaxPrime := StringToInteger(MaxPrime); end if;
if not assigned OnlyCuts then OnlyCuts := false; end if;
if Type(OnlyCuts) eq MonStgElt then OnlyCuts := OnlyCuts in {"true", "True", "1", "yes"}; end if;

PrimeList := [11,13,19,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199];
if assigned MaxPrime then
    PrimeList := [p : p in PrimesUpTo(MaxPrime) | p gt 7];
end if;
if assigned Primes then
    PrimeList := [StringToInteger(s) : s in Split(Primes, ",")];
end if;

P<m> := PolynomialRing(Q);
Gamma := 1024*m^4 - 865600*m^2 + 231800625;
K := Q!15225/32;

alphaP := 29/100*m^4 - 16907/64*m^2 + 268888725/4096; // alpha_B*m^2
hP := 2/525*m^4 - 29/8*m^2 + 441525/512;              // h_B*m^2

E := EllipticCurve([Q!0, Q!-1, Q!0, Q!-468441, Q!122191641]);
TgensQ := [E![361,0,1], E![429,0,1]];
FgensQ := [E![545,-5336,1], E![-489,-15300,1], E![Q!8289/25,Q!224112/125,1]];

// Base point on the quartic fiber corresponding to the elliptic origin.
okBase, gBase := IsSquare(Evaluate(Gamma, Q!-25/2));
assert okBase;
mBase := Q!-25/2;

function RedRat(q, F)
    q := Q!q;
    den := Denominator(q);
    if den mod (Z!Characteristic(F)) eq 0 then
        return false, F!0;
    end if;
    return true, F!Numerator(q) / F!den;
end function;

function RedPoint(Pt, EF)
    F := BaseRing(EF);
    if Pt eq Parent(Pt)!0 then
        return true, EF!0;
    end if;
    xq := Pt[1]/Pt[3];
    yq := Pt[2]/Pt[3];
    okx, xr := RedRat(xq, F);
    oky, yr := RedRat(yq, F);
    if not (okx and oky) then
        return false, EF!0;
    end if;
    return true, EF![xr, yr, F!1];
end function;

function AddGeneratedSubgroup(gens, EF)
    S := {EF!0};
    for g in gens do
        if g eq EF!0 then
            continue;
        end if;
        T := {};
        og := Order(g);
        for P0 in S do
            for k in [0..og-1] do
                Include(~T, P0 + k*g);
            end for;
        end for;
        S := T;
    end for;
    return S;
end function;

function SquareOrZero(a)
    F := Parent(a);
    if a eq F!0 then
        return true;
    end if;
    return a^((#F-1) div 2) eq F!1;
end function;

function EvalPolyQ(poly, val)
    F := Parent(val);
    out := F!0;
    for i in [0..Degree(poly)] do
        c := Coefficient(poly, i);
        ok, cr := RedRat(c, F);
        if not ok then
            return false, F!0;
        end if;
        out +:= cr*val^i;
    end for;
    return true, out;
end function;

function QuarticCoordsFromE(Qp)
    PS := Parent(Qp);
    F := Ring(PS);
    if Qp eq PS!0 then
        okm, mr := RedRat(mBase, F);
        okg, gr := RedRat(gBase, F);
        return okm and okg, mr, gr;
    end if;

    x := Qp[1]/Qp[3];
    y := Qp[2]/Qp[3];

    // Inverse of the minimal-model isomorphism, followed by the
    // inverse map from the original elliptic model to the quartic fiber.
    c1 := F!62500 / F!194481;
    c2 := F!2432977562500 / F!37822859361;
    c3 := F!15625000 / F!85766121;
    c4 := F!410439062500 / F!37822859361;
    c5 := F!4582629687500 / F!1400846643;

    x0 := c1*x - c2;
    y0 := c3*y + c4*x - c5;

    zc := 2*x0 - (F!2/F!25)*y0;
    if zc eq F!0 then
        return false, F!0, F!0;
    end if;

    mval := y0/zc;
    gnum := (F!3528/F!25)*x0^3
            - (F!5734246356100/F!85766121)*x0^2
            + (F!2101448/F!441)*x0*y0
            - (F!1764/F!25)*y0^2;
    gval := gnum/zc^2;

    return true, mval, gval;
end function;

function SBValueSquareAt(Qp)
    F := Ring(Parent(Qp));
    ok, mval, gval := QuarticCoordsFromE(Qp);
    if not ok then
        // Conservative: a bad map denominator should not eliminate a
        // rational class.
        return true;
    end if;

    den := mval^2 - (F!15225/F!32);
    if den eq F!0 then
        return true;
    end if;

    oka, av := EvalPolyQ(alphaP, mval);
    okh, hv := EvalPolyQ(hP, mval);
    if not (oka and okh) then
        return true;
    end if;

    // gval is the reduced fiber coordinate.  The unreduced square root
    // in alpha^2 - d/4 = G*h_B^2 is (609/256)*gval.
    val := 2*av + (F!609/F!256)*(2*hv/den)*gval;
    return SquareOrZero(val);
end function;

function ClassPoint(cls, gens, tors)
    P0 := Parent(gens[1])!0;
    if cls[1] eq 1 then P0 +:= tors[1]; end if;
    if cls[2] eq 1 then P0 +:= tors[2]; end if;
    P0 +:= cls[3]*gens[1] + cls[4]*gens[2] + cls[5]*gens[3];
    return P0;
end function;

TorsRange := (N mod 2 eq 0) select [0..1] else [0];
classes := [ <t1,t2,a,b,c> : t1 in TorsRange, t2 in TorsRange,
                            a in [0..N-1], b in [0..N-1], c in [0..N-1] ];

printf "R=-25/4 S_B MW sieve, N=%o, initial classes=%o\n", N, #classes;

for pp in PrimeList do
    F := GF(pp);
    try
        EF := EllipticCurve([F!0, F!-1, F!0, F!-468441, F!122191641]);
    catch e
        if not OnlyCuts then printf "p=%o skipped: singular model\n", pp; end if;
        continue;
    end try;
    if Discriminant(EF) eq F!0 then
        if not OnlyCuts then printf "p=%o skipped: bad elliptic reduction\n", pp; end if;
        continue;
    end if;

    redT := [];
    redF := [];
    good := true;
    for Tq in TgensQ do
        ok, Tp := RedPoint(Tq, EF);
        good and:= ok;
        Append(~redT, Tp);
    end for;
    for Gq in FgensQ do
        ok, Gp := RedPoint(Gq, EF);
        good and:= ok;
        Append(~redF, Gp);
    end for;
    if not good then
        if not OnlyCuts then printf "p=%o skipped: generator denominator\n", pp; end if;
        continue;
    end if;

    HN := AddGeneratedSubgroup([N*g : g in redF] cat [N*t : t in redT], EF);
    allowed := {Pp : Pp in Points(EF) | SBValueSquareAt(Pp)};

    allowedResidues := {};
    for A in allowed do
        for H in HN do
            Include(~allowedResidues, A - H);
        end for;
    end for;

    before := #classes;
    classes := [cls : cls in classes | ClassPoint(cls, redF, redT) in allowedResidues];
    if (not OnlyCuts) or (#classes lt before) then
        printf "p=%o #E(Fp)=%o allowed_points=%o #Nimage=%o classes %o -> %o\n",
            pp, #Points(EF), #allowed, #HN, before, #classes;
    end if;
    if (not OnlyCuts) and #classes le 20 then
        printf "surviving classes: %o\n", classes;
    end if;
    if #classes eq 0 then
        printf "SIEVE CERTIFICATE: no classes survive for N=%o\n", N;
        break;
    end if;
end for;

printf "FINAL surviving classes for N=%o: %o\n", N, #classes;
if #classes gt 0 and #classes le 200 then
    printf "%o\n", classes;
end if;

quit;
