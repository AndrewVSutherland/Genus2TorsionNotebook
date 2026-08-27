//////////////////////////////////////////////////////////////////////
//  Elkies/contact 5-torsion versus 4/6-torsion hybrid searches.
//
//  Modes:
//
//    finite_all
//        Finite-field existence/density checks for the three rational
//        4/8 or 6 bases and for contact-5 plus halving.
//
//    qsplit_search
//        Exact search in the split M_1(8) chart
//            y^2 = q(x)*(x^4 + q(x)), q = b*(x-b)*(x-c).
//
//    m14_search
//        Exact search on the M_1(8,4) tangent cover, filtered by
//        5-divisibility of reductions.
//
//    m2226_search
//        Exact search in the split M(2,2,2,6) family, filtered by
//        5-divisibility of reductions.
//
//    contact4_search
//        Exact search in the quintic-contact 5-torsion family, asking
//        whether the rational Weierstrass class at x=1 is 2-divisible.
//
//  Typical runs:
//      magma -b mode:="finite_all" code/elkies5_4_6_hybrid_search.m
//      magma -b mode:="qsplit_search" height:=40 code/elkies5_4_6_hybrid_search.m
//      magma -b mode:="m14_search" height:=18 code/elkies5_4_6_hybrid_search.m
//      magma -b mode:="m2226_search" height:=35 code/elkies5_4_6_hybrid_search.m
//      magma -b mode:="contact4_search" height:=30 code/elkies5_4_6_hybrid_search.m
//////////////////////////////////////////////////////////////////////

if not assigned mode then
    mode := "finite_all";
end if;
if not assigned height then
    height := 20;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned max_hits then
    max_hits := 20;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;
if not assigned max_tests then
    max_tests := 200;
elif Type(max_tests) eq MonStgElt then
    max_tests := StringToInteger(max_tests);
end if;

Q := Rationals();
Z := Integers();
Qx<x> := PolynomialRing(Q);

five_filter_primes := [3,7,11,13,17,19,23,29,31,37,41,43];
finite_primes := [3,7,11,13,17,19,23,29,31,37];
contact_halving_primes := [3,5,7,11,13,17,19,23];

function AbsGCD(vals)
    nz := [ Abs(Z!v) : v in vals | v ne 0 ];
    if #nz eq 0 then
        return 0;
    end if;
    return GCD(nz);
end function;

function IsPrimitivePair(a, b)
    return not (a eq 0 and b eq 0) and AbsGCD([a,b]) eq 1;
end function;

function IsPrimitiveTriple(a, b, c)
    return not (a eq 0 and b eq 0 and c eq 0) and AbsGCD([a,b,c]) eq 1;
end function;

function NormalizedTriple(a, b, c)
    if not IsPrimitiveTriple(a,b,c) then
        return false;
    end if;
    for z in [a,b,c] do
        if z ne 0 then
            return z gt 0;
        end if;
    end for;
    return false;
end function;

function RationalParametersOfHeight(B)
    vals := [];
    seen := {};
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num, den) ne 1 then
                continue;
            end if;
            r := Q!num/den;
            key := Sprint(r);
            if key notin seen then
                Include(~seen, key);
                Append(~vals, r);
            end if;
        end for;
    end for;
    return vals;
end function;

function IsSquareQ(q)
    q := Q!q;
    if q lt 0 then
        return false;
    end if;
    return IsSquare(Numerator(q)) and IsSquare(Denominator(q));
end function;

function IntegralModelPolynomial(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return Qx!(L^2*f), L;
end function;

function GoodHyperellipticPolynomial(f)
    return Degree(f) in {5,6} and Discriminant(f) ne 0;
end function;

function HasFiveInvariants(invs)
    return &or [ n mod 5 eq 0 : n in invs ];
end function;

function PassesFiveReduction(f, primes)
    for p in primes do
        try
            fp := ChangeRing(f, GF(p));
        catch e
            continue;
        end try;
        if not GoodHyperellipticPolynomial(fp) then
            continue;
        end if;
        C := HyperellipticCurve(fp);
        if (#Jacobian(C) mod 5) ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;

function QSplitM18Polynomial(K, b, c)
    P<X> := PolynomialRing(K);
    q := b*(X-b)*(X-c);
    return q*(X^4 + q);
end function;

function M14FamilyPolynomial(R, w)
    m := R;
    n := Q!1;
    t := (2*m^2 + (1-w^2)*m*n - 2*w^2*n^2)/(4*(w^2-1));
    A := n^4*x^2
         + (m^3*n + 4*m^2*t + m*n^3 - 8*m*n*t + 4*n^2*t)*x
         + m^4;
    B := (m*n + 2*n^2 + 4*t)*x^2
         + (m^2 + 4*m*n + n^2 + 8*t)*x
         + (2*m^2 + m*n + 4*t);
    return x*A*B, t, A, B;
end function;

function M14FamilyPolynomialFinite(K, R, w)
    P<X> := PolynomialRing(K);
    m := R;
    n := K!1;
    t := (2*m^2 + (1-w^2)*m*n - 2*w^2*n^2)/(4*(w^2-1));
    A := n^4*X^2
         + (m^3*n + 4*m^2*t + m*n^3 - 8*m*n*t + 4*n^2*t)*X
         + m^4;
    B := (m*n + 2*n^2 + 4*t)*X^2
         + (m^2 + 4*m*n + n^2 + 8*t)*X
         + (2*m^2 + m*n + 4*t);
    return X*A*B, t;
end function;

function PlusDisc(R, w)
    return -4*(w-R)*(R-1)^2*(R+1)*(R+w)*(w+1)
            *(R*w - 3*R + 3*w - 1);
end function;

function MinusDisc(R, w)
    return 4*(w-1)*(R+1)*(R*w + 3*R + 3*w + 1)
           *(R^4 - 2*R^3 + R^2*w^2 - R^2 + 2*R*w^2 - w^2);
end function;

function M14CoverPossible(R, w)
    return IsSquareQ(PlusDisc(R,w)) or IsSquareQ(MinusDisc(R,w));
end function;

function M2226Polynomial(K, s, m, n)
    P<X> := PolynomialRing(K);
    L := [
        K!1 + (2*s^2 - s*n)*X,
        K!1 + (2*s^2 + s*m - 2*s*n - m*n)*X,
        K!1 + (2*s^2 + s*m - s*n - m*n)*X,
        K!2 + (-m*n)*X,
        K!2 + (4*s^2 - 4*s*n - m*n)*X
    ];
    return &*L;
end function;

function ContactFivePolynomial(K, a, b)
    P<X> := PolynomialRing(K);
    h := 1 + a*X + b*X^2;
    k := (1+a+b)^2;
    return h^2 - k*X^5;
end function;

function ProjectiveTriples(K)
    pts := [];
    seen := {};
    for s in K do
        for m in K do
            for n in K do
                if s eq 0 and m eq 0 and n eq 0 then
                    continue;
                end if;
                if s ne 0 then
                    lam := s^-1;
                elif m ne 0 then
                    lam := m^-1;
                else
                    lam := n^-1;
                end if;
                vals := [s*lam, m*lam, n*lam];
                key := Sprint(vals);
                if key notin seen then
                    Include(~seen, key);
                    Append(~pts, vals);
                end if;
            end for;
        end for;
    end for;
    return pts;
end function;

procedure FiniteQSplit()
    print "FINITE qsplit M_1(8) chart: q=b*(x-b)*(x-c), need 5 | #J(Fp)";
    for p in finite_primes do
        K := GF(p);
        total := 0;
        good := 0;
        five := 0;
        samples := [];
        for b in K do
            for c in K do
                if b eq 0 then
                    continue;
                end if;
                total +:= 1;
                f := QSplitM18Polynomial(K, b, c);
                if not GoodHyperellipticPolynomial(f) then
                    continue;
                end if;
                good +:= 1;
                C := HyperellipticCurve(f);
                if (#Jacobian(C) mod 5) eq 0 then
                    five +:= 1;
                    if #samples lt 8 then
                        Append(~samples, <Z!b,Z!c>);
                    end if;
                end if;
            end for;
        end for;
        print "p", p, "total", total, "good", good, "five_possible", five,
              "samples", samples;
    end for;
end procedure;

procedure FiniteM14()
    print "FINITE M_1(8,4) tangent-cover chart: cover condition and 5 | #J(Fp)";
    print "  note: Magma IsDivisibleBy is Q-only, so finite verification is by cover equations";
    for p in finite_primes do
        K := GF(p);
        checked := 0;
        cover := 0;
        good := 0;
        five := 0;
        samples := [];
        for R in K do
            for w in K do
                if R eq 0 or w eq 0 or w eq 1 or w eq -1 then
                    continue;
                end if;
                checked +:= 1;
                if not (IsSquare(K!PlusDisc(R,w)) or IsSquare(K!MinusDisc(R,w))) then
                    continue;
                end if;
                cover +:= 1;
                f, t := M14FamilyPolynomialFinite(K, R, w);
                if not GoodHyperellipticPolynomial(f) then
                    continue;
                end if;
                good +:= 1;
                C := HyperellipticCurve(f);
                if (#Jacobian(C) mod 5) eq 0 then
                    five +:= 1;
                    if #samples lt 8 then
                        Append(~samples, <Z!R,Z!w>);
                    end if;
                end if;
            end for;
        end for;
        print "p", p, "checked", checked, "cover", cover, "good", good,
              "five_possible", five, "samples", samples;
    end for;
end procedure;

procedure FiniteM2226()
    print "FINITE split M(2,2,2,6): need 5 | #J(Fp)";
    for p in finite_primes do
        K := GF(p);
        pts := ProjectiveTriples(K);
        good := 0;
        five := 0;
        samples := [];
        for vals in pts do
            f := M2226Polynomial(K, vals[1], vals[2], vals[3]);
            if not GoodHyperellipticPolynomial(f) then
                continue;
            end if;
            good +:= 1;
            C := HyperellipticCurve(f);
            if (#Jacobian(C) mod 5) eq 0 then
                five +:= 1;
                if #samples lt 8 then
                    Append(~samples, <Z!vals[1], Z!vals[2], Z!vals[3]>);
                end if;
            end if;
        end for;
        print "p", p, "projective", #pts, "good", good,
              "five_possible", five, "samples", samples;
    end for;
end procedure;

procedure FiniteContactHalving()
    print "FINITE contact-5 family";
    print "  note: finite 2-divisibility is skipped because IsDivisibleBy is Q-only";
    for p in contact_halving_primes do
        K := GF(p);
        total := 0;
        good := 0;
        j5 := 0;
        samples := [];
        for a in K do
            for b in K do
                if 1+a+b eq 0 then
                    continue;
                end if;
                total +:= 1;
                f := ContactFivePolynomial(K, a, b);
                if not GoodHyperellipticPolynomial(f) then
                    continue;
                end if;
                good +:= 1;
                C := HyperellipticCurve(f);
                if (#Jacobian(C) mod 5) eq 0 then
                    j5 +:= 1;
                    if #samples lt 8 then
                        Append(~samples, <Z!a,Z!b>);
                    end if;
                end if;
            end for;
        end for;
        print "p", p, "total", total, "good", good, "j5", j5,
              "samples", samples;
    end for;
end procedure;

procedure SearchQSplit()
    print "SEARCH qsplit M_1(8) chart";
    print "height", height, "five_filter_primes", five_filter_primes;
    checked := 0;
    smooth := 0;
    survivors := 0;
    torsion_tests := 0;
    hits := [];
    stop := false;
    for bb in [1..height] do
        if stop then
            break;
        end if;
        for cc in [-height..height] do
            if not IsPrimitivePair(bb, cc) then
                continue;
            end if;
            checked +:= 1;
            f := QSplitM18Polynomial(Q, Q!bb, Q!cc);
            if not GoodHyperellipticPolynomial(f) then
                continue;
            end if;
            smooth +:= 1;
            if not PassesFiveReduction(f, five_filter_primes) then
                continue;
            end if;
            survivors +:= 1;
            fI, L := IntegralModelPolynomial(f);
            C := HyperellipticCurve(fI);
            J := Jacobian(C);
            torsion_tests +:= 1;
            G, phi := TorsionSubgroup(J);
            invs := Invariants(G);
            print "SURVIVOR", [bb,cc], "torsion", invs, "f", fI;
            if HasFiveInvariants(invs) then
                Append(~hits, <[bb,cc], invs, fI>);
                print "HIT", [bb,cc], invs;
                if #hits ge max_hits then
                    stop := true;
                    break;
                end if;
            end if;
            if torsion_tests ge max_tests then
                stop := true;
                break;
            end if;
        end for;
    end for;
    print "DONE qsplit height", height, "checked", checked, "smooth", smooth,
          "five_survivors", survivors, "torsion_tests", torsion_tests,
          "hits", #hits;
end procedure;

procedure SearchM14()
    print "SEARCH M_1(8,4) tangent cover plus 5";
    print "height", height, "five_filter_primes", five_filter_primes;
    params := RationalParametersOfHeight(height);
    checked := 0;
    cover := 0;
    smooth := 0;
    five_survivors := 0;
    verified := 0;
    torsion_tests := 0;
    hits := [];
    stop := false;
    for R in params do
        if stop then
            break;
        end if;
        for w in params do
            if R eq 0 or w in {Q!-1, Q!0, Q!1} then
                continue;
            end if;
            checked +:= 1;
            if not M14CoverPossible(R,w) then
                continue;
            end if;
            cover +:= 1;
            f, t, A, B := M14FamilyPolynomial(R,w);
            if not GoodHyperellipticPolynomial(f) then
                continue;
            end if;
            smooth +:= 1;
            if not PassesFiveReduction(f, five_filter_primes) then
                continue;
            end if;
            five_survivors +:= 1;
            fI, L := IntegralModelPolynomial(f);
            C := HyperellipticCurve(fI);
            J := Jacobian(C);
            D := J![x, Q!0];
            ok, half := IsDivisibleBy(D, 2);
            if not ok then
                continue;
            end if;
            verified +:= 1;
            torsion_tests +:= 1;
            G, phi := TorsionSubgroup(J);
            invs := Invariants(G);
            print "SURVIVOR", "R", R, "w", w, "t", t,
                  "half_order", Order(half), "torsion", invs;
            if HasFiveInvariants(invs) then
                Append(~hits, <R,w,t,invs,fI>);
                print "HIT", "R", R, "w", w, "torsion", invs, "f", fI;
                if #hits ge max_hits then
                    stop := true;
                    break;
                end if;
            end if;
            if torsion_tests ge max_tests then
                stop := true;
                break;
            end if;
        end for;
    end for;
    print "DONE m14 height", height, "checked", checked, "cover", cover,
          "smooth", smooth, "five_survivors", five_survivors,
          "verified", verified, "torsion_tests", torsion_tests, "hits", #hits;
end procedure;

procedure SearchM2226()
    print "SEARCH split M(2,2,2,6) plus 5";
    print "height", height, "five_filter_primes", five_filter_primes;
    checked := 0;
    smooth := 0;
    five_survivors := 0;
    torsion_tests := 0;
    hits := [];
    stop := false;
    for ss in [-height..height] do
        if stop then
            break;
        end if;
        for mm in [-height..height] do
            if stop then
                break;
            end if;
            for nn in [-height..height] do
                if not NormalizedTriple(ss,mm,nn) then
                    continue;
                end if;
                checked +:= 1;
                f := M2226Polynomial(Q, Q!ss, Q!mm, Q!nn);
                if not GoodHyperellipticPolynomial(f) then
                    continue;
                end if;
                smooth +:= 1;
                if not PassesFiveReduction(f, five_filter_primes) then
                    continue;
                end if;
                five_survivors +:= 1;
                fI, L := IntegralModelPolynomial(f);
                C := HyperellipticCurve(fI);
                J := Jacobian(C);
                torsion_tests +:= 1;
                G, phi := TorsionSubgroup(J);
                invs := Invariants(G);
                print "SURVIVOR", [ss,mm,nn], "torsion", invs, "f", fI;
                if HasFiveInvariants(invs) then
                    Append(~hits, <[ss,mm,nn], invs, fI>);
                    print "HIT", [ss,mm,nn], invs;
                    if #hits ge max_hits then
                        stop := true;
                        break;
                    end if;
                end if;
                if torsion_tests ge max_tests then
                    stop := true;
                    break;
                end if;
            end for;
        end for;
    end for;
    print "DONE m2226 height", height, "checked", checked, "smooth", smooth,
          "five_survivors", five_survivors, "torsion_tests", torsion_tests,
          "hits", #hits;
end procedure;

procedure SearchContact4()
    print "SEARCH contact-5 family with halved visible 2-torsion at x=1";
    print "height", height;
    checked := 0;
    smooth := 0;
    halved := 0;
    torsion_tests := 0;
    hits := [];
    stop := false;
    for aa in [-height..height] do
        if stop then
            break;
        end if;
        for bb in [-height..height] do
            if 1+aa+bb eq 0 then
                continue;
            end if;
            checked +:= 1;
            f := ContactFivePolynomial(Q, Q!aa, Q!bb);
            if not GoodHyperellipticPolynomial(f) then
                continue;
            end if;
            smooth +:= 1;
            fI, L := IntegralModelPolynomial(f);
            C := HyperellipticCurve(fI);
            J := Jacobian(C);
            D := J![x - Q!1, Q!0];
            ok, half := IsDivisibleBy(D, 2);
            if not ok then
                continue;
            end if;
            halved +:= 1;
            torsion_tests +:= 1;
            G, phi := TorsionSubgroup(J);
            invs := Invariants(G);
            print "HALVED", [aa,bb], "half_order", Order(half),
                  "torsion", invs, "f", fI;
            if #invs gt 0 and &*[ n : n in invs ] mod 20 eq 0 then
                Append(~hits, <[aa,bb], invs, fI>);
                print "HIT", [aa,bb], invs;
                if #hits ge max_hits then
                    stop := true;
                    break;
                end if;
            end if;
            if torsion_tests ge max_tests then
                stop := true;
                break;
            end if;
        end for;
    end for;
    print "DONE contact4 height", height, "checked", checked, "smooth", smooth,
          "halved", halved, "torsion_tests", torsion_tests, "hits", #hits;
end procedure;

if mode eq "finite_all" then
    FiniteQSplit();
    FiniteM14();
    FiniteM2226();
    FiniteContactHalving();
elif mode eq "qsplit_search" then
    SearchQSplit();
elif mode eq "m14_search" then
    SearchM14();
elif mode eq "m2226_search" then
    SearchM2226();
elif mode eq "contact4_search" then
    SearchContact4();
else
    error "unknown mode";
end if;

quit;
