//////////////////////////////////////////////////////////////////////
//  Focused p=3 exceptional charts for cyclic [49] in the contact-7
//  family.  This deliberately excludes the ordinary h(1)=0 residues
//  (a,b)=(0,1),(2,2), which have their own search.
//
//  Modes (default "all"):
//    q5            split ordinary-node branch (a,b)=(1,1) mod 3;
//                  search the first possible thickness v3(Q5)=7.
//    intersection  exceptional h(1)=Q5=0 residue (1,0); search the
//                  first contact-depth slice h(1)=3^7*s.
//    infinity      first three common parameter-pole layers
//                  a=u/3^e, b=v/3^e, e=1,2,3.
//    all           run all three after building away-from-3 masks once.
//
//  Every rational candidate is subjected to the exact necessary marked
//  condition D7_bar in 7*J(F_p), p=5,...,43.  Nonintegral or singular
//  displayed reductions are conservatively retained at that prime.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned mode then
    mode := "all";
end if;
if not assigned height then
    height := 12;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned prime_bound then
    prime_bound := 43;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;
if not assigned depth then
    depth := 7;
elif Type(depth) eq MonStgElt then
    depth := StringToInteger(depth);
end if;
if not assigned pole_depth then
    pole_depth := 3;
elif Type(pole_depth) eq MonStgElt then
    pole_depth := StringToInteger(pole_depth);
end if;
if not assigned max_exact then
    max_exact := 20;
elif Type(max_exact) eq MonStgElt then
    max_exact := StringToInteger(max_exact);
end if;
if not assigned max_record then
    max_record := 12;
elif Type(max_record) eq MonStgElt then
    max_record := StringToInteger(max_record);
end if;

Z := Integers();
Q := Rationals();
P<x> := PolynomialRing(Q);
filter_primes := [p : p in PrimesUpTo(prime_bound) | p ge 5 and p ne 7];

function Q5Component(a,b)
    return 432*a^4 - 64*a^3*b^3 + 1008*a^3*b + 3024*a^3
        - 448*a^2*b^3 + 224*a^2*b^2 + 21168*a^2*b - 32536*a^2
        - 2016*a*b^4 + 4480*a*b^3 + 38416*a*b^2 - 109760*a*b
        + 78890*a - 864*b^5 + 5936*b^4 + 7056*b^3
        - 96040*b^2 + 120050*b - 60025;
end function;

function Contact7Polynomial(a,b)
    h := 1 - (Q!7/2)*x + a*x^2 + b*x^3;
    num := h^2 + (x-1)^7;
    if Coefficient(num,0) ne 0 or Coefficient(num,1) ne 0 then
        return false,P!0,P!0;
    end if;
    return true,ExactQuotient(num,x^2),h;
end function;

function Contact7PolynomialFinite(F,a,b)
    PF<X> := PolynomialRing(F);
    h := 1 - (F!7/F!2)*X + a*X^2 + b*X^3;
    num := h^2 + (X-1)^7;
    if Coefficient(num,0) ne 0 or Coefficient(num,1) ne 0 then
        return false,PF!0,PF!0;
    end if;
    f := ExactQuotient(num,X^2);
    if Degree(f) ne 5 or Discriminant(f) eq 0 or Evaluate(h,F!1) eq 0 then
        return false,f,h;
    end if;
    return true,f,h;
end function;

function MarkedDivisibleBy7(D,G,phi)
    g := D @@ phi;
    coords := Eltseq(g);
    invs := Invariants(G);
    return &and[(Z!coords[i]) mod GCD(7,Z!invs[i]) eq 0 :
                i in [1..#coords]];
end function;

function BuildMarkedKillMask(p)
    F := GF(p);
    PF<X> := PolynomialRing(F);
    mask := [false : i in [1..p^2]];
    good := 0;
    pass := 0;
    for ai in [0..p-1] do
        for bi in [0..p-1] do
            ok,f,h := Contact7PolynomialFinite(F,F!ai,F!bi);
            if not ok then
                continue;
            end if;
            J := Jacobian(HyperellipticCurve(f));
            D7 := J![X-1,Evaluate(h,F!1)];
            if Order(D7) ne 7 then
                continue;
            end if;
            good +:= 1;
            G,phi := AbelianGroup(J);
            if MarkedDivisibleBy7(D7,G,phi) then
                pass +:= 1;
            else
                mask[ai*p+bi+1] := true;
            end if;
        end for;
    end for;
    return mask,good,pass;
end function;

function RationalParameters3Integral(B)
    vals := [];
    seen := {};
    for den in [1..B] do
        if den mod 3 eq 0 then
            continue;
        end if;
        for num in [-B..B] do
            if GCD(num,den) ne 1 then
                continue;
            end if;
            q := Q!num/Q!den;
            key := Sprint(q);
            if key notin seen then
                Include(~seen,key);
                Append(~vals,q);
            end if;
        end for;
    end for;
    return vals;
end function;

function V3(q)
    if q eq 0 then
        return 10^9;
    end if;
    return Valuation(Numerator(q),3)-Valuation(Denominator(q),3);
end function;

function RationalResidue(q,p)
    d := Denominator(q);
    if d mod p eq 0 then
        return false,0;
    end if;
    R := Integers(p);
    return true,Z!((R!Numerator(q))/(R!d));
end function;

function PassesAwayMasks(a,b,masks)
    used := 0;
    for j in [1..#filter_primes] do
        p := filter_primes[j];
        oka,ai := RationalResidue(a,p);
        okb,bi := RationalResidue(b,p);
        if not oka or not okb then
            continue;
        end if;
        used +:= 1;
        if masks[j][ai*p+bi+1] then
            return false,p,used;
        end if;
    end for;
    return true,0,used;
end function;

function ExactMarkedDivision(a,b)
    ok,f,h := Contact7Polynomial(a,b);
    if not ok or Degree(f) ne 5 or Discriminant(f) eq 0 then
        return false,false,[],"singular";
    end if;
    yp := Evaluate(h,Q!1);
    if yp eq 0 then
        return false,false,[],"marked boundary";
    end if;
    J := Jacobian(HyperellipticCurve(f));
    D7 := J![x-1,yp];
    if Order(D7) ne 7 then
        return false,false,[],"marked order not 7";
    end if;
    try
        isdiv,Q49 := IsDivisibleBy(D7,7);
        if isdiv then
            return true,true,[Order(Q49)],"";
        end if;
        return true,false,[],"";
    catch e
        return true,false,[],e`Object;
    end try;
end function;

// Lift the simple Q5 root a=1 mod 3 at fixed b=1 mod 3 to mod 3^N.
function HenselARootForB(b,N)
    modulus := 3^N;
    R := Integers(modulus);
    br := Z!((R!Numerator(b))/(R!Denominator(b)));
    ar := 1;
    for k in [1..N-1] do
        step := 3^k;
        nextmod := 3^(k+1);
        found := false;
        for t in [0..2] do
            cand := ar+t*step;
            if (Z!Q5Component(Z!cand,Z!br)) mod nextmod eq 0 then
                ar := cand;
                found := true;
                break;
            end if;
        end for;
        assert found;
    end for;
    ar mod:= modulus;
    if ar gt modulus div 2 then
        ar -:= modulus;
    end if;
    return ar;
end function;

procedure PrintKillCounts(label,kills)
    print label,[<filter_primes[j],kills[j]> : j in [1..#filter_primes]];
end procedure;

procedure RunQ5Chart(params,masks)
    modulus := 3^depth;
    checked := 0;
    depth7 := 0;
    smooth := 0;
    away_survivors := 0;
    exact_tests := 0;
    exact_hits := 0;
    kills := [0 : j in [1..#filter_primes]];
    records := [];

    brows := [b : b in params |
        (Integers(3)!Numerator(b))/(Integers(3)!Denominator(b)) eq 1];
    roots := [HenselARootForB(b,depth) : b in brows];
    for ib in [1..#brows] do
        b := brows[ib];
        a0 := roots[ib];
        for u in params do
            checked +:= 1;
            a := Q!a0 + modulus*u;
            q5 := Q5Component(a,b);
            vq := V3(q5);
            if vq ne depth then
                continue;
            end if;
            depth7 +:= 1;
            ok,f,h := Contact7Polynomial(a,b);
            if not ok or Degree(f) ne 5 or Discriminant(f) eq 0 then
                continue;
            end if;
            smooth +:= 1;
            pass,pbad,used := PassesAwayMasks(a,b,masks);
            if not pass then
                j := Index(filter_primes,pbad);
                kills[j] +:= 1;
                continue;
            end if;
            away_survivors +:= 1;
            if #records lt max_record then
                Append(~records,<a,b,vq>);
            end if;
            if exact_tests lt max_exact then
                valid,isdiv,extra,msg := ExactMarkedDivision(a,b);
                exact_tests +:= 1;
                if valid and isdiv then
                    exact_hits +:= 1;
                    print "Q5_EXACT_HIT",<a,b,extra>;
                end if;
            end if;
        end for;
    end for;
    print "Q5_CHART_SUMMARY","height",height,"brows",#brows,
          "urows",#params,"checked",checked,"vQ5_eq_depth",depth7,
          "smooth",smooth,"away_survivors",away_survivors,
          "exact_tests",exact_tests,"exact_hits",exact_hits;
    PrintKillCounts("Q5_FIRST_KILL",kills);
    print "Q5_SURVIVOR_RECORDS",records;
end procedure;

procedure RunIntersectionChart(params,masks)
    scale := 3^depth;
    checked := 0;
    smooth := 0;
    away_survivors := 0;
    exact_tests := 0;
    exact_hits := 0;
    kills := [0 : j in [1..#filter_primes]];
    q5vals := AssociativeArray();
    records := [];
    sunits := [s : s in params | V3(s) eq 0];

    for u in params do
        a := 1+3*u;
        for s in sunits do
            checked +:= 1;
            b := Q!5/2-a+scale*s;
            qv := V3(Q5Component(a,b));
            key := IntegerToString(qv);
            if IsDefined(q5vals,key) then
                q5vals[key] +:= 1;
            else
                q5vals[key] := 1;
            end if;
            ok,f,h := Contact7Polynomial(a,b);
            if not ok or Degree(f) ne 5 or Discriminant(f) eq 0 then
                continue;
            end if;
            smooth +:= 1;
            pass,pbad,used := PassesAwayMasks(a,b,masks);
            if not pass then
                j := Index(filter_primes,pbad);
                kills[j] +:= 1;
                continue;
            end if;
            away_survivors +:= 1;
            if #records lt max_record then
                Append(~records,<a,b,qv>);
            end if;
            if exact_tests lt max_exact then
                valid,isdiv,extra,msg := ExactMarkedDivision(a,b);
                exact_tests +:= 1;
                if valid and isdiv then
                    exact_hits +:= 1;
                    print "INTERSECTION_EXACT_HIT",<a,b,extra>;
                end if;
            end if;
        end for;
    end for;
    print "INTERSECTION_CHART_SUMMARY","height",height,"urows",#params,
          "sunit_rows",#sunits,"checked",checked,"smooth",smooth,
          "away_survivors",away_survivors,"exact_tests",exact_tests,
          "exact_hits",exact_hits;
    print "INTERSECTION_Q5_VALUATIONS",q5vals;
    PrintKillCounts("INTERSECTION_FIRST_KILL",kills);
    print "INTERSECTION_SURVIVOR_RECORDS",records;
end procedure;

procedure RunInfinityChart(params,masks)
    checked := 0;
    smooth := 0;
    away_survivors := 0;
    exact_tests := 0;
    exact_hits := 0;
    kills := [0 : j in [1..#filter_primes]];
    layer_counts := [0 : e in [1..pole_depth]];
    layer_survivors := [0 : e in [1..pole_depth]];
    records := [];
    for e in [1..pole_depth] do
        scale := 3^e;
        for u in params do
            for v in params do
                if Min(V3(u),V3(v)) ne 0 then
                    continue;
                end if;
                checked +:= 1;
                layer_counts[e] +:= 1;
                a := u/scale;
                b := v/scale;
                ok,f,h := Contact7Polynomial(a,b);
                if not ok or Degree(f) ne 5 or Discriminant(f) eq 0 then
                    continue;
                end if;
                smooth +:= 1;
                pass,pbad,used := PassesAwayMasks(a,b,masks);
                if not pass then
                    j := Index(filter_primes,pbad);
                    kills[j] +:= 1;
                    continue;
                end if;
                away_survivors +:= 1;
                layer_survivors[e] +:= 1;
                if #records lt max_record then
                    Append(~records,<e,a,b>);
                end if;
                if exact_tests lt max_exact then
                    valid,isdiv,extra,msg := ExactMarkedDivision(a,b);
                    exact_tests +:= 1;
                    if valid and isdiv then
                        exact_hits +:= 1;
                        print "INFINITY_EXACT_HIT",<e,a,b,extra>;
                    end if;
                end if;
            end for;
        end for;
    end for;
    print "INFINITY_CHART_SUMMARY","height",height,"pole_depth",pole_depth,
          "checked",checked,"smooth",smooth,"away_survivors",away_survivors,
          "exact_tests",exact_tests,"exact_hits",exact_hits;
    print "INFINITY_LAYER_CHECKED",layer_counts;
    print "INFINITY_LAYER_SURVIVORS",layer_survivors;
    PrintKillCounts("INFINITY_FIRST_KILL",kills);
    print "INFINITY_SURVIVOR_RECORDS",records;
end procedure;

print "Z49 EXCEPTIONAL p=3 CHARTS";
print "mode",mode,"height",height,"depth",depth,
      "pole_depth",pole_depth,"primes",filter_primes;
print "Q5 special fiber","f=x^2*(x^3+2*x+1), split node x=0";
print "Q5 normalization","E(F3)=Z/7, marked projection (1,1) has order 7";
print "Q5 necessary condition","7 divides node thickness v3(Q5); first depth=7";
print "intersection special fiber","f=x*(x-1)^4; non-ordinary four-root cluster";
print "intersection local h","h(1+z)=H/2+(H+b+3/2)z+(H/2+5/2+2b)z^2+bz^3";
print "infinity first scaled fiber","3^(2e)f -> x^2*(u+v*x)^2; highly degenerate";

masks := [];
for p in filter_primes do
    mask,good,pass := BuildMarkedKillMask(p);
    Append(~masks,mask);
    print "MASK",p,"good",good,"marked_div7",pass;
end for;
params := RationalParameters3Integral(height);
print "three_integral_parameters",#params;

if mode in {"q5","all"} then
    RunQ5Chart(params,masks);
end if;
if mode in {"intersection","all"} then
    RunIntersectionChart(params,masks);
end if;
if mode in {"infinity","all"} then
    RunInfinityChart(params,masks);
end if;

quit;
