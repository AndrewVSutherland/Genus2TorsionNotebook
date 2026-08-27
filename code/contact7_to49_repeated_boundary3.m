//////////////////////////////////////////////////////////////////////
// Repeated-root Q=2P-K search on the forced p=3 boundary.
//
// The relation tested away from 3 is the exact marked-point condition
//     14*[P-infinity] = +/-D7,
// not merely D7 in 7*J(F_p).  The first feasible charts are:
//   ordinary:    h(1)=3^7*s at residues (0,1),(2,2), v3(s)=0;
//   q5:          v3(Q5(a,b))=7 at residue (1,1);
//   intersection h(1)=3^7*s at the four-root residue (1,0);
//   poles:       a=u/3^e,b=v/3^e (bounded unresolved probe).
// Integral point abscissas and the first x(P)-pole layer are tested.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned mode then mode := "all"; end if;
if not assigned height then height := 6;
elif Type(height) eq MonStgElt then height := StringToInteger(height); end if;
if not assigned prime_bound then prime_bound := 43;
elif Type(prime_bound) eq MonStgElt then prime_bound := StringToInteger(prime_bound); end if;
if not assigned depth then depth := 7;
elif Type(depth) eq MonStgElt then depth := StringToInteger(depth); end if;
if not assigned pole_depth then pole_depth := 1;
elif Type(pole_depth) eq MonStgElt then pole_depth := StringToInteger(pole_depth); end if;
if not assigned point_pole_depth then point_pole_depth := 1;
elif Type(point_pole_depth) eq MonStgElt then point_pole_depth := StringToInteger(point_pole_depth); end if;
if not assigned max_record then max_record := 12;
elif Type(max_record) eq MonStgElt then max_record := StringToInteger(max_record); end if;

Z := Integers(); Q := Rationals();
P<x> := PolynomialRing(Q);
filter_primes := [p : p in PrimesUpTo(prime_bound) | p ge 5 and p ne 7];

function V3(q)
    if q eq 0 then return 10^9; end if;
    return Valuation(Numerator(q),3)-Valuation(Denominator(q),3);
end function;

function Parameters3Integral(B)
    vals := []; seen := {};
    for d in [1..B] do
        if d mod 3 eq 0 then continue; end if;
        for n in [-B..B] do
            if GCD(n,d) ne 1 then continue; end if;
            q := Q!n/Q!d; key := Sprint(q);
            if key notin seen then Include(~seen,key); Append(~vals,q); end if;
        end for;
    end for;
    return vals;
end function;

function ContactData(k,a,b)
    R<X> := PolynomialRing(k);
    h := 1-(k!7/k!2)*X+a*X^2+b*X^3;
    num := h^2+(X-1)^7;
    if Coefficient(num,0) ne 0 or Coefficient(num,1) ne 0 then
        return false,R!0,R!0;
    end if;
    return true,ExactQuotient(num,X^2),h;
end function;

function Q5Component(a,b)
    return 432*a^4-64*a^3*b^3+1008*a^3*b+3024*a^3
      -448*a^2*b^3+224*a^2*b^2+21168*a^2*b-32536*a^2
      -2016*a*b^4+4480*a*b^3+38416*a*b^2-109760*a*b
      +78890*a-864*b^5+5936*b^4+7056*b^3-96040*b^2
      +120050*b-60025;
end function;

function PairKey(a,b)
    return Sprintf("%o,%o",Z!a,Z!b);
end function;
function PointKey(a,b,pv,yv)
    return Sprintf("%o,%o,%o,%o",Z!a,Z!b,Z!pv,Z!yv);
end function;

function BuildPointMask(ell)
    F := GF(ell); R<X> := PolynomialRing(F);
    good := {}; allowed := {}; curves_hit := 0;
    for a in F do for b in F do
        ok,f,h := ContactData(F,a,b);
        if not ok or Degree(f) ne 5 or Discriminant(f) eq 0 or Evaluate(h,F!1) eq 0 then
            continue;
        end if;
        J := Jacobian(HyperellipticCurve(f));
        D := J![X-1,Evaluate(h,F!1)];
        if Order(D) ne 7 then continue; end if;
        Include(~good,PairKey(a,b)); hit := false;
        for pv in F do
            rhs := Evaluate(f,pv);
            if not IsSquare(rhs) then continue; end if;
            y0 := SquareRoot(rhs);
            ys := (y0 eq 0) select [F!0] else [y0,-y0];
            for yv in ys do
                T := J![X-pv,yv];
                if 14*T eq D or 14*T eq -D then
                    Include(~allowed,PointKey(a,b,pv,yv)); hit := true;
                end if;
            end for;
        end for;
        if hit then curves_hit +:= 1; end if;
    end for; end for;
    printf "MASK p=%o good_curves=%o curves_with_P=%o points=%o\n",
        ell,#good,curves_hit,#allowed;
    return good,allowed;
end function;

function RationalResidueTuple(ell,vals)
    F := GF(ell); ans := [];
    for q in vals do
        if Denominator(q) mod ell eq 0 then return false,[]; end if;
        Append(~ans,F!q);
    end for;
    return true,ans;
end function;

function RationalSquareRoot(q)
    if q lt 0 then return false,Q!0; end if;
    n := Numerator(q); d := Denominator(q);
    if not IsSquare(n) or not IsSquare(d) then return false,Q!0; end if;
    return true,Q!Isqrt(n)/Q!Isqrt(d);
end function;

function PassPointMasks(a,b,pv,yv,masks)
    used := 0;
    for j in [1..#masks] do
        ell := filter_primes[j];
        ok,r := RationalResidueTuple(ell,[a,b,pv,yv]);
        if not ok then continue; end if;
        if PairKey(r[1],r[2]) notin masks[j][1] then continue; end if;
        used +:= 1;
        if PointKey(r[1],r[2],r[3],r[4]) notin masks[j][2] then
            return false,ell,used;
        end if;
    end for;
    return true,0,used;
end function;

function ExactRelation(f,h,pv,yv)
    L := LCM([Denominator(Coefficient(f,i)) : i in [0..Degree(f)]]);
    fI := P!(L^2*f);
    try
        J := Jacobian(HyperellipticCurve(fI));
        D := J![x-1,L*Evaluate(h,Q!1)];
        T := J![x-pv,L*yv];
        rel := 14*T eq D or 14*T eq -D;
        return rel,Order(D),rel select Order(2*T) else 0,"";
    catch e
        return false,0,0,Sprint(e`Object);
    end try;
end function;

procedure ScanCurve(label,a,b,xvals,masks,~stats,~kills,~records)
    stats[1] +:= 1;
    ok,f,h := ContactData(Q,a,b);
    if not ok or Degree(f) ne 5 or Discriminant(f) eq 0 then return; end if;
    stats[2] +:= 1;
    for pv in xvals do
        stats[3] +:= 1;
        rhs := Evaluate(f,pv);
        square,y0 := RationalSquareRoot(rhs);
        if not square then continue; end if;
        ys := (y0 eq 0) select [Q!0] else [y0,-y0];
        for yv in ys do
            stats[4] +:= 1;
            pass,pbad,used := PassPointMasks(a,b,pv,yv,masks);
            if not pass then
                kills[Index(filter_primes,pbad)] +:= 1;
                continue;
            end if;
            stats[5] +:= 1; stats[6] +:= 1;
            rel,ordD,ordQ,msg := ExactRelation(f,h,pv,yv);
            if rel then
                stats[7] +:= 1;
                print "BOUNDARY_HIT49",label,<a,b,pv,yv,ordD,ordQ>;
            elif #records lt max_record then
                Append(~records,<label,a,b,pv,yv,used,msg>);
            end if;
        end for;
    end for;
end procedure;

procedure PrintSummary(label,stats,kills,records)
    print label,"curves",stats[1],"smooth",stats[2],"x_checks",stats[3],
          "rational_points",stats[4],"mask_survivors",stats[5],
          "exact_tests",stats[6],"hits49",stats[7];
    print label cat "_FIRST_KILL",
        [<filter_primes[j],kills[j]> : j in [1..#filter_primes] | kills[j] gt 0];
    print label cat "_SURVIVORS",records;
end procedure;

function HenselARootForB(b,N)
    modulus := 3^N; R := Integers(modulus);
    br := Z!((R!Numerator(b))/(R!Denominator(b))); ar := 1;
    for k in [1..N-1] do
        step := 3^k; nextmod := 3^(k+1); found := false;
        for t in [0..2] do
            cand := ar+t*step;
            if (Z!Q5Component(cand,br)) mod nextmod eq 0 then
                ar := cand; found := true; break;
            end if;
        end for;
        assert found;
    end for;
    ar mod:= modulus;
    if ar gt modulus div 2 then ar -:= modulus; end if;
    return ar;
end function;

procedure RunOrdinary(params,xvals,masks)
    sunits := [s : s in params | V3(s) eq 0];
    for a0 in [0,2] do
        stats := [0 : i in [1..7]]; kills := [0 : p in filter_primes]; records := [];
        for u in params do for s in sunits do
            a := Q!a0+3*u; b := Q!5/2-a+3^depth*s;
            ScanCurve("ordinary_" cat IntegerToString(a0),a,b,xvals,masks,
                      ~stats,~kills,~records);
        end for; end for;
        PrintSummary("ORDINARY_" cat IntegerToString(a0),stats,kills,records);
    end for;
end procedure;

procedure RunQ5(params,xvals,masks)
    stats := [0 : i in [1..7]]; kills := [0 : p in filter_primes]; records := [];
    brows := [b : b in params | V3(b-1) gt 0];
    roots := [HenselARootForB(b,depth) : b in brows];
    for i in [1..#brows] do for u in params do
        b := brows[i]; a := Q!roots[i]+3^depth*u;
        if V3(Q5Component(a,b)) ne depth then continue; end if;
        ScanCurve("q5",a,b,xvals,masks,~stats,~kills,~records);
    end for; end for;
    PrintSummary("Q5",stats,kills,records);
end procedure;

procedure RunIntersection(params,xvals,masks)
    stats := [0 : i in [1..7]]; kills := [0 : p in filter_primes]; records := [];
    sunits := [s : s in params | V3(s) eq 0];
    for u in params do for s in sunits do
        a := 1+3*u; b := Q!5/2-a+3^depth*s;
        ScanCurve("intersection",a,b,xvals,masks,~stats,~kills,~records);
    end for; end for;
    PrintSummary("INTERSECTION",stats,kills,records);
end procedure;

procedure RunPoles(params,xvals,masks)
    for e in [1..pole_depth] do
        stats := [0 : i in [1..7]]; kills := [0 : p in filter_primes]; records := [];
        for u in params do for v in params do
            if Minimum(V3(u),V3(v)) ne 0 then continue; end if;
            a := u/3^e; b := v/3^e;
            ScanCurve("poles_" cat IntegerToString(e),a,b,xvals,masks,
                      ~stats,~kills,~records);
        end for; end for;
        PrintSummary("POLES_" cat IntegerToString(e),stats,kills,records);
    end for;
end procedure;

params := Parameters3Integral(height);
xvals := params; seenx := {Sprint(t) : t in xvals};
for e in [1..point_pole_depth] do
    for w in params do
        if V3(w) ne 0 then continue; end if;
        t := w/3^e; key := Sprint(t);
        if key notin seenx then Include(~seenx,key); Append(~xvals,t); end if;
    end for;
end for;
print "BOUNDARY_REPEATED_ROOT","height",height,"params",#params,
      "xvals",#xvals,"depth",depth,"point_pole_depth",point_pole_depth;
masks := [* *];
for ell in filter_primes do
    good,allowed := BuildPointMask(ell); Append(~masks,<good,allowed>);
end for;

if mode in {"ordinary","all"} then RunOrdinary(params,xvals,masks); end if;
if mode in {"q5","all"} then RunQ5(params,xvals,masks); end if;
if mode in {"intersection","all"} then RunIntersection(params,xvals,masks); end if;
if mode in {"poles","all"} then RunPoles(params,xvals,masks); end if;
quit;
