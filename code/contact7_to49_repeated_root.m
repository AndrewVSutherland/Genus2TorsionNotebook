//////////////////////////////////////////////////////////////////////
// Repeated-root (2P-K) sublocus of the contact-7 to 49 norm cover.
//
// U=(x-p)^2 gives E=A^2-B^2*f+(x-1)*(x-p)^14=0.
// E14,...,E10 solve b4,...,b0 globally with coefficient -2.
// On a7 != 0, E9,E8,E7 then solve a2,a1,a0.  The residual
// chart has seven equations in alpha,beta,p,a3,...,a7, so its
// expected dimension is one.  No Groebner basis is used.
//
// Modes: symbolic (default), local, search, residual.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned mode then mode := "symbolic"; end if;
if not assigned primes then primes := "3,5,11,13,17,19,23,29"; end if;
if not assigned height then
    height := 4;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned max_print then
    max_print := 20;
elif Type(max_print) eq MonStgElt then
    max_print := StringToInteger(max_print);
end if;
if not assigned progress_interval then
    progress_interval := 50000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;

Q := Rationals(); Z := Integers();
PQ<xQ> := PolynomialRing(Q);
rr_input := [];
if assigned rr_alpha and assigned rr_beta and assigned rr_p and
        assigned rr_a3 and assigned rr_a4 and assigned rr_a5 and
        assigned rr_a6 and assigned rr_a7 then
    rr_input := [Q!StringToInteger(rr_alpha),Q!StringToInteger(rr_beta),Q!StringToInteger(rr_p),Q!StringToInteger(rr_a3),Q!StringToInteger(rr_a4),Q!StringToInteger(rr_a5),
                 Q!StringToInteger(rr_a6),Q!StringToInteger(rr_a7)];
end if;

function ParseIntegerList(S)
    if Type(S) eq SeqEnum then return [Z!z : z in S]; end if;
    return [StringToInteger(t) : t in Split(S, ",")];
end function;

function Contact7Data(k,av,bv)
    P<x> := PolynomialRing(k);
    h := 1-(k!7/k!2)*x+av*x^2+bv*x^3;
    num := h^2+(x-1)^7;
    if Coefficient(num,0) ne 0 or Coefficient(num,1) ne 0 then
        return false,P!0,P!0;
    end if;
    return true,ExactQuotient(num,x^2),h;
end function;

// Straight-line elimination evaluator. highA=[a3,a4,a5,a6,a7].
function ForcedRepeatedResiduals(k,av,bv,pv,highA)
    if #highA ne 5 then error "highA must be [a3,a4,a5,a6,a7]"; end if;
    P<x> := PolynomialRing(k);
    ok,f,h := Contact7Data(k,av,bv);
    if not ok or highA[5] eq 0 then
        return false,P!0,P!0,P!0,[],"bad contact data or a7=0";
    end if;
    A := &+[highA[i]*x^(i+2) : i in [1..5]];
    B := x^5; G := (x-1)*(x-pv)^14;
    for j in [4..0 by -1] do
        E := A^2-B^2*f+G;
        B +:= (Coefficient(E,j+10)/(k!2))*x^j;
        assert Coefficient(A^2-B^2*f+G,j+10) eq 0;
    end for;
    for j in [2..0 by -1] do
        E := A^2-B^2*f+G;
        A +:= (-Coefficient(E,j+7)/(k!2*highA[5]))*x^j;
        assert Coefficient(A^2-B^2*f+G,j+7) eq 0;
    end for;
    E := A^2-B^2*f+G;
    assert &and[Coefficient(E,j) eq 0 : j in [7..14]];
    return true,A,B,f,[Coefficient(E,j) : j in [0..6]],"ok";
end function;

procedure SymbolicReport()
    R<aa,bb,pp,a0,a1,a2,a3,a4,a5,a6,a7,b0,b1,b2,b3,b4> :=
        PolynomialRing(Q,16,"grevlex");
    P<x> := PolynomialRing(R);
    h := 1-(R!7/R!2)*x+aa*x^2+bb*x^3;
    f := ExactQuotient(h^2+(x-1)^7,x^2);
    A := a7*x^7+a6*x^6+a5*x^5+a4*x^4+a3*x^3+a2*x^2+a1*x+a0;
    B := x^5+b4*x^4+b3*x^3+b2*x^2+b1*x+b0;
    E := A^2-B^2*f+(x-1)*(x-pp)^14;
    assert Degree(E) le 14;
    print "# repeated-root compact cover: 16 variables, 15 equations";
    print "# expected dimension 1";
    print "TOP14",Coefficient(E,14);
    for j in [0..4] do
        assert Derivative(Coefficient(E,j+10),12+j) eq -2;
    end for;
    print "# E14..E10 solve b4..b0; E9..E7 solve a2..a0 if a7 != 0";
    print "# residual E0..E6 in alpha,beta,p,a3..a7";
    for i in [0..14] do
        ei := Coefficient(E,i);
        printf "E%o total_degree=%o terms=%o\n",i,TotalDegree(ei),#Terms(ei);
    end for;
end procedure;

function PairKey(a,b)
    return Sprintf("%o,%o",Z!a,Z!b);
end function;
function PointKey(a,b,pv,yv)
    return Sprintf("%o,%o,%o,%o",Z!a,Z!b,Z!pv,Z!yv);
end function;

function BuildLocalMask(ell,print_hits)
    F := GF(ell); P<x> := PolynomialRing(F);
    goodpairs := {}; allowed := {};
    smooth := 0; marked := 0; curves_hit := 0; shown := 0;
    for av in F do for bv in F do
        ok,f,h := Contact7Data(F,av,bv);
        if not ok or Degree(f) ne 5 or Discriminant(f) eq 0 or
                Evaluate(h,F!1) eq 0 then
            continue;
        end if;
        smooth +:= 1;
        J := Jacobian(HyperellipticCurve(f));
        D7 := J![x-1,Evaluate(h,F!1)];
        if Order(D7) ne 7 then continue; end if;
        marked +:= 1; Include(~goodpairs,PairKey(av,bv));
        this_hit := false;
        for pv in F do
            rhs := Evaluate(f,pv);
            if not IsSquare(rhs) then continue; end if;
            y0 := SquareRoot(rhs);
            ys := (y0 eq 0) select [F!0] else [y0,-y0];
            for yv in ys do
                T := J![x-pv,yv];
                if 14*T eq D7 or 14*T eq -D7 then
                    Include(~allowed,PointKey(av,bv,pv,yv));
                    this_hit := true;
                    if print_hits and shown lt max_print then
                        printf "LOCAL_POINT p=%o alpha=%o beta=%o xP=%o yP=%o #J=%o sign=%o\n",
                            ell,Z!av,Z!bv,Z!pv,Z!yv,#J,
                            (14*T eq D7) select 1 else -1;
                        shown +:= 1;
                    end if;
                end if;
            end for;
        end for;
        if this_hit then curves_hit +:= 1; end if;
    end for; end for;
    printf "LOCAL_SUMMARY p=%o smooth=%o marked7=%o good_pairs=%o curves_with_P=%o repeated_root_points=%o\n",
        ell,smooth,marked,#goodpairs,curves_hit,#allowed;
    return goodpairs,allowed;
end function;

procedure LocalScan()
    for ell in ParseIntegerList(primes) do
        if not IsPrime(ell) or ell in {2,7} then error "bad local prime"; end if;
        good,allowed := BuildLocalMask(ell,true);
    end for;
end procedure;

function RationalParametersOfHeight(B)
    vals := []; seen := {};
    for den in [1..B] do for num in [-B..B] do
        if GCD(num,den) ne 1 then continue; end if;
        q := Q!num/Q!den; key := Sprint(q);
        if key notin seen then Include(~seen,key); Append(~vals,q); end if;
    end for; end for;
    return vals;
end function;

function ReduceRationals(ell,vals)
    F := GF(ell); out := [];
    for q in vals do
        if Denominator(q) mod ell eq 0 then return false,[]; end if;
        Append(~out,F!q);
    end for;
    return true,out;
end function;

function PassMasks(av,bv,pv,yv,masks)
    used := 0;
    for M in masks do
        ell := M[1]; ok,r := ReduceRationals(ell,[av,bv,pv,yv]);
        if not ok then continue; end if;
        if PairKey(r[1],r[2]) notin M[2] then continue; end if;
        used +:= 1;
        if PointKey(r[1],r[2],r[3],r[4]) notin M[3] then
            return false,used,ell;
        end if;
    end for;
    return true,used,0;
end function;

function ExactRepeatedPoint(av,bv,pv,yv)
    ok,f,h := Contact7Data(Q,av,bv);
    if not ok or Degree(f) ne 5 or Discriminant(f) eq 0 then
        return false,false,0,0,"singular";
    end if;
    if yv^2 ne Evaluate(f,pv) then
        return true,false,0,0,"point identity failed";
    end if;
    L := LCM([Denominator(Coefficient(f,i)) : i in [0..Degree(f)]]);
    fI := PQ!(L^2*f);
    try
        J := Jacobian(HyperellipticCurve(fI));
        D7 := J![xQ-1,L*Evaluate(h,Q!1)];
        T := J![xQ-pv,L*yv];
        rel := (14*T eq D7) or (14*T eq -D7);
        return true,rel,Order(D7),rel select Order(2*T) else 0,"ok";
    catch e
        return true,false,0,0,Sprint(e`Object);
    end try;
end function;

procedure Increment(~A,key)
    if not IsDefined(A,key) then A[key] := 0; end if;
    A[key] +:= 1;
end procedure;

procedure RationalSearch()
    masks := [* *];
    print "# building exact repeated-root local masks";
    for ell in ParseIntegerList(primes) do
        if not IsPrime(ell) or ell in {2,7} then error "bad search prime"; end if;
        good,allowed := BuildLocalMask(ell,false);
        Append(~masks,<ell,good,allowed>);
    end for;
    vals := RationalParametersOfHeight(height);
    print "# rational parameter count",#vals,"height",height;
    checked := 0; checked_zero := 0; survivors := 0;
    exact_tests := 0; hits := 0; shown := 0;
    kills := AssociativeArray(); used_counts := AssociativeArray();

    // Generic p != 0,1 chart:
    // Y=(r+(p-1)^7/r)/2, H=((p-1)^7/r-r)/2.
    for pv in vals do
        if pv in {Q!0,Q!1} then continue; end if;
        c := (pv-1)^7;
        for rv in vals do
            if rv eq 0 then continue; end if;
            Y := (rv+c/rv)/2; H := (c/rv-rv)/2;
            for av in vals do
                checked +:= 1;
                if progress_interval gt 0 and checked mod progress_interval eq 0 then
                    print "PROGRESS",checked,"survivors",survivors,"hits",hits;
                end if;
                bv := (H-1+(Q!7/2)*pv-av*pv^2)/pv^3;
                yv := Y/pv;
                pass,used,kp := PassMasks(av,bv,pv,yv,masks);
                if not pass then Increment(~kills,kp); continue; end if;
                survivors +:= 1; Increment(~used_counts,used);
                smooth,rel,ordD,ordQ,msg := ExactRepeatedPoint(av,bv,pv,yv);
                exact_tests +:= 1;
                if rel then
                    hits +:= 1;
                    print "HIT49_GENERIC","alpha",av,"beta",bv,"xP",pv,"yP",yv,
                          "D7order",ordD,"Qorder",ordQ;
                elif shown lt max_print then
                    print "LOCAL_SURVIVOR_GENERIC","alpha",av,"beta",bv,"xP",pv,
                          "yP",yv,"used",used,"smooth",smooth,"message",msg;
                    shown +:= 1;
                end if;
            end for;
        end for;
    end for;

    // Separate p=0 chart: f(0)=2*alpha-35/4=y(P)^2.
    for yv in vals do
        av := yv^2/2+Q!35/8;
        for bv in vals do
            checked_zero +:= 1;
            pass,used,kp := PassMasks(av,bv,Q!0,yv,masks);
            if not pass then Increment(~kills,kp); continue; end if;
            survivors +:= 1; Increment(~used_counts,used);
            smooth,rel,ordD,ordQ,msg := ExactRepeatedPoint(av,bv,Q!0,yv);
            exact_tests +:= 1;
            if rel then
                hits +:= 1;
                print "HIT49_XZERO","alpha",av,"beta",bv,"yP",yv,
                      "D7order",ordD,"Qorder",ordQ;
            elif shown lt max_print then
                print "LOCAL_SURVIVOR_XZERO","alpha",av,"beta",bv,"yP",yv,
                      "used",used,"smooth",smooth,"message",msg;
                shown +:= 1;
            end if;
        end for;
    end for;

    print "SEARCH_DONE";
    print "generic_checked",checked; print "xzero_checked",checked_zero;
    print "local_survivors",survivors; print "exact_tests",exact_tests;
    print "hits49",hits; print "FIRST_KILL";
    for key in Sort([k : k in Keys(kills)]) do print key,kills[key]; end for;
    print "GOOD_PRIMES_USED_BY_SURVIVORS";
    for key in Sort([k : k in Keys(used_counts)]) do
        print key,used_counts[key];
    end for;
end procedure;

procedure ResidualEvaluation()
    if #rr_input ne 8 then
        error "residual mode needs rr_alpha,rr_beta,rr_p,rr_a3,rr_a4,rr_a5,rr_a6,rr_a7";
    end if;
    v := rr_input;
    ok,A,B,f,res,msg := ForcedRepeatedResiduals(Q,v[1],v[2],v[3],v[4..8]);
    print "OK",ok,"message",msg; print "A",A; print "B",B; print "f",f;
    print "residuals_E0_to_E6",res;
end procedure;

if mode eq "symbolic" then SymbolicReport();
elif mode eq "local" then LocalScan();
elif mode eq "search" then RationalSearch();
elif mode eq "residual" then ResidualEvaluation();
else error "unknown mode";
end if;
quit;
