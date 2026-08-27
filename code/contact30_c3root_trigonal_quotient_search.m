//////////////////////////////////////////////////////////////////////
// Reproducible quotient and bounded rational searches for the
// degree-3 C3-root cover of the simultaneous contact-(5,6) family.
//
// The companion geometry script proves that the original trigonal
// curve has genus 12.  Here we:
//   * certify its involution and derive the exact genus-6 quotient;
//   * sieve R=a/b and test the surviving cubics exactly over Q;
//   * sieve the quotient coordinate v and test which quotient points
//     lift through 9v^2-60v+84 = square.
//
// Typical runs (all default searches use primes 7,...,83):
//   magma -b mode:=quotient eps:=-1 \
//       code/contact30_c3root_trigonal_quotient_search.m
//   magma -b mode:=search eps:=-1 height:=5000 \
//       code/contact30_c3root_trigonal_quotient_search.m
//   magma -b mode:=quotient-search eps:=-1 height:=5000 \
//       code/contact30_c3root_trigonal_quotient_search.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);
SetMemoryLimit(4*10^9);

if not assigned mode then mode := "quotient"; end if;
if not assigned eps then eps := -1;
elif Type(eps) eq MonStgElt then eps := StringToInteger(eps);
end if;
assert eps in {-1,1};
if not assigned height then height := 5000;
elif Type(height) eq MonStgElt then height := StringToInteger(height);
end if;
if assigned primes then
    if Type(primes) eq MonStgElt then
        prime_list := [ StringToInteger(a) : a in Split(primes,",")
                        | #a gt 0 ];
    else
        prime_list := primes;
    end if;
else
    prime_list := [7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,
                   67,71,73,79,83];
end if;
assert mode in {"quotient","search","quotient-search"};
assert height ge 1;
assert &and[ IsPrime(p) and p gt 5 : p in prime_list ];

function PrimitiveCubic(k,branch)
    kR<R> := PolynomialRing(k);
    K := FieldOfFractions(kR);
    KRho<rho> := PolynomialRing(K);

    t := (5*R^2-20*R+19)/(R^2-5);
    Y := -2*(5*R^2-22*R+25)/(R^2-5);
    u := t^3;
    s := t^5+t^4+(k!5/2)*t^3+(k!1/2)*t
       + branch*t*(t-k!1/2)*(t+1)*Y;
    C := (u^2+1)/(2*u);
    D := u^6+6*u^4*s-2*u^4+15*u^3*s-u*s^3+u^2;
    N := 15*u^5+90*u^4+20*u^3*s-6*u^2*s^2+231*u^3
       + 2*u^2*s-15*u*s^2+90*u^2-20*u*s+15*u-2*s;
    q := N/D;
    A := (s+q)/2;
    B := (15-s*q)/2;
    f := 2*rho^3+(A-3)*rho^2+(B+3)*rho+(C-1);

    den := &*[ Denominator(Coefficient(f,i)) : i in [0..3] ];
    cs := [ kR!(den*Coefficient(f,i)) : i in [0..3] ];
    cont := cs[1];
    for i in [2..4] do cont := GCD(cont,cs[i]); end for;
    cs := [ ExactQuotient(c,cont) : c in cs ];
    if k cmpeq Rationals() then
        zcont := GCD([ GCD([ Abs(Numerator(cc))
                            : cc in Coefficients(c) ])
                       : c in cs | c ne 0 ]);
        if zcont gt 1 then cs := [ c/zcont : c in cs ]; end if;
    end if;
    kRRho<z> := PolynomialRing(kR);
    F := &+[ kRRho!cs[i+1]*z^i : i in [0..3] ];
    return F,R,t,Y,u,s,D,N;
end function;

function PrimitiveIntegralNested(f)
    ZZ := Integers(); QQ := Rationals();
    cs := [];
    for i in [0..Degree(f)] do
        c := Coefficient(f,i);
        if c ne 0 then
            cs cat:= [ Coefficient(c,j) : j in [0..Degree(c)] ];
        end if;
    end for;
    den := LCM([ Denominator(c) : c in cs ]);
    nums := [ ZZ!(den*c) : c in cs | c ne 0 ];
    cont := GCD([ Abs(n) : n in nums ]);
    scale := QQ!den/QQ!cont;
    B := BaseRing(Parent(f));
    return Parent(f)![ B!(scale*Coefficient(f,i))
                        : i in [0..Degree(f)] ],scale;
end function;

function QuotientPolynomial(F)
    QQ := Rationals();
    A<r0,v0,z0> := PolynomialRing(QQ,3);
    fm := A!0;
    for i in [0..Degree(F)] do
        c := Coefficient(F,i);
        for j in [0..Degree(c)] do
            fm +:= QQ!Coefficient(c,j)*r0^j*z0^i;
        end for;
    end for;
    rel := 3*r0^2-3*v0*r0+5*v0-7;
    res := Resultant(fm,rel,r0);
    fs := Factorization(res);
    cand := [ fe[1] : fe in fs |
              fe[2] eq 2 and Degree(fe[1],v0) eq 10
              and Degree(fe[1],z0) eq 3 ];
    assert #cand eq 1;
    assert Degree(cand[1],r0) eq 0;
    scalar := ExactQuotient(res,cand[1]^2);
    assert TotalDegree(scalar) eq 0;

    B<v,z> := PolynomialRing(QQ,2);
    mp := hom<A -> B | 0,v,z>;
    qmv := mp(cand[1]);
    Qv<V> := PolynomialRing(QQ);
    QvZ<Z> := PolynomialRing(Qv);
    q := QvZ!0;
    mons := Monomials(qmv);
    coeffs := Coefficients(qmv);
    for i in [1..#mons] do
        ee := Exponents(mons[i]);
        q +:= Qv!(coeffs[i])*V^ee[1]*Z^ee[2];
    end for;
    qI,scale := PrimitiveIntegralNested(q);
    term_count := &+[ #Terms(Coefficient(q,i)) : i in [0..Degree(q)] ];
    return q,qI,scale,term_count;
end function;

function MobiusTransform(h,a,b,c,d)
    P := Parent(h); X := P.1;
    K := FieldOfFractions(P);
    y := (a*X+b)/(c*X+d);
    out := P!Numerator(K!Evaluate(h,y));
    return out/LeadingCoefficient(out);
end function;

function BuildAllowedMask(f,p)
    // Degree drops are kept: they include affine roots reducing to
    // z=infinity.  Code p denotes the parameter point at infinity.
    k := GF(p);
    kT<T> := PolynomialRing(k);
    kz<z> := PolynomialRing(k);
    coeffs := [];
    for i in [0..Degree(f)] do
        c := Coefficient(f,i);
        if c eq 0 then
            Append(~coeffs,kT!0);
        else
            Append(~coeffs,kT![ k!Coefficient(c,j)
                                : j in [0..Degree(c)] ]);
        end if;
    end for;
    allowed := { Integers() | };
    for tv in k do
        spec := kz![ Evaluate(c,tv) : c in coeffs ];
        if Degree(spec) lt Degree(f) or #Roots(spec) gt 0 then
            Include(~allowed,Integers()!tv);
        end if;
    end for;
    d := Max([ Degree(c) : c in coeffs ]);
    inf := kz![ Coefficient(c,d) : c in coeffs ];
    if Degree(inf) lt Degree(f) or #Roots(inf) gt 0 then
        Include(~allowed,Integers()!p);
    end if;
    return allowed;
end function;

function SearchTrigonal(fQ,fI,height,prime_list,label)
    QQ := Rationals(); ZZ := Integers();
    Qz<z> := PolynomialRing(QQ);
    masks := [ BuildAllowedMask(fI,p) : p in prime_list ];
    tested := 0; survivors := 0; hits := [];
    for b in [1..height] do
        for a in [-height..height] do
            if GCD(Abs(a),b) ne 1 then continue; end if;
            tested +:= 1;
            passed := true;
            for j in [1..#prime_list] do
                p := prime_list[j];
                code := (b mod p eq 0) select p
                        else ZZ!((GF(p)!a)/(GF(p)!b));
                if not code in masks[j] then
                    passed := false;
                    break;
                end if;
            end for;
            if not passed then continue; end if;
            survivors +:= 1;
            tv := QQ!a/QQ!b;
            spec := Qz![ Evaluate(Coefficient(fQ,i),tv)
                         : i in [0..Degree(fQ)] ];
            for rt in Roots(spec) do
                Append(~hits,<tv,rt[1]>);
            end for;
        end for;
    end for;
    d := Max([ Degree(Coefficient(fQ,i)) : i in [0..Degree(fQ)] ]);
    inf := Qz![ Coefficient(Coefficient(fQ,i),d)
                : i in [0..Degree(fQ)] ];
    print label,"HEIGHT",height,"PRIMES",prime_list;
    print "TESTED",tested,"SURVIVORS",survivors,"HITS",hits;
    print "INFINITY_CUBIC",inf,"INFINITY_ROOTS",Roots(inf);
    return hits;
end function;

function IsSquareQ(a)
    if a lt 0 then return false; end if;
    return IsSquare(Integers()!Numerator(a))
           and IsSquare(Integers()!Denominator(a));
end function;

function IsOpenR(rv,R,t,D)
    rden := R^2-5;
    tnum := 5*R^2-20*R+19;
    tminus := Numerator(t-1);
    tplus := Numerator(t+1);
    Dnum := Numerator(D);
    return Evaluate(rden,rv) ne 0 and Evaluate(tnum,rv) ne 0
           and Evaluate(tminus,rv) ne 0 and Evaluate(tplus,rv) ne 0
           and Evaluate(Dnum,rv) ne 0;
end function;

QQ := Rationals();
F,R,t,Y,u,s,D,N := PrimitiveCubic(QQ,eps);
FI,Fscale := PrimitiveIntegralNested(F);

if mode eq "search" then
    print "CONTACT30_C3ROOT_TRIGONAL_SEARCH","eps",eps;
    hits := SearchTrigonal(F,FI,height,prime_list,"R_SEARCH");
    print "CONTACT30_C3ROOT_TRIGONAL_SEARCH_DONE";
    quit;
end if;

q,qI,qscale,qterms := QuotientPolynomial(F);

if mode eq "quotient" then
    QR := Parent(R);
    K := FieldOfFractions(QR);
    Kz<z> := PolynomialRing(K);
    sigma := (5*R-7)/(3*R-5);
    vfun := (3*R^2-7)/(3*R-5);
    FK := Kz![ K!Coefficient(F,i) : i in [0..Degree(F)] ];
    Fsigma := Kz![ Evaluate(Coefficient(F,i),sigma)
                    : i in [0..Degree(F)] ];

    print "CONTACT30_C3ROOT_TRIGONAL_QUOTIENT","eps",eps;
    print "SIGMA_SQUARE",Evaluate(sigma,sigma) eq R;
    print "SIGMA_V_FIXED",Evaluate(vfun,sigma) eq vfun;
    print "SIGMA_RHO_FIXED_MONIC_CUBIC",
          Fsigma/LeadingCoefficient(Fsigma) eq FK/LeadingCoefficient(FK);
    max_v_degree := Max([ Degree(Coefficient(q,i))
                          : i in [0..Degree(q)] ]);
    print "QUOTIENT_BIDEGREE",<max_v_degree,Degree(q)>,
          "MAX_V_DEGREE",max_v_degree,
          "TERMS",qterms,"IRREDUCIBLE",IsIrreducible(q);
    print "QUOTIENT_PRIMITIVE_INTEGRAL_SCALE",qscale;

    Qv := BaseRing(Parent(q));
    discq := Discriminant(q);
    print "QUOTIENT_DISCRIMINANT_DEGREE",Degree(discq);
    print "QUOTIENT_DISCRIMINANT_FACTORIZATION",Factorization(discq);

    Kv<T> := FunctionField(QQ);
    Kvz<w> := PolynomialRing(Kv);
    qK := Kvz![ Kv!Coefficient(q,i) : i in [0..Degree(q)] ];
    assert IsIrreducible(qK);
    time L<a> := FunctionField(qK : Check := true);
    time g := Genus(L);
    time Diff := DifferentDivisor(L);
    Ps,ns := Support(Diff);
    print "QUOTIENT_GENUS",g,"DIFFERENT_DEGREE",Degree(Diff),
          "RH_GENUS",1-3+Degree(Diff) div 2;
    print "QUOTIENT_DIFFERENT_SUPPORT_PROFILE",
          [ <Degree(Ps[i]),RamificationIndex(Ps[i]),ns[i]>
            : i in [1..#Ps] ];
    print "ORIGINAL_BRANCH_LEDGER",
          "R=3:e2; R=7/3:e3; R=2:e2; Q4(R):e2; Q20(R):e2";
    print "QUOTIENT_BRANCH_LEDGER",
          "v=5:e2; v=14/3:e3; v=2:e2; Q2(v):e2; Q10(v):e2";
    print "LIFT_DISCRIMINANT","9*v^2-60*v+84";

    // A tempting Mobius symmetry of the three rational branch values
    // is not an automorphism: it fails the nonrational branch factors.
    tau_factors := [ fe[1] : fe in Factorization(discq)
                     | (Degree(fe[1]) eq 2 and fe[2] eq 5)
                       or Degree(fe[1]) eq 10 ];
    for h in tau_factors do
        hn := h/LeadingCoefficient(h);
        ht := MobiusTransform(h,106,-532,21,-106);
        print "CANDIDATE_TAU_FACTOR_DEGREE",Degree(h),
              "PRESERVED",ht eq hn;
    end for;
    print "CONTACT30_C3ROOT_TRIGONAL_QUOTIENT_DONE";
    quit;
end if;

assert mode eq "quotient-search";
print "CONTACT30_C3ROOT_TRIGONAL_QUOTIENT_SEARCH","eps",eps;
hits := SearchTrigonal(q,qI,height,prime_list,"V_SEARCH");
QRR<rr> := PolynomialRing(QQ);
Qz<z> := PolynomialRing(QQ);
rows := [];
open_lifts := [];
for hit in hits do
    vv := hit[1]; rho := hit[2];
    delta := 9*vv^2-60*vv+84;
    rroots := [];
    if IsSquareQ(delta) then
        rroots := [ rt[1] : rt in Roots(3*rr^2-3*vv*rr+5*vv-7) ];
    end if;
    good := [];
    for rv in rroots do
        spec := Qz![ Evaluate(Coefficient(F,i),rv)
                     : i in [0..Degree(F)] ];
        if Evaluate(spec,rho) eq 0 and IsOpenR(rv,R,t,D) then
            Append(~good,rv);
            Append(~open_lifts,<vv,rho,rv>);
        end if;
    end for;
    Append(~rows,<vv,rho,delta,rroots,good>);
end for;
print "QUOTIENT_LIFT_ROWS",rows;
print "OPEN_LIFTS",open_lifts;
print "CONTACT30_C3ROOT_TRIGONAL_QUOTIENT_SEARCH_DONE";
quit;
