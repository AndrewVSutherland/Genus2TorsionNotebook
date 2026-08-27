//////////////////////////////////////////////////////////////////////
// Target [2,2,2,24]: halve the marked order-12 class on the direct
// A(2,2,2,12) presentation.
//
// This script supplies:
//   1. a global polynomial incidence presentation of the 2-division cover;
//   2. an exact recovery of the new order-96 record in the direct chart;
//   3. the saturated square-quartic fiber over the record;
//   4. finite-field diagnostics and a bounded rational search on the
//      transverse q-square slice through the record.
//
// Run from the repository root:
//   magma -b code/target_22224_order12_halving.m
//
// Options:
//   height:=20
//   PrimeList:=29,31,37,41
//   run_search:=false
//   run_finite:=false
//   write_equations:=false
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);

if not assigned height then
    height := 20;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;

if not assigned PrimeList then
    // The record's q-square coordinates have poles at every prime <= 23
    // except 2 and 3; these are the first useful affine slice primes.
    PrimeList := [29,31,37,41];
elif Type(PrimeList) eq MonStgElt then
    PrimeList := [StringToInteger(s) : s in Split(PrimeList,",")];
end if;

if not assigned run_search then
    run_search := true;
elif Type(run_search) eq MonStgElt then
    run_search := run_search in {"true","True","1","yes"};
end if;

if not assigned run_finite then
    run_finite := true;
elif Type(run_finite) eq MonStgElt then
    run_finite := run_finite in {"true","True","1","yes"};
end if;

if not assigned write_equations then
    write_equations := true;
elif Type(write_equations) eq MonStgElt then
    write_equations := write_equations in {"true","True","1","yes"};
end if;

if not assigned max_exact_tests then
    max_exact_tests := 500;
elif Type(max_exact_tests) eq MonStgElt then
    max_exact_tests := StringToInteger(max_exact_tests);
end if;

if not assigned log_file then
    log_file := "results/target_22224_order12_halving.log";
end if;
if not assigned equations_file then
    equations_file := "results/target_22224_global_cover_equations.txt";
end if;

SetLogFile(log_file : Overwrite := true);

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

function RationalSquareRoot(r)
    if r lt 0 then return false,Q!0; end if;
    n := Z!Numerator(r); d := Z!Denominator(r);
    okn,sn := IsSquare(n); okd,sd := IsSquare(d);
    if not okn or not okd then return false,Q!0; end if;
    return true,Q!sn/sd;
end function;

function MakeMonic(g)
    return g/LeadingCoefficient(g);
end function;

function IntegralSquareModel(f)
    den := LCM([Denominator(Coefficient(f,i)) : i in [0..Degree(f)]]);
    return Parent(f)!(den^2*f),den;
end function;

function DirectData(K,aa,bb,cc,dd,uu,tt,vv)
    PK<X> := PolynomialRing(K);
    f := X*(X+aa^2)*(X+bb^2)*(X+cc^2)*(X+dd^2);
    q := X^2+uu*X+tt^2;
    h := X^3+(1+3*uu)/2*X^2+vv*X+tt^3;
    return PK,f,q,h;
end function;

function MarkedD12Data(K,aa,bb,cc,dd,uu,tt,vv)
    PK,Xf,q,h := DirectData(K,aa,bb,cc,dd,uu,tt,vv);
    X := PK.1;
    if h^2-Xf-q^3 ne 0 then
        return false,PK!0,PK!0,PK!0,PK!0,PK!0,PK!0,PK!0;
    end if;

    g4 := (X-aa*bb)*(X-cc*dd)-X*(aa+bb)*(cc+dd);
    L4 := (X-aa*bb)*(cc+dd)+(aa+bb)*(X-cc*dd);
    v3 := h mod q;
    v4 := (-X*L4) mod g4;

    // CRT: ell=v3+q*(mc*X+nc), with ell=v4 mod g4.
    Rm := (q*X) mod g4;
    Rn := q mod g4;
    Rt := (v4-v3) mod g4;
    A0 := Coefficient(Rm,0); A1 := Coefficient(Rm,1);
    B0 := Coefficient(Rn,0); B1 := Coefficient(Rn,1);
    T0 := Coefficient(Rt,0); T1 := Coefficient(Rt,1);
    det := A0*B1-A1*B0;
    if det eq 0 then
        return false,PK!0,PK!0,PK!0,PK!0,PK!0,PK!0,PK!0;
    end if;
    mc := (T0*B1-T1*B0)/det;
    nc := (A0*T1-A1*T0)/det;
    ell := v3+q*(mc*X+nc);
    if (ell^2-Xf) mod (q*g4) ne 0 then
        return false,PK!0,PK!0,PK!0,PK!0,PK!0,PK!0,PK!0;
    end if;
    Uraw := ExactQuotient(ell^2-Xf,q*g4);
    U := MakeMonic(Uraw);
    V := (-ell) mod U;
    if (V^2-Xf) mod U ne 0 then
        return false,PK!0,PK!0,PK!0,PK!0,PK!0,PK!0,PK!0;
    end if;
    return true,Xf,q,h,g4,L4,U,V;
end function;

function HalvingEquationsQ(U,V,f)
    A<M,N> := PolynomialRing(Q,2);
    AX<X> := PolynomialRing(A);
    phi := hom<Parent(f)->AX|X>;
    UX := phi(U); VX := phi(V); fX := phi(f);
    ell := VX+UX*(M*X+N);
    if (ell^2-fX) mod UX ne 0 then
        return false,A!0,A!0,AX!0;
    end if;
    S := ExactQuotient(ell^2-fX,UX);
    s4:=Coefficient(S,4); s3:=Coefficient(S,3);
    s2:=Coefficient(S,2); s1:=Coefficient(S,1); s0:=Coefficient(S,0);
    E1:=8*s4^2*s1-s3*(4*s4*s2-s3^2);
    E0:=64*s4^3*s0-(4*s4*s2-s3^2)^2;
    return true,E1,E0,S;
end function;

function HalvingEquationsFinite(U,V,f)
    K := BaseRing(Parent(f));
    A<M,N> := PolynomialRing(K,2);
    AX<X> := PolynomialRing(A);
    phi := hom<Parent(f)->AX|X>;
    UX:=phi(U); VX:=phi(V); fX:=phi(f);
    ell:=VX+UX*(M*X+N);
    if (ell^2-fX) mod UX ne 0 then
        return false,A!0,A!0,AX!0;
    end if;
    S:=ExactQuotient(ell^2-fX,UX);
    s4:=Coefficient(S,4); s3:=Coefficient(S,3);
    s2:=Coefficient(S,2); s1:=Coefficient(S,1); s0:=Coefficient(S,0);
    E1:=8*s4^2*s1-s3*(4*s4*s2-s3^2);
    E0:=64*s4^3*s0-(4*s4*s2-s3^2)^2;
    return true,E1,E0,S;
end function;

function CountHalvingFiberFinite(E1,E0,S)
    K:=BaseRing(Parent(E1));
    s4:=Coefficient(S,4);
    n:=0; boundary:=0;
    for mm in K do for nn in K do
        vals:=[mm,nn];
        if Evaluate(E1,vals) eq 0 and Evaluate(E0,vals) eq 0 then
            if Evaluate(s4,vals) eq 0 then boundary+:=1;
            else n+:=1;
            end if;
        end if;
    end for; end for;
    return n,boundary;
end function;

function LocalHalvingCountQ(E1,E0,S,p)
    try
        Fp:=GF(p); Ap<m,n>:=PolynomialRing(Fp,2);
        e1:=Ap!E1; e0:=Ap!E0;
        // Rebuild only the leading coefficient of S in the same ring.
        s4:=Ap!Coefficient(S,4);
        count:=0; boundary:=0;
        for mm in Fp do for nn in Fp do
            vals:=[mm,nn];
            if Evaluate(e1,vals) eq 0 and Evaluate(e0,vals) eq 0 then
                if Evaluate(s4,vals) eq 0 then boundary+:=1;
                else count+:=1;
                end if;
            end if;
        end for; end for;
        return true,count,boundary;
    catch e
        return false,0,0;
    end try;
end function;

function RootPowerWitness(f)
    C:=HyperellipticCurve(f);
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97] do
        try
            fp:=ChangeRing(f,GF(p));
            if Degree(fp) notin {5,6} or Discriminant(fp) eq 0 then continue; end if;
            chi:=LPolynomial(ChangeRing(C,GF(p)));
            if not IsIrreducible(chi) then continue; end if;
            R<T>:=PolynomialRing(Q); chQ:=R!chi;
            K<pi>:=NumberField(chQ);
            degs:=[Degree(MinimalPolynomial(pi^n)):n in [2..12]];
            if &and[d eq 4:d in degs] then return true,p,chQ,degs; end if;
        catch e
            continue;
        end try;
    end for;
    R<T>:=PolynomialRing(Q);
    return false,0,R!0,[];
end function;

//////////////////////////////////////////////////////////////////////
// Global polynomial incidence presentation.
//////////////////////////////////////////////////////////////////////

procedure WriteGlobalCover(filename)
    R<a,b,c,d,u,t,v,mc,nc,A0,B0,M,N,s4,s3,s2,s1,s0> :=
        PolynomialRing(Q,18,"grevlex");
    RX<X>:=PolynomialRing(R);
    f:=X*(X+a^2)*(X+b^2)*(X+c^2)*(X+d^2);
    q:=X^2+u*X+t^2;
    h:=X^3+(1+3*u)/2*X^2+v*X+t^3;
    g4:=(X-a*b)*(X-c*d)-X*(a+b)*(c+d);
    L4:=(X-a*b)*(c+d)+(a+b)*(X-c*d);
    v3:=h mod q;
    v4:=(-X*L4) mod g4;
    ell:=v3+q*(mc*X+nc);
    crt:=(ell-v4) mod g4;
    U:=X^2+A0*X+B0;
    d12identity:=ell^2-f-mc^2*q*g4*U;
    VD:=-(ell mod U);
    ellH:=VD+U*(M*X+N);
    Sq:=s4*X^4+s3*X^3+s2*X^2+s1*X+s0;
    halfidentity:=ellH^2-f-U*Sq;
    E1:=8*s4^2*s1-s3*(4*s4*s2-s3^2);
    E0:=64*s4^3*s0-(4*s4*s2-s3^2)^2;
    contact:=[Coefficient(h^2-q^3-f,i):i in [1..4]];
    crteqs:=[Coefficient(crt,i):i in [0..1]];
    // d12identity has degree at most 5.  Its x^6 coefficient is identically
    // zero, so do not print a spurious equation of total degree -1.
    d12eqs:=[Coefficient(d12identity,i):i in [0..5]];
    halfeqs:=[Coefficient(halfidentity,i):i in [0..6]];

    print "GLOBAL_COVER_VARIABLES",["a","b","c","d","u","t","v","mc","nc","A0","B0","M","N","s4","s3","s2","s1","s0"];
    print "GLOBAL_COVER_EQUATION_COUNTS","contact",#contact,"crt",#crteqs,
          "D12_identity",#d12eqs,"halving_identity",#halfeqs,"square_quartic",2;
    print "GLOBAL_COVER_TOTAL_DEGREES",
          "contact",[TotalDegree(e):e in contact],
          "crt",[TotalDegree(e):e in crteqs],
          "D12",[TotalDegree(e):e in d12eqs],
          "half",[TotalDegree(e):e in halfeqs],
          "E",[TotalDegree(E1),TotalDegree(E0)];
    print "GLOBAL_COVER_OPEN","mc*det(q,g4)*s4*disc(f) != 0 (det is the CRT determinant)";

    if write_equations then
        out:=Open(filename,"w");
        fprintf out,"# Global A(2,2,2,24) incidence cover\n";
        fprintf out,"# variables a,b,c,d,u,t,v,mc,nc,A0,B0,M,N,s4,s3,s2,s1,s0\n";
        fprintf out,"# Base: y^2=f, q=x^2+ux+t^2, h=x^3+(1+3u)/2*x^2+v*x+t^3\n";
        fprintf out,"# D3=[q,h mod q], D4=[g4,-xL4 mod g4]\n";
        fprintf out,"# ell=v3+q*(mc*x+nc), D12=[x^2+A0*x+B0, -ell mod U]\n";
        fprintf out,"# Half chart: ellH=VD+U*(M*x+N), (ellH^2-f)/U=S square quartic\n";
        for i in [1..#contact] do fprintf out,"contact_%o = %o\n",i,contact[i]; end for;
        for i in [1..#crteqs] do fprintf out,"crt_%o = %o\n",i-1,crteqs[i]; end for;
        for i in [1..#d12eqs] do fprintf out,"d12_%o = %o\n",i-1,d12eqs[i]; end for;
        fprintf out,"VD = %o\n",VD;
        for i in [1..#halfeqs] do fprintf out,"half_%o = %o\n",i-1,halfeqs[i]; end for;
        fprintf out,"E1 = %o\nE0 = %o\n",E1,E0;
        fprintf out,"# open: discriminant(f)*mc*s4*CRTdet != 0\n";
        delete out;
        print "GLOBAL_COVER_EQUATIONS_WRITTEN",filename;
    end if;
end procedure;

//////////////////////////////////////////////////////////////////////
// Recover the record in the direct chart and analyze its halving fiber.
//////////////////////////////////////////////////////////////////////

procedure AnalyzeRecord()
    fs:=3027600*x^6+2950382280*x^5+602288814361*x^4
       -63417934304484*x^3-2122595910966478*x^2
       +128056619498204124*x+3322970988364151397;
    roots:=[Q!(-75593/145),-519,-39,Q!(-299/12),59,Q!(123109/1740)];
    print "RECORD_ROOTS",roots;

    // Mobius orientation: -519 -> 0, -39 -> infinity; normalize the
    // root -75593/145 to -1.  The other finite roots are negative squares.
    rzero:=Q!(-519); rinf:=Q!(-39); rbase:=Q!(-75593/145);
    zbase:=(rbase-rzero)/(rbase-rinf);
    zs:=[-((r-rzero)/(r-rinf))/zbase:r in roots|r notin {rzero,rinf}];
    a_mid:=Q!1; b_mid:=Q!(14399/169); c_mid:=Q!(3179/91); d_mid:=Q!(189431/5681);
    print "RECORD_MOBIUS","zero",rzero,"infinity",rinf,"zbase",zbase,
          "finite_images",zs,"negative_square_parameters",[a_mid,b_mid,c_mid,d_mid];
    assert Sort(zs) eq Sort([-a_mid^2,-b_mid^2,-c_mid^2,-d_mid^2]);
    fmid:=x*(x+a_mid^2)*(x+b_mid^2)*(x+c_mid^2)*(x+d_mid^2);
    iso,mp:=IsIsomorphic(HyperellipticCurve(fs),HyperellipticCurve(fmid));
    print "RECORD_INTERMEDIATE_ISOMORPHIC",iso,mp;
    assert iso;

    // The order-3 Mumford polynomial and contact identity on this model.
    qmid:=x^2+69938/169*x+1222830961/28561;
    vmid:=(472377235240565294640*x-28973884423844905905360)/7633184722873801;
    contact_scale:=Q!(1183/209760);
    contact_n:=Q!(22289207089/248146080);
    Hmid:=vmid+qmid*(contact_scale*x+contact_n);
    print "RECORD_INTERMEDIATE_CONTACT","q",qmid,"vmod",vmid,
          "a3",contact_scale,"n",contact_n,"H",Hmid,
          "identity",Hmid^2-fmid eq contact_scale^2*qmid^3;
    assert Hmid^2-fmid eq contact_scale^2*qmid^3;

    // Scale X=a3^2*x and Y=a3^5*y so that the contact cubic has a3=1.
    aa:=Q!(1183/209760);
    bb:=Q!(100793/209760);
    cc:=Q!(41327/209760);
    dd:=Q!(17238221/91665120);
    uu:=Q!(289578289/21999628800);
    tt:=Q!(289578289/43999257600);
    vv:=Q!(11494656094062061921/645311556450385920000);
    print "RECORD_DIRECT_COORDINATES","a",aa,"b",bb,"c",cc,"d",dd,
          "u",uu,"t",tt,"v",vv;
    ok,f,q,h,g4,L4,U,V:=MarkedD12Data(Q,aa,bb,cc,dd,uu,tt,vv);
    assert ok;
    print "RECORD_DIRECT_IDENTITIES","contact",h^2-f eq q^3,
          "q_square",q eq (x+tt)^2,
          "D12_u",U,"D12_v",V;
    iso2,mp2:=IsIsomorphic(HyperellipticCurve(fs),HyperellipticCurve(f));
    print "RECORD_DIRECT_ISOMORPHIC",iso2,mp2;
    assert iso2;

    J:=Jacobian(HyperellipticCurve(f)); O:=J!0;
    D3:=J![q,h mod q];
    D4:=J![g4,(-x*L4) mod g4];
    D12:=J![U,V];
    print "RECORD_MARKED_ORDERS",Order(D3),Order(D4),Order(D12),
          "cantor_match",D12 eq D3+D4;
    assert D12 eq D3+D4 and Order(D12) eq 12;

    // q-square slice coordinates through this point.
    ss:=Q!(17017/209760);
    Arec:=Q!(13/187); Brec:=Q!(77/13); Crec:=Q!(17/7);
    rhorec:=Q!(104880/17017); sigmarec:=Q!(243120/17017);
    print "RECORD_QSQUARE_SLICE","s",ss,"A",Arec,"B",Brec,"C",Crec,
          "rho",rhorec,"sigma",sigmarec,
          "ABC",Arec*Brec*Crec,
          "Rcheck",Arec^2+Brec^2+Crec^2-3 eq rhorec^2,
          "Scheck",1/Arec^2+1/Brec^2+1/Crec^2-3 eq sigmarec^2;
    assert tt eq ss^2 and uu eq 2*tt;

    okH,E1,E0,S:=HalvingEquationsQ(U,V,f);
    assert okH;
    A2:=Parent(E1); M:=A2.1; N:=A2.2;
    s4:=Coefficient(S,4);
    print "RECORD_HALVING_EQUATIONS","E1_degrees",
          [TotalDegree(E1),Degree(E1,M),Degree(E1,N)],
          "E0_degrees",[TotalDegree(E0),Degree(E0,M),Degree(E0,N)],
          "s4",s4,"gcd",GCD(E1,E0);

    for p in PrimeList do
        good,n,boundary:=LocalHalvingCountQ(E1,E0,S,p);
        print "RECORD_HALVING_MOD",p,"good_prime_for_equations",good,
              "saturated_points",n,"boundary_points",boundary;
    end for;

    RN:=Resultant(E1,E0,N);
    fac:=Factorization(RN);
    affine:=[<Degree(fe[1],M),fe[2]>:fe in fac|GCD(fe[1],M) eq 1];
    boundary:=[<Degree(fe[1],M),fe[2]>:fe in fac|GCD(fe[1],M) ne 1];
    print "RECORD_HALVING_RESULTANT","degree",Degree(RN,M),
          "boundary_factors",boundary,"saturated_factor_degrees",affine,
          "saturated_degree",&+[z[1]*z[2]:z in affine];

    fI,scale:=IntegralSquareModel(f);
    JI:=Jacobian(HyperellipticCurve(fI));
    DI:=JI![U,scale*V];
    divisible,half:=IsDivisibleBy(DI,2);
    TI,tmp:=TorsionSubgroup(JI);
    print "RECORD_EXACT_INTEGRAL","scale",scale,
          "torsion",Invariants(TI),"D12_order",Order(DI),
          "D12_divisible_by_2",divisible;
    assert Invariants(TI) eq [2,2,2,12] and not divisible;
end procedure;

//////////////////////////////////////////////////////////////////////
// The q-square slice and bounded transverse search.
//
// Freely vary A,B; put C=1/(AB).  The constrained double cover is
//   A^2+B^2+C^2-3=rho^2,
//   A^-2+B^-2+C^-2-3=sigma^2.
// Then s=1/(2rho), t=s^2, u=2t,
//   (a,b,c) = s*(A,B,C), d=2s^2*sigma,
//   v=3t^2+d^2/2.
//////////////////////////////////////////////////////////////////////

function PositiveHeightRationals(H)
    vals:=[]; seen:={};
    for den in [1..H] do for num in [1..H] do
        if GCD(num,den) ne 1 then continue; end if;
        z:=Q!num/den; key:=Sprint(z);
        if key notin seen then Include(~seen,key); Append(~vals,z); end if;
    end for; end for;
    return Sort(vals);
end function;

function ExactDivisibilityTest(f,U,V)
    fI,scale:=IntegralSquareModel(f);
    try
        JI:=Jacobian(HyperellipticCurve(fI));
        D:=JI![U,scale*V];
        if Order(D) ne 12 then return false,false,[],fI; end if;
        divisible,half:=IsDivisibleBy(D,2);
        if not divisible then return true,false,[],fI; end if;
        T,mp:=TorsionSubgroup(JI);
        return true,true,Invariants(T),fI;
    catch e
        return false,false,[],fI;
    end try;
end function;

procedure RationalQSquareSearch(H)
    vals:=PositiveHeightRationals(H);
    pairs:=0; double_square:=0; collision_bases:=0;
    smooth_bases:=0; marked_classes:=0;
    local_survivors:=0; exact_tests:=0; hits:=0;
    seen_bases:={};
    print "QSQUARE_SEARCH_START","height",H,"positive_values",#vals;
    for A0 in vals do for B0 in vals do
        C0:=1/(A0*B0);
        if not (A0 le B0 and B0 le C0) then continue; end if;
        pairs+:=1;
        R0:=A0^2+B0^2+C0^2-3;
        S0:=1/A0^2+1/B0^2+1/C0^2-3;
        okr,rho:=RationalSquareRoot(R0);
        oks,sigma:=RationalSquareRoot(S0);
        if not okr or not oks or rho eq 0 or sigma eq 0 then continue; end if;
        double_square+:=1;
        ss:=1/(2*rho); tt:=ss^2; uu:=2*tt;
        mags:=[ss*A0,ss*B0,ss*C0];
        ddmag:=2*ss^2*sigma;
        // This catches the rational component {A,B,C}={r,1,1/r}:
        // it satisfies both square equations identically but gives d^2=a_i^2
        // and hence a singular curve.
        if #{z^2:z in mags cat [ddmag]} lt 4 then
            collision_bases+:=1;
            continue;
        end if;
        basekey:=Join([Sprint(z^2):z in Sort(mags cat [ddmag])],",");
        if basekey in seen_bases then continue; end if;
        Include(~seen_bases,basekey);

        // Three choices for which cubic root is paired with the linear root
        // in D4, all sign patterns on the cubic roots, and both signs of d.
        choices:=[<mags[1],mags[2],mags[3]>,
                  <mags[1],mags[3],mags[2]>,
                  <mags[2],mags[3],mags[1]>];
        base_smooth:=false;
        for ch in choices do
          for sa in [-1,1] do for sb in [-1,1] do for sc in [-1,1] do
            for sd in [-1,1] do
                aa:=sa*ch[1]; bb:=sb*ch[2]; cc:=sc*ch[3]; dd:=sd*ddmag;
                vv:=3*tt^2+dd^2/2;
                ok,f,q,h,g4,L4,U,V:=MarkedD12Data(Q,aa,bb,cc,dd,uu,tt,vv);
                if not ok or Discriminant(f) eq 0 then continue; end if;
                if not base_smooth then smooth_bases+:=1; base_smooth:=true; end if;
                marked_classes+:=1;
                okH,E1,E0,S:=HalvingEquationsQ(U,V,f);
                if not okH then continue; end if;
                local_ok:=true;
                for p in PrimeList do
                    good,n,boundary:=LocalHalvingCountQ(E1,E0,S,p);
                    if good and n eq 0 then local_ok:=false; break; end if;
                end for;
                if not local_ok then continue; end if;
                local_survivors+:=1;
                if exact_tests ge max_exact_tests then continue; end if;
                exact_tests+:=1;
                tested,divisible,inv,fI:=ExactDivisibilityTest(f,U,V);
                if tested and divisible then
                    hits+:=1;
                    simple,pcert,chi,degs:=RootPowerWitness(fI);
                    print "TARGET_22224_HIT","A",A0,"B",B0,"C",C0,
                          "rho",rho,"sigma",sigma,
                          "a",aa,"b",bb,"c",cc,"d",dd,
                          "u",uu,"t",tt,"v",vv,
                          "torsion",inv,"simple",simple,"pcert",pcert,
                          "chi",chi,"degrees",degs,"curve",fI;
                end if;
            end for;
          end for; end for; end for;
        end for;
    end for; end for;
    print "QSQUARE_SEARCH_DONE","height",H,"pairs",pairs,
          "double_square_points",double_square,"unique_bases",#seen_bases,
          "branch_collision_bases",collision_bases,
          "smooth_bases",smooth_bases,"marked_classes",marked_classes,
          "local_survivors",local_survivors,"exact_tests",exact_tests,
          "hits",hits;
end procedure;

procedure FiniteQSquareStats(p)
    K:=GF(p);
    bases:=0; collision_bases:=0; smooth:=0; classes:=0;
    contact_fail:=0; chart_fail:=0; target_classes:=0; target_points:=0;
    for A0 in K do
      if A0 eq 0 then continue; end if;
      for B0 in K do
        if B0 eq 0 then continue; end if;
        C0:=1/(A0*B0);
        R0:=A0^2+B0^2+C0^2-3;
        S0:=1/A0^2+1/B0^2+1/C0^2-3;
        if R0 eq 0 or S0 eq 0 or not IsSquare(R0) or not IsSquare(S0) then continue; end if;
        rho:=SquareRoot(R0); sigma:=SquareRoot(S0);
        ss:=1/(2*rho); tt:=ss^2; uu:=2*tt;
        mags:=[ss*A0,ss*B0,ss*C0]; ddmag:=2*ss^2*sigma;
        bases+:=1;
        if #{z^2:z in mags cat [ddmag]} lt 4 then
            collision_bases+:=1;
            continue;
        end if;
        base_smooth:=false;
        choices:=[<mags[1],mags[2],mags[3]>,
                  <mags[1],mags[3],mags[2]>,
                  <mags[2],mags[3],mags[1]>];
        // This is a labelled six-class sample: the three choices of the
        // distinguished cubic-root pairing and the two signs of d.  The
        // rational search below is the exhaustive sign-pattern search.
        for ch in choices do for sd in [-1,1] do
            aa:=ch[1];bb:=ch[2];cc:=ch[3];dd:=K!sd*ddmag;
            vv:=3*tt^2+dd^2/2;
            PK,f0,q0,h0:=DirectData(K,aa,bb,cc,dd,uu,tt,vv);
            if h0^2-f0-q0^3 ne 0 then contact_fail+:=1; continue; end if;
            if Discriminant(f0) eq 0 then continue; end if;
            ok,f,q,h,g4,L4,U,V:=MarkedD12Data(K,aa,bb,cc,dd,uu,tt,vv);
            if not ok then chart_fail+:=1; continue; end if;
            if not base_smooth then smooth+:=1;base_smooth:=true;end if;
            classes+:=1;
            okH,E1,E0,S:=HalvingEquationsFinite(U,V,f);
            if not okH then continue; end if;
            n,boundary:=CountHalvingFiberFinite(E1,E0,S);
            if n gt 0 then target_classes+:=1;target_points+:=n;end if;
        end for; end for;
      end for;
    end for;
    print "QSQUARE_FINITE",p,"double_square_oriented_bases",bases,
          "branch_collision_bases",collision_bases,
          "smooth_bases",smooth,"marked_classes_sampled",classes,
          "contact_failures",contact_fail,"CRT_chart_failures",chart_fail,
          "target_classes",target_classes,"target_halves",target_points;
end procedure;

print "TARGET_22224_ORDER12_HALVING_START";
print "CONFIG","height",height,"PrimeList",PrimeList,
      "run_search",run_search,"run_finite",run_finite,
      "write_equations",write_equations;
WriteGlobalCover(equations_file);
AnalyzeRecord();
if run_finite then
    for p in PrimeList do FiniteQSquareStats(p); end for;
end if;
if run_search then RationalQSquareSearch(height); end if;
print "TARGET_22224_ORDER12_HALVING_DONE";
UnsetLogFile();
quit;
