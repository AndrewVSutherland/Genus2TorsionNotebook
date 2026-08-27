//////////////////////////////////////////////////////////////////////
// Exact genus computation for the orthogonal relative-3 support cover
// on the rational P8 component, and for its signed L^2=M double cover.
//
// This is deliberately independent of the exact reconstruction driver:
// it reconstructs only f12 and M(v), then performs the P8 base change.
//
// Usage:
//   magma -b code/contact6_m612_p8_relative3_exact_genus.m
//   magma -b do_signed:=false code/contact6_m612_p8_relative3_exact_genus.m
//////////////////////////////////////////////////////////////////////

SetColumns(0); SetSeed(1);
if not assigned do_signed then do_signed:=true;
elif Type(do_signed) eq MonStgElt then
    do_signed:=do_signed in {"true","True","1","yes"};
end if;
if not assigned relative_only then relative_only:=false;
elif Type(relative_only) eq MonStgElt then
    relative_only:=relative_only in {"true","True","1","yes"};
end if;
if not assigned discriminant_only then discriminant_only:=false;
elif Type(discriminant_only) eq MonStgElt then
    discriminant_only:=discriminant_only in {"true","True","1","yes"};
end if;
if not assigned ramification_only then ramification_only:=false;
elif Type(ramification_only) eq MonStgElt then
    ramification_only:=ramification_only in {"true","True","1","yes"};
end if;
if not assigned m_divisor_only then m_divisor_only:=false;
elif Type(m_divisor_only) eq MonStgElt then
    m_divisor_only:=m_divisor_only in {"true","True","1","yes"};
end if;
if not assigned automorphism_only then automorphism_only:=false;
elif Type(automorphism_only) eq MonStgElt then
    automorphism_only:=automorphism_only in {"true","True","1","yes"};
end if;
if not assigned p8_rh_only then p8_rh_only:=false;
elif Type(p8_rh_only) eq MonStgElt then
    p8_rh_only:=p8_rh_only in {"true","True","1","yes"};
end if;
if p8_rh_only then relative_only:=true; end if;

Q:=Rationals(); Qz<z>:=PolynomialRing(Q);

function FixedRelative12AndM(ee)
    R0<M0,U0,v0>:=PolynomialRing(Q,3,"grevlex");
    N0:=(3*U0+6*M0)/2;
    Rc0:=(3*U0^2+3*v0^2+(2/ee-15)*M0-N0^2)/2;
    Fs0:=[2*v0^3+2*N0*Rc0-U0^3-6*U0*v0^2-22*M0,
          Rc0^2+2*N0*v0^3-3*U0^2*v0^2-3*v0^4
              -(1/ee^2-15)*M0,
          2*Rc0*v0^3-3*U0*v0^4-(2/ee+6)*M0];
    I0:=Saturation(ideal<R0|Fs0>,ideal<R0|M0>);
    assert Dimension(I0) eq 0 and Dimension(quo<R0|I0>) eq 40;
    G0:=GroebnerBasis(ChangeOrder(I0,"lex"));
    univ:=[g:g in G0|Degree(g,1) eq 0 and Degree(g,2) eq 0 and
                       Degree(g,3) gt 0];
    assert #univ eq 1;
    rv:=Qz!UnivariatePolynomial(univ[1]);
    rv/:=LeadingCoefficient(rv);
    f12s:=[fe[1]:fe in Factorization(rv)|
                    Degree(fe[1]) eq 12 and fe[2] eq 1];
    assert #f12s eq 1;
    f12:=f12s[1]/LeadingCoefficient(f12s[1]);
    gm:=[g:g in G0|Degree(g,1) eq 1 and Degree(g,2) eq 0][1];
    cm:=Q!Coefficient(gm,1,1);
    mrec:=(-Qz!UnivariatePolynomial(
                Evaluate(gm,[R0!0,R0!0,v0]))/cm) mod f12;
    return Qz!f12,Qz!mrec;
end function;

sample_e:=[Q!j:j in [1..20]];
samples:=[]; msamples:=[];
for ee in sample_e do
    f,mf:=FixedRelative12AndM(ee);
    Append(~samples,f); Append(~msamples,mf);
end for;

K<e>:=FunctionField(Q); Kz<zz>:=PolynomialRing(K);
function ReconstructPolynomial(poly_samples,degree,bds,sample_e,Q,K,e,Kz)
    coeffs:=[];
    for i in [0..degree] do
        nd:=bds[i+1][1]; dd:=bds[i+1][2];
        rows:=[];
        for j in [1..#sample_e] do
            x:=sample_e[j]; y:=Coefficient(poly_samples[j],i);
            Append(~rows,[x^h:h in [0..nd]] cat
                         [-y*x^h:h in [0..dd]]);
        end for;
        A:=Matrix(Q,#rows,nd+dd+2,&cat rows);
        W:=Nullspace(Transpose(A)); assert Dimension(W) eq 1;
        w:=Basis(W)[1];
        num:=&+[w[h+1]*e^h:h in [0..nd]];
        den:=&+[w[nd+2+h]*e^h:h in [0..dd]];
        assert den ne 0;
        c:=num/den;
        assert &and[Evaluate(Numerator(c),sample_e[j])/
                    Evaluate(Denominator(c),sample_e[j]) eq
                    Coefficient(poly_samples[j],i)
                    :j in [1..#sample_e]];
        Append(~coeffs,c);
    end for;
    return Kz!coeffs;
end function;

fbounds:=[<4,4>,<3,4>,<4,4>,<3,3>,<3,3>,<3,3>,<2,3>,
          <3,3>,<3,3>,<2,2>,<1,1>,<0,1>,<0,0>];
mbounds:=[<6,3>,<6,3>,<6,3>,<6,3>,<7,4>,<6,4>,
          <7,5>,<8,5>,<8,5>,<8,5>,<8,5>,<8,5>];
f12:=ReconstructPolynomial(samples,12,fbounds,sample_e,Q,K,e,Kz);
mrec:=ReconstructPolynomial(msamples,11,mbounds,sample_e,Q,K,e,Kz)
      mod f12;
assert Degree(f12) eq 12 and IsMonic(f12) and IsIrreducible(f12);

if relative_only then
    print "CONTACT6_M612_RELATIVE3_EXACT_GENUS";
    if discriminant_only then
        dr:=Discriminant(f12);
        print "RELATIVE_PRIMITIVE_DISCRIMINANT_DEGREES",
              <Degree(Numerator(dr)),Degree(Denominator(dr))>;
        fdrn:=Factorization(Numerator(dr));
        fdrd:=Factorization(Denominator(dr));
        print "RELATIVE_PRIMITIVE_DISCRIMINANT_NUMERATOR_FACTORS",
              [<Degree(fe[1]),fe[2]>:fe in fdrn],fdrn;
        print "RELATIVE_PRIMITIVE_DISCRIMINANT_DENOMINATOR_FACTORS",
              [<Degree(fe[1]),fe[2]>:fe in fdrd],fdrd;
        print "CONTACT6_M612_RELATIVE3_EXACT_DISCRIMINANT_DONE";
        quit;
    end if;
    Frel<wrel>:=ext<K|f12>;
    if p8_rh_only then
        print "CONTACT6_M612_P8_RELATIVE3_EXACT_RH";
        KuRH<uRH>:=FunctionField(Q);
        tRH:=4*(uRH^2+uRH-6)/(uRH^2+6);
        eRH:=-(Q!25/3)*tRH^2/
              (tRH^4-25*tRH^2+Q!1250/3);

        Drel:=DifferentDivisor(Frel); Pd,nd:=Support(Drel);
        print "RELATIVE_DIFFERENT_DEGREE",Degree(Drel);
        totaldiff:=0;
        for i in [1..#Pd] do
            q:=Minimum(Pd[i]); r:=RamificationIndex(Pd[i]);
            fpull:=Factorization(Numerator(Evaluate(q,eRH)));
            c:=Degree(Pd[i]) div Degree(q) *
               &+[Degree(fe[1])*(r-GCD(r,fe[2])):fe in fpull];
            totaldiff+:=c;
            print "P8_DIFFERENT_PULLBACK",q,
                  "SOURCE_PLACE",<Degree(Pd[i]),r,nd[i]>,
                  "BASE_FACTORS",[<Degree(fe[1]),fe[2]>:fe in fpull],
                  "CONTRIBUTION",c;
        end for;
        g12RH:=1-12+totaldiff div 2;
        print "P8_TOTAL_DIFFERENT_DEGREE",totaldiff;
        print "P8_SUPPORT_GENUS_BY_RH",g12RH;
        assert Degree(Drel) eq 32 and totaldiff eq 168 and g12RH eq 73;

        mrelRH:=Evaluate(mrec,wrel);
        DmRH:=Divisor(mrelRH); PmRH,nmRH:=Support(DmRH);
        signram:=0;
        for i in [1..#PmRH] do
            q:=Minimum(PmRH[i]); r:=RamificationIndex(PmRH[i]);
            fpull:=Factorization(Numerator(Evaluate(q,eRH)));
            ci:=0;
            newvals:=[];
            for fe in fpull do
                gg:=GCD(r,fe[2]); nv:=nmRH[i]*fe[2] div gg;
                Append(~newvals,<Degree(fe[1]),fe[2],gg,nv>);
                if IsOdd(nv) then
                    ci+:=Degree(PmRH[i]) div Degree(q)*Degree(fe[1])*gg;
                end if;
            end for;
            signram+:=ci;
            print "P8_M_DIVISOR_PULLBACK",q,
                  "SOURCE_PLACE",<Degree(PmRH[i]),r,nmRH[i]>,
                  "PULLBACK_DATA",newvals,"ODD_DEGREE_CONTRIBUTION",ci;
        end for;
        g24RH:=2*g12RH-1+signram div 2;
        print "P8_SIGN_QUADRATIC_RAMIFICATION_DEGREE",signram;
        print "P8_SIGNED_GENUS_BY_RH",g24RH;
        assert signram eq 0 and g24RH eq 145;

        // Exact connectedness check for the signed cover.  This is the
        // norm polynomial of L^2-M(v), pulled back from Q(e) to Q(u).
        KTc<Tc>:=PolynomialRing(K); KTcV<Vc>:=PolynomialRing(KTc);
        fcV:=&+[(KTc!Coefficient(f12,i))*Vc^i:i in [0..12]];
        mcV:=&+[(KTc!Coefficient(mrec,i))*Vc^i
                :i in [0..Degree(mrec)]];
        p24c:=Resultant(fcV,Tc^2-mcV);
        p24c/:=LeadingCoefficient(p24c);
        phiRH:=hom<K -> KuRH | eRH>;
        KuRHT<TT>:=PolynomialRing(KuRH);
        p24RH:=KuRHT![phiRH(Coefficient(p24c,i)):i in [0..24]];
        tirr:=Cputime(); irr24RH:=IsIrreducible(p24RH);
        print "P8_SIGNED_COVER_CONNECTED",irr24RH,
              "IRREDUCIBILITY_SECONDS",Cputime(tirr);
        assert irr24RH;
        print "CONTACT6_M612_P8_RELATIVE3_EXACT_RH_DONE";
        quit;
    end if;
    if automorphism_only then
        print "RELATIVE_SIGNED_AUTOMORPHISM_START";
        KT<T>:=PolynomialRing(K); KTv<VV>:=PolynomialRing(KT);
        f12V:=&+[(KT!Coefficient(f12,i))*VV^i:i in [0..12]];
        mrecV:=&+[(KT!Coefficient(mrec,i))*VV^i
                  :i in [0..Degree(mrec)]];
        tp:=Cputime(); p24:=Resultant(f12V,T^2-mrecV);
        p24/:=LeadingCoefficient(p24);
        print "RELATIVE_SIGNED_PRIMITIVE_DEGREE",Degree(p24),
              "RESULTANT_SECONDS",Cputime(tp);
        assert Degree(p24) eq 24 and IsIrreducible(p24);
        F24<ell>:=ext<K|p24>; F24Y<Y>:=PolynomialRing(F24);
        p24F:=&+[(F24!Coefficient(p24,i))*Y^i:i in [0..24]];
        tf:=Cputime(); fac24:=Factorization(p24F);
        print "RELATIVE_SIGNED_INTERNAL_FACTOR_DEGREES",
              [<Degree(fe[1]),fe[2]>:fe in fac24],
              "SECONDS",Cputime(tf);
        linroots:=[-Coefficient(fe[1],0)/Coefficient(fe[1],1)
                  :fe in fac24|Degree(fe[1]) eq 1];
        print "RELATIVE_SIGNED_DECK_ROOT_COUNT",#linroots;
        for j in [1..#linroots] do
            hj:=hom<F24 -> F24 | linroots[j]>;
            r2:=hj(linroots[j]); r3:=hj(r2);
            ord:=linroots[j] eq ell select 1 else
                 (r2 eq ell select 2 else
                 (r3 eq ell select 3 else 0));
            print "RELATIVE_SIGNED_DECK_MAP",j,"ORDER_UP_TO_3",ord;
            if ord eq 3 then
                s3:=ell+linroots[j]+r2;
                mp8:=MinimalPolynomial(s3);
                print "RELATIVE_C3_INVARIANT_MINPOLY_DEGREE",Degree(mp8);
                FY<YY>:=PolynomialRing(F24);
                hneg:=hom<F24 -> F24 | -ell>;
                s6:=s3+hneg(s3);
                mp4:=MinimalPolynomial(s6);
                print "RELATIVE_S3_INVARIANT_MINPOLY_DEGREE",Degree(mp4);
                E8<a8>:=ext<K|mp8>;
                E4<a4>:=ext<K|mp4>;
                print "RELATIVE_C3_QUOTIENT_GENUS",Genus(E8);
                print "RELATIVE_S3_QUOTIENT_GENUS",Genus(E4);
                break;
            end if;
        end for;
        print "CONTACT6_M612_RELATIVE3_SIGNED_AUTOMORPHISM_DONE";
        quit;
    end if;
    if ramification_only then
        branchpolys:=[fe[1]:fe in Factorization(Numerator(
                                              Discriminant(f12)))];
        Append(~branchpolys,Parent(branchpolys[1]).1);
        totaldiff:=0;
        for q in branchpolys do
            P:=Zeros(K!q)[1];
            tb:=Cputime(); dec:=Decomposition(Frel,P);
            dat:=[<Degree(QQ),RamificationIndex(QQ),
                   InertiaDegree(QQ)>:QQ in dec];
            contr:=&+[Degree(QQ)*(RamificationIndex(QQ)-1):QQ in dec];
            totaldiff+:=contr;
            print "RELATIVE_FINITE_BRANCH",q,"BASE_DEGREE",Degree(P),
                  "DECOMPOSITION",dat,"DIFFERENT_CONTRIBUTION",contr,
                  "SECONDS",Cputime(tb);
        end for;
        Pinf:=Poles(K!e)[1];
        tbi:=Cputime(); decinf:=Decomposition(Frel,Pinf);
        datinf:=[<Degree(QQ),RamificationIndex(QQ),
                  InertiaDegree(QQ)>:QQ in decinf];
        contrinf:=&+[Degree(QQ)*(RamificationIndex(QQ)-1)
                     :QQ in decinf];
        totaldiff+:=contrinf;
        print "RELATIVE_INFINITE_BRANCH_DECOMPOSITION",datinf,
              "DIFFERENT_CONTRIBUTION",contrinf,
              "SECONDS",Cputime(tbi);
        print "RELATIVE_TOTAL_TAME_DIFFERENT",totaldiff;
        print "RELATIVE_RH_GENUS",1-12+totaldiff div 2;
        print "CONTACT6_M612_RELATIVE3_EXACT_RAMIFICATION_DONE";
        quit;
    end if;
    if m_divisor_only then
        print "RELATIVE_M_DIVISOR_START";
        tm:=Cputime(); mrel:=Evaluate(mrec,wrel);
        Dm:=Divisor(mrel); Pm,nm:=Support(Dm);
        print "RELATIVE_M_DIVISOR_DEGREE",Degree(Dm),
              "SUPPORT_SIZE",#Pm,"SECONDS",Cputime(tm);
        print "RELATIVE_M_DIVISOR_DATA",
              [<nm[i],Degree(Pm[i]),RamificationIndex(Pm[i]),
                 InertiaDegree(Pm[i]),Minimum(Pm[i])>:i in [1..#Pm]];
        ramdeg:=&+[Degree(Pm[i]):i in [1..#Pm]|IsOdd(nm[i])];
        print "RELATIVE_SIGN_QUADRATIC_RAMIFICATION_DEGREE",ramdeg;
        grel:=Genus(Frel);
        gsigned:=2*grel-1+ramdeg div 2;
        print "RELATIVE_SUPPORT_GENUS",grel,
              "RELATIVE_SIGNED_GENUS_BY_RH",gsigned;
        assert ramdeg eq 2 and grel eq 5 and gsigned eq 10;
        print "CONTACT6_M612_RELATIVE3_M_DIVISOR_DONE";
        quit;
    end if;
    print "RELATIVE_SUPPORT_GENUS_START";
    tr0:=Cputime(); grel:=Genus(Frel);
    print "RELATIVE_SUPPORT_GENUS",grel,"SECONDS",Cputime(tr0);
    mrel:=Evaluate(mrec,wrel);
    msqrel:=IsSquare(mrel);
    print "RELATIVE_M_IS_SQUARE",msqrel;
    if do_signed and not msqrel then
        FrelL<Lrel>:=PolynomialRing(Frel);
        Frels<ellrel>:=ext<Frel|Lrel^2-mrel>;
        print "RELATIVE_SIGNED_GENUS_START";
        tr1:=Cputime(); grels:=Genus(Frels);
        print "RELATIVE_SIGNED_GENUS",grels,"SECONDS",Cputime(tr1);
        print "RELATIVE_SIGNED_QUADRATIC_RAMIFICATION_DEGREE",
              2*grels-4*grel+2;
    end if;
    print "CONTACT6_M612_RELATIVE3_EXACT_GENUS_DONE";
    quit;
end if;

Ku<u>:=FunctionField(Q);
t:=4*(u^2+u-6)/(u^2+6);
eP:=-(Q!25/3)*t^2/(t^4-25*t^2+Q!1250/3);
phi:=hom<K -> Ku | eP>;
Kuz<Z>:=PolynomialRing(Ku);
f12P:=Kuz![phi(Coefficient(f12,i)):i in [0..12]];
mrecP:=Kuz![phi(Coefficient(mrec,i)):i in [0..Degree(mrec)]];

print "CONTACT6_M612_P8_RELATIVE3_EXACT_GENUS";
print "P8_SUPPORT_POLYNOMIAL_DEGREE",Degree(f12P);
if discriminant_only then
    print "P8_PRIMITIVE_DISCRIMINANT_START";
    td:=Cputime(); dP:=Discriminant(f12P);
    print "P8_PRIMITIVE_DISCRIMINANT_DEGREES",
          <Degree(Numerator(dP)),Degree(Denominator(dP))>,
          "SECONDS",Cputime(td);
    print "P8_PRIMITIVE_DISCRIMINANT_NUMERATOR_FACTOR_START";
    tdn:=Cputime(); fdn:=Factorization(Numerator(dP));
    print "P8_PRIMITIVE_DISCRIMINANT_NUMERATOR_FACTORS",
          [<Degree(fe[1]),fe[2]>:fe in fdn],"SECONDS",Cputime(tdn);
    print "P8_PRIMITIVE_DISCRIMINANT_DENOMINATOR_FACTOR_START";
    tdd:=Cputime(); fdd:=Factorization(Denominator(dP));
    print "P8_PRIMITIVE_DISCRIMINANT_DENOMINATOR_FACTORS",
          [<Degree(fe[1]),fe[2]>:fe in fdd],"SECONDS",Cputime(tdd);
    print "CONTACT6_M612_P8_RELATIVE3_EXACT_DISCRIMINANT_DONE";
    quit;
end if;
print "P8_SUPPORT_IRREDUCIBILITY_START";
t0:=Cputime(); irr12:=IsIrreducible(f12P);
print "P8_SUPPORT_IRREDUCIBLE",irr12,"SECONDS",Cputime(t0);
assert irr12;

F12<w>:=ext<Ku|f12P>;
print "P8_SUPPORT_GENUS_START";
t1:=Cputime(); g12:=Genus(F12);
print "P8_SUPPORT_GENUS",g12,"SECONDS",Cputime(t1);

m12:=Evaluate(mrecP,w);
print "P8_SIGN_SQUARE_TEST_START";
t2:=Cputime(); msquare,mroot:=IsSquare(m12);
print "P8_M_IS_SQUARE",msquare,"SECONDS",Cputime(t2);
assert not msquare;

if do_signed then
    F12L<L>:=PolynomialRing(F12);
    F24<ell>:=ext<F12|L^2-m12>;
    print "P8_SIGNED_GENUS_START";
    t3:=Cputime(); g24:=Genus(F24);
    print "P8_SIGNED_GENUS",g24,"SECONDS",Cputime(t3);
    ramdeg:=2*g24-4*g12+2;
    print "P8_SIGNED_QUADRATIC_RAMIFICATION_DEGREE",ramdeg;
    assert ramdeg ge 0 and IsEven(ramdeg);
end if;
print "CONTACT6_M612_P8_RELATIVE3_EXACT_GENUS_DONE";
quit;
