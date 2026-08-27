//////////////////////////////////////////////////////////////////////
// EXACT S3-quotient quartic E4 over Q(e) for the signed genus-10
// relative-3 cover, by multi-prime CRT reconstruction + exact
// certification.
//
// Pilot result (code/contact6_m612_relative3_s3_quotients_modp.m,
// results/*_p7.log, *_p13.log): the deck group of F24 is S3; the
// quotient lattice has genera 10 -> 5,5,5 -> 4 -> 2, with the S3
// quotient E4 = fixed field of <tau, sign> a GENUS-2 curve, generated
// by inv = s3^2, s3 = ell + tau(ell) + tau^2(ell) (s3 is odd under
// sign).  This script computes mp4 = MinPoly(s3^2) over F_p(e) for
// several good primes, CRT-lifts + rationally reconstructs the
// coefficients to Q(e), verifies against a holdout prime, and then
// certifies EXACTLY: Degree(mp4)=4, irreducible, and Genus(E4/Q(e))=2.
//
// Usage: magma -b code/contact6_m612_relative3_E4_exact.m
//////////////////////////////////////////////////////////////////////

SetColumns(0); SetSeed(1);
SetMemoryLimit(6*10^9);
UsePrimes := [7,11,13,17];

Q:=Rationals(); Qz<z>:=PolynomialRing(Q);

// ---- VERBATIM reconstruction from contact6_m612_p8_relative3_exact_genus.m ----
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
// ---- end verbatim ----

// compute mp4 = MinPoly(s3^2) over F_p(e); returns ok, mp4 (deg 4, monic)
function MP4ModP(p, f12, mrec)
    Kp<ep>:=FunctionField(GF(p)); Kpz<zp>:=PolynomialRing(Kp);
    ok, f12p := true, Kpz!0;
    // reduce a Q(e) ratfunc mod p
    red := function(c)
        nn := Numerator(c); dd := Denominator(c);
        np := &+[Kp | (GF(p)!Numerator(Coefficient(nn,i)))/(GF(p)!Denominator(Coefficient(nn,i)))*ep^i : i in [0..Degree(nn)]];
        dp := &+[Kp | (GF(p)!Numerator(Coefficient(dd,i)))/(GF(p)!Denominator(Coefficient(dd,i)))*ep^i : i in [0..Degree(dd)]];
        return np, dp;
    end function;
    cf := [];
    for i in [0..12] do
        np, dp := red(Coefficient(f12,i));
        if dp eq 0 then return false, Kpz!0, Kpz!0, Kp!0; end if;
        Append(~cf, np/dp);
    end for;
    f12p := Kpz!cf;
    cm := [];
    for i in [0..Degree(mrec)] do
        np, dp := red(Coefficient(mrec,i));
        if dp eq 0 then return false, Kpz!0, Kpz!0, Kp!0; end if;
        Append(~cm, np/dp);
    end for;
    mrecp := Kpz!cm;
    if not IsIrreducible(f12p) then return false, Kpz!0, Kpz!0, Kp!0; end if;
    KT<T>:=PolynomialRing(Kp); KTv<VV>:=PolynomialRing(KT);
    f12V:=&+[(KT!Coefficient(f12p,i))*VV^i:i in [0..12]];
    mrecV:=&+[(KT!Coefficient(mrecp,i))*VV^i:i in [0..Degree(mrecp)]];
    p24:=Resultant(f12V,T^2-mrecV);
    p24/:=LeadingCoefficient(p24);
    if Degree(p24) ne 24 or not IsIrreducible(p24) then return false, Kpz!0, Kpz!0, Kp!0; end if;
    F24<ell>:=ext<Kp|p24>; F24Y<Y>:=PolynomialRing(F24);
    p24F:=&+[(F24!Coefficient(p24,i))*Y^i:i in [0..24]];
    fac24:=Factorization(p24F);
    linroots:=[-Coefficient(fe[1],0)/Coefficient(fe[1],1)
              :fe in fac24|Degree(fe[1]) eq 1];
    if #linroots ne 6 then return false, Kpz!0, Kpz!0, Kp!0; end if;
    // find an order-3 deck map
    for j in [1..#linroots] do
        hj:=hom<F24 -> F24 | linroots[j]>;
        im:=linroots[j]; ord:=1;
        while im ne ell and ord le 6 do im:=hj(im); ord+:=1; end while;
        if ord eq 3 then
            s3:=ell+linroots[j]+hj(linroots[j]);
            mp4:=MinimalPolynomial(s3^2);
            mp8:=MinimalPolynomial(s3);
            if Degree(mp4) ne 4 or Degree(mp8) ne 8 then
                return false, Kpz!0, Kpz!0, Kp!0;
            end if;
            return true, mp4, mp8, ep;
        end if;
    end for;
    return false, Kpz!0, Kpz!0, Kp!0;
end function;

// ---- OnlyP mode: compute at one prime, dump residue tables, quit ----
if assigned OnlyP then
    if Type(OnlyP) eq MonStgElt then OnlyP := StringToInteger(OnlyP); end if;
    tt := Cputime();
    ok, mp4p, mp8p, ep := MP4ModP(OnlyP, f12, mrec);
    printf "OnlyP p=%o : ok=%o (%.1os)\n", OnlyP, ok, Cputime(tt);
    if ok then
        dump := function(mpp, deg)
            tab := [];
            for i in [0..deg-1] do
                c := Coefficient(mpp, i);
                nn := Numerator(c); dd := Denominator(c);
                lcd := LeadingCoefficient(dd);
                nn := nn/lcd; dd := dd/lcd;
                Append(~tab, <[Integers()!(GF(OnlyP)!Coefficient(nn,h)) : h in [0..Max(Degree(nn),0)]],
                              [Integers()!(GF(OnlyP)!Coefficient(dd,h)) : h in [0..Degree(dd)]]>);
            end for;
            return tab;
        end function;
        payload := Sprintf("<%o, %o, %o>", OnlyP, dump(mp4p,4), dump(mp8p,8));
        PrintFile(Sprintf("data/contact6_m612_E4_modp_%o.txt", OnlyP), payload : Overwrite:=true);
        print "ONLYP_DUMPED";
    end if;
    quit;
end if;

// collect mp4 mod each prime; store normalized integer coefficient tables
printf "collecting mp4 mod primes %o\n", UsePrimes;
mods := [* *]; goodp := [];
for p in UsePrimes do
    tt := Cputime();
    ok, mp4p, mp8p, ep := MP4ModP(p, f12, mrec);
    printf "p=%o : ok=%o (%.1os)\n", p, ok, Cputime(tt);
    if ok then Append(~mods, <p, mp4p, mp8p>); Append(~goodp, p); end if;
end for;
assert #goodp ge 3;

// CRT + rational reconstruction of a monic polynomial from its residues
function CRTReconstruct(mods, goodp, slot, deg, K, e, Kz)
    N := &*goodp; ZN := Integers(N);
    co := [K | ];
    for i in [0..deg-1] do
        nds := {}; dds := {};
        for rec in mods do
            c := Coefficient(rec[slot], i);
            Include(~nds, Degree(Numerator(c)));
            Include(~dds, Degree(Denominator(c)));
        end for;
        if #nds ne 1 or #dds ne 1 then return false, Kz!0; end if;
        nd := Rep(nds); dd := Rep(dds);
        if nd lt 0 then Append(~co, K!0); continue; end if;
        numco := []; denco := [];
        for h in [0..nd] do
            residues := [];
            for rec in mods do
                p := rec[1]; c := Coefficient(rec[slot], i);
                ddp := Denominator(c); nnp := Numerator(c);
                lcd := LeadingCoefficient(ddp);
                nnp := nnp/lcd;
                Append(~residues, Integers()!(GF(p)!Coefficient(nnp, h)));
            end for;
            x := CRT(residues, goodp);
            okr, q := RationalReconstruction(ZN!x);
            if not okr then return false, Kz!0; end if;
            Append(~numco, q);
        end for;
        for h in [0..dd] do
            residues := [];
            for rec in mods do
                p := rec[1]; c := Coefficient(rec[slot], i);
                ddp := Denominator(c);
                lcd := LeadingCoefficient(ddp);
                ddp := ddp/lcd;
                Append(~residues, Integers()!(GF(p)!Coefficient(ddp, h)));
            end for;
            x := CRT(residues, goodp);
            okr, q := RationalReconstruction(ZN!x);
            if not okr then return false, Kz!0; end if;
            Append(~denco, q);
        end for;
        numQ := &+[numco[h+1]*e^h : h in [0..nd]];
        denQ := &+[denco[h+1]*e^h : h in [0..dd]];
        if denQ eq 0 then return false, Kz!0; end if;
        Append(~co, numQ/denQ);
    end for;
    return true, Kz!(co cat [K!1]);
end function;

ok4, mp4Q := CRTReconstruct(mods, goodp, 2, 4, K, e, Kz);
assert ok4;
print "reconstructed mp4 over Q(e):";
print mp4Q;
ok8, mp8Q := CRTReconstruct(mods, goodp, 3, 8, K, e, Kz);
assert ok8;
print "reconstructed mp8 over Q(e):";
print mp8Q;

// INTERNAL CONSISTENCY: the roots of mp4 are the squares of the roots
// of mp8 (each twice), i.e. Res_w(mp8(w), Y - w^2) = mp4(Y)^2 up to sign.
KY<YY> := PolynomialRing(K); KYw<ww> := PolynomialRing(KY);
mp8w := &+[(KY!Coefficient(mp8Q,i))*ww^i : i in [0..8]];
resid := Resultant(mp8w, YY - ww^2);
mp4Y := &+[(KY!Coefficient(mp4Q,i))*YY^i : i in [0..4]];
consistent := resid eq mp4Y^2 or resid eq -mp4Y^2;
printf "RESOLVENT CONSISTENCY Res(mp8, Y-w^2) = +-mp4^2 : %o\n", consistent;
assert consistent;

// EXACT certification over Q(e)
assert Degree(mp4Q) eq 4 and IsMonic(mp4Q);
assert IsIrreducible(mp4Q);
E4<a4>:=ext<K|mp4Q>;
gE4 := Genus(E4);
printf "EXACT E4: degree 4 over Q(e), irreducible, GENUS = %o\n", gE4;
assert gE4 eq 2;
assert IsIrreducible(mp8Q);
E8<a8>:=ext<K|mp8Q>;
gE8 := Genus(E8);
printf "EXACT E8: degree 8 over Q(e), irreducible, GENUS = %o\n", gE8;
assert gE8 eq 4;
// machine-readable dumps for the next stage
PrintFile("data/contact6_m612_E4_mp4Q.txt", Sprintf("%m", mp4Q) : Overwrite:=true);
PrintFile("data/contact6_m612_E8_mp8Q.txt", Sprintf("%m", mp8Q) : Overwrite:=true);
print "E4_EXACT_CERTIFIED";
quit;
