//////////////////////////////////////////////////////////////////////
// S3-quotient genera of the signed genus-10 relative-3 cover, MOD p.
//
// Good-reduction pilot for the automorphism_only mode of
// code/contact6_m612_p8_relative3_exact_genus.m (whose exact Q(e)
// factorization step is slow).  Rebuilds f12, m(v) exactly over Q(e)
// (fast; same FixedRelative12AndM + reconstruction as the parent script),
// reduces mod p into F_p(e), forms the signed degree-24 primitive
// polynomial, factors it over its own root field (cheap in char p),
// finds the C3 deck maps, and computes the genera of the C3- and
// S3-invariant subfields E8, E4.  At good p these equal the
// characteristic-zero quotient genera (the technique used for 73/145).
//
// Usage: magma -b p:=13 code/contact6_m612_relative3_s3_quotients_modp.m
//////////////////////////////////////////////////////////////////////

SetColumns(0); SetSeed(1);
if not assigned p then p := 13; elif Type(p) eq MonStgElt then p := StringToInteger(p); end if;
SetMemoryLimit(4*10^9);

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
// ---- end verbatim reconstruction ----

print "S3_QUOTIENTS_MODP p =", p;

// reduce coefficients of f12, mrec into F_p(e)
Kp<ep>:=FunctionField(GF(p)); Kpz<zp>:=PolynomialRing(Kp);
function RedC(c, Kp, ep)   // c in Q(e) -> F_p(e); fails if p | denominator content
    nn := Numerator(c); dd := Denominator(c);
    np := &+[Kp | (GF(p)!Numerator(Coefficient(nn,i)))/(GF(p)!Denominator(Coefficient(nn,i)))*ep^i : i in [0..Degree(nn)]];
    dp := &+[Kp | (GF(p)!Numerator(Coefficient(dd,i)))/(GF(p)!Denominator(Coefficient(dd,i)))*ep^i : i in [0..Degree(dd)]];
    assert dp ne 0;
    return np/dp;
end function;
f12p := Kpz![RedC(Coefficient(f12,i),Kp,ep) : i in [0..12]];
mrecp := Kpz![RedC(Coefficient(mrec,i),Kp,ep) : i in [0..Degree(mrec)]];
assert IsIrreducible(f12p);
print "f12 irreducible mod p: true";

// signed degree-24 primitive polynomial over F_p(e)
KT<T>:=PolynomialRing(Kp); KTv<VV>:=PolynomialRing(KT);
f12V:=&+[(KT!Coefficient(f12p,i))*VV^i:i in [0..12]];
mrecV:=&+[(KT!Coefficient(mrecp,i))*VV^i:i in [0..Degree(mrecp)]];
tp:=Cputime(); p24:=Resultant(f12V,T^2-mrecV);
p24/:=LeadingCoefficient(p24);
printf "p24 degree %o, resultant %.2os\n", Degree(p24), Cputime(tp);
assert Degree(p24) eq 24;
assert IsIrreducible(p24);
print "p24 irreducible mod p: true (connected signed cover)";

F24<ell>:=ext<Kp|p24>; F24Y<Y>:=PolynomialRing(F24);
p24F:=&+[(F24!Coefficient(p24,i))*Y^i:i in [0..24]];
tf:=Cputime(); fac24:=Factorization(p24F);
printf "internal factorization %.2os : degrees %o\n", Cputime(tf),
    [<Degree(fe[1]),fe[2]>:fe in fac24];
linroots:=[-Coefficient(fe[1],0)/Coefficient(fe[1],1)
          :fe in fac24|Degree(fe[1]) eq 1];
printf "deck (linear) roots: %o\n", #linroots;

// classify all deck maps by exact order and compute EVERY subgroup quotient
orders:=[];
for j in [1..#linroots] do
    hj:=hom<F24 -> F24 | linroots[j]>;
    im:=linroots[j]; ord:=1;
    while im ne ell and ord le 6 do im:=hj(im); ord+:=1; end while;
    Append(~orders, ord);
    printf "deck map %o: exact order %o (root = -ell? %o)\n", j, ord, linroots[j] eq -ell;
end for;
printf "deck group orders: %o (S3 iff [1,2,2,2,3,3] as multiset)\n", Sort(orders);

// helper: genus of the fixed field generated by a given invariant
function QuotGenus(inv, expectdeg, Kp)
    mp:=MinimalPolynomial(inv);
    if Degree(mp) ne expectdeg then return -Degree(mp); end if;  // signal wrong degree
    E:=ext<Kp|mp>;
    return Genus(E);
end function;

// C3 quotient
j3:=[j:j in [1..#linroots]|orders[j] eq 3][1];
h3:=hom<F24 -> F24 | linroots[j3]>;
s3:=ell+linroots[j3]+h3(linroots[j3]);
g8:=QuotGenus(s3, 8, Kp);
printf "C3_QUOTIENT_GENUS (E8, deg 8) = %o\n", g8;

// S3 quotient: s3 is odd under sign (check), so use s3^2 (or product fallback)
hneg:=hom<F24 -> F24 | -ell>;
printf "hneg(s3) = -s3? %o\n", hneg(s3) eq -s3;
invS3 := hneg(s3) eq -s3 select s3^2 else s3*hneg(s3);
g4:=QuotGenus(invS3, 4, Kp);
printf "S3_QUOTIENT_GENUS (E4, deg 4) = %o\n", g4;

// the three C2 quotients (degree-12 subfields)
for j in [j:j in [1..#linroots]|orders[j] eq 2] do
    hj:=hom<F24 -> F24 | linroots[j]>;
    isSigma := linroots[j] eq -ell;
    inv2 := ell + linroots[j];
    if inv2 eq 0 or MinimalPolynomial(inv2) eq Parent(MinimalPolynomial(inv2)).1 then
        inv2 := ell * linroots[j];   // fallback: product invariant
    end if;
    mp:=MinimalPolynomial(inv2);
    if Degree(mp) eq 12 then
        E12:=ext<Kp|mp>;
        printf "C2_QUOTIENT_GENUS (deck %o%o, deg 12) = %o\n",
            j, isSigma select " = sign flip -> F12" else "", Genus(E12);
    else
        printf "C2 quotient deck %o: invariant degree %o (not primitive; skipped)\n",
            j, Degree(mp);
    end if;
end for;

// reference
F12p<w>:=ext<Kp|f12p>;
printf "SUPPORT_GENUS (F12, deg 12) = %o\n", Genus(F12p);
printf "SIGNED_GENUS (F24, deg 24) = %o\n", Genus(F24);
print "S3_QUOTIENTS_MODP_DONE";
quit;
