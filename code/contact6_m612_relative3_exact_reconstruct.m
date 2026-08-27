//////////////////////////////////////////////////////////////////////
// Reconstruct and certify the exact degree-12 orthogonal support cover
// over Q(e) from small fixed-e fibers.  Modular computations provide
// small interpolation bounds; direct symbolic identities in the exact
// quotient certify the reconstructed polynomial and recovery maps.
//
// The resulting component is then pulled back by the rational P8
// parameter e=e(u).  Adjoining L with L^2=M gives the actual-point cover.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned do_primitive24 then do_primitive24:=true;
elif Type(do_primitive24) eq MonStgElt then
    do_primitive24:=do_primitive24 in {"true","True","1","yes"};
end if;
if not assigned print_maps then print_maps:=false;
elif Type(print_maps) eq MonStgElt then
    print_maps:=print_maps in {"true","True","1","yes"};
end if;
Q:=Rationals(); Qz<z>:=PolynomialRing(Q);

function FixedRelative12(ee)
    R0<M0,U0,v0>:=PolynomialRing(Q,3,"grevlex");
    N0:=(3*U0+6*M0)/2;
    Rc0:=(3*U0^2+3*v0^2+(2/ee-15)*M0-N0^2)/2;
    Fs0:=[2*v0^3+2*N0*Rc0-U0^3-6*U0*v0^2-22*M0,
          Rc0^2+2*N0*v0^3-3*U0^2*v0^2-3*v0^4
              -(1/ee^2-15)*M0,
          2*Rc0*v0^3-3*U0*v0^4-(2/ee+6)*M0];
    I0:=Saturation(ideal<R0|Fs0>,ideal<R0|M0>);
    if Dimension(I0) ne 0 or Dimension(quo<R0|I0>) ne 40 then
        return false,Qz!0,Qz!0,Qz!0;
    end if;
    G0:=GroebnerBasis(ChangeOrder(I0,"lex"));
    univ:=[g:g in G0|Degree(g,1) eq 0 and Degree(g,2) eq 0 and
                       Degree(g,3) gt 0];
    if #univ ne 1 then return false,Qz!0,Qz!0,Qz!0; end if;
    rv:=Qz!UnivariatePolynomial(univ[1]);
    rv/:=LeadingCoefficient(rv);
    f12s:=[fe[1]:fe in Factorization(rv)|Degree(fe[1]) eq 12 and fe[2] eq 1];
    if #f12s ne 1 then return false,Qz!0,Qz!0,Qz!0; end if;
    f12:=f12s[1]/LeadingCoefficient(f12s[1]);
    gm:=[g:g in G0|Degree(g,1) eq 1 and Degree(g,2) eq 0][1];
    gu:=[g:g in G0|Degree(g,1) eq 0 and Degree(g,2) eq 1][1];
    cm:=Q!Coefficient(gm,1,1); cu:=Q!Coefficient(gu,2,1);
    mrec:=(-Qz!UnivariatePolynomial(Evaluate(gm,[R0!0,R0!0,v0]))/cm)
          mod f12;
    urec:=(-Qz!UnivariatePolynomial(Evaluate(gu,[R0!0,R0!0,v0]))/cu)
          mod f12;
    return true,Qz!f12,Qz!mrec,Qz!urec;
end function;

// Bounds (numerator degree, denominator degree) stabilized across good
// modular factors; the exact identities below certify the reconstruction.
bounds:=[<4,4>,<3,4>,<4,4>,<3,3>,<3,3>,<3,3>,<2,3>,
         <3,3>,<3,3>,<2,2>,<1,1>,<0,1>,<0,0>];
sample_e:=[Q!j:j in [1..20]];
samples:=[]; msamples:=[]; usamples:=[];
for ee in sample_e do
    ok,f,mf,uf:=FixedRelative12(ee);
    assert ok;
    Append(~samples,f);
    Append(~msamples,mf);
    Append(~usamples,uf);
end for;

K<e>:=FunctionField(Q); Kz<zz>:=PolynomialRing(K);
function ReconstructPolynomial(poly_samples,degree,bds,Q,K,e,Kz)
  coeffs:=[];
  for i in [0..degree] do
    nd:=bds[i+1][1]; dd:=bds[i+1][2];
    ncols:=nd+dd+2;
    rows:=[];
    for j in [1..#sample_e] do
        x:=sample_e[j]; y:=Coefficient(poly_samples[j],i);
        Append(~rows,[x^h:h in [0..nd]] cat
                      [-y*x^h:h in [0..dd]]);
    end for;
    A:=Matrix(Q,#rows,ncols,&cat rows);
    W:=Nullspace(Transpose(A));
    assert Dimension(W) eq 1;
    w:=Basis(W)[1];
    num:=&+[w[h+1]*e^h:h in [0..nd]];
    den:=&+[w[nd+2+h]*e^h:h in [0..dd]];
    assert den ne 0;
    c:=num/den;
    assert &and[Evaluate(Numerator(c),sample_e[j])/
                 Evaluate(Denominator(c),sample_e[j]) eq
                 Coefficient(poly_samples[j],i):j in [1..#sample_e]];
    Append(~coeffs,c);
  end for;
  return Kz!coeffs;
end function;

f12:=ReconstructPolynomial(samples,12,bounds,Q,K,e,Kz);
assert Degree(f12) eq 12 and IsMonic(f12);

mbounds:=[<6,3>,<6,3>,<6,3>,<6,3>,<7,4>,<6,4>,
          <7,5>,<8,5>,<8,5>,<8,5>,<8,5>,<8,5>];
ubounds:=[<7,4>,<7,4>,<6,3>,<7,4>,<7,4>,<6,4>,
          <7,5>,<8,5>,<8,5>,<9,6>,<9,6>,<9,6>];
mrec:=ReconstructPolynomial(msamples,11,mbounds,Q,K,e,Kz) mod f12;
urec:=ReconstructPolynomial(usamples,11,ubounds,Q,K,e,Kz) mod f12;

print "CONTACT6_M612_RELATIVE3_EXACT_RECONSTRUCT";
print "SAMPLE_COUNT",#samples,"SAMPLE_E",sample_e;
print "RELATIVE12_COEFFICIENT_SHAPES",
      [<Degree(Numerator(Coefficient(f12,i))),
         Degree(Denominator(Coefficient(f12,i)))>:i in [0..12]];
print "RELATIVE12_POLYNOMIAL",f12;

// Direct exact symbolic certificate: the recovered maps satisfy every
// contact equation in K[z]/(f12).  This avoids a costly characteristic-zero
// primary decomposition or FGLM conversion.
nrec:=(3*urec+6*mrec)/2;
rrec:=(3*urec^2+3*zz^2+(2/e-15)*mrec-nrec^2)/2;
F3rec:=2*zz^3+2*nrec*rrec-urec^3-6*urec*zz^2-22*mrec;
F2rec:=rrec^2+2*nrec*zz^3-3*urec^2*zz^2-3*zz^4
       -(1/e^2-15)*mrec;
F1rec:=2*rrec*zz^3-3*urec*zz^4-(2/e+6)*mrec;
assert &and[(g mod f12) eq 0:g in [F3rec,F2rec,F1rec]];
assert GCD(f12,mrec) eq 1 and Coefficient(f12,0) ne 0;
drec:=(urec^2-4*zz^2) mod f12;
assert GCD(f12,drec) eq 1;
print "EXACT_COMPONENT_LENGTH",12,"DIRECT_CONTACT_IDENTITIES",true;
print "M_RECOVERY_DEGREE",Degree(mrec),"U_RECOVERY_DEGREE",Degree(urec);
if print_maps then
    print "M_RECOVERY",mrec;
    print "U_RECOVERY",urec;
end if;

assert IsIrreducible(f12);
F12<w>:=ext<K|f12>;
melt:=Evaluate(mrec,w);
msquare,mroot:=IsSquare(melt);
print "M_IS_SQUARE_ON_EXACT_RELATIVE12",msquare;
assert not msquare;

KTLplaceholder<TLplaceholder>:=PolynomialRing(K);
p24:=KTLplaceholder!0;
if do_primitive24 then
    // The norm of T^2-M(v) is a candidate minimal polynomial for L.
    // Degree 24 and irreducibility prove that L itself is primitive.
    KT<T>:=PolynomialRing(K); KTv<VV>:=PolynomialRing(KT);
    f12V:=&+[(KT!Coefficient(f12,i))*VV^i:i in [0..12]];
    mrecV:=&+[(KT!Coefficient(mrec,i))*VV^i:i in [0..Degree(mrec)]];
    p24:=Resultant(f12V,T^2-mrecV);
    p24/:=LeadingCoefficient(p24);
    print "L_PRIMITIVE_POLYNOMIAL_DEGREE",Degree(p24);
    print "L_PRIMITIVE_COEFFICIENT_SHAPES",
          [<Degree(Numerator(Coefficient(p24,i))),
             Degree(Denominator(Coefficient(p24,i)))>:i in [0..Degree(p24)]];
    assert Degree(p24) eq 24 and IsIrreducible(p24);

    // Recover v by the degree-one gcd over K(L); all remaining contact
    // coordinates then follow from the maps already certified above.
    F24<ell>:=ext<K|p24>; F24V<X>:=PolynomialRing(F24);
    f24:=&+[(F24!Coefficient(f12,i))*X^i:i in [0..12]];
    m24:=&+[(F24!Coefficient(mrec,i))*X^i:i in [0..Degree(mrec)]];
    gv:=GCD(f24,ell^2-m24);
    assert Degree(gv) eq 1;
    v_of_L:=-Coefficient(gv,0)/Coefficient(gv,1);
    u_of_L:=Evaluate(F24V!urec,v_of_L);
    print "L_IS_PRIMITIVE",true,
          "V_RECOVERY_BASIS_LENGTH",#Eltseq(v_of_L),
          "U_RECOVERY_BASIS_LENGTH",#Eltseq(u_of_L);
end if;

// Exact P8 base change.  Irreducibility and nonsquareness are also
// certified by the good F_7(u) reduction in the modular companion script.
Ku<u>:=FunctionField(Q);
t:=4*(u^2+u-6)/(u^2+6);
eP:=-(Q!25/3)*t^2/(t^4-25*t^2+Q!1250/3);
phi:=hom<K -> Ku | eP>;
Kuz<Z>:=PolynomialRing(Ku);
f12P:=Kuz![phi(Coefficient(f12,i)):i in [0..12]];
mrecP:=Kuz![phi(Coefficient(mrec,i)):i in [0..Degree(mrec)]];
urecP:=Kuz![phi(Coefficient(urec,i)):i in [0..Degree(urec)]];
print "P8_RELATIVE12_DEGREE",Degree(f12P);
print "P8_RELATIVE12_COEFFICIENT_SHAPES",
      [<Degree(Numerator(Coefficient(f12P,i))),
         Degree(Denominator(Coefficient(f12P,i)))>:i in [0..12]];
if do_primitive24 then
    p24P:=Kuz![phi(Coefficient(p24,i)):i in [0..24]];
    print "P8_L_PRIMITIVE_POLYNOMIAL_DEGREE",Degree(p24P);
    print "P8_L_PRIMITIVE_COEFFICIENT_SHAPES",
          [<Degree(Numerator(Coefficient(p24P,i))),
             Degree(Denominator(Coefficient(p24P,i)))>:i in [0..24]];
end if;
print "ACTUAL_POINT_COVER","f12P(v)=0, L^2=mrecP(v)",
      "TOTAL_DEGREE",24;

print "CONTACT6_M612_RELATIVE3_EXACT_RECONSTRUCT_DONE";
quit;
