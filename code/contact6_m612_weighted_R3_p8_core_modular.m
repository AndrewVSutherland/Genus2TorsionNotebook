//////////////////////////////////////////////////////////////////////
// Modular generic geometry of the cubic-contact [3,3] core pulled back
// to the rational P8 endpoint R3 component.
//
// The P8 normalization parameter is u, with
//   t=4(u^2+u-6)/(u^2+6),
//   e=-(25/3)t^2/(t^4-25t^2+1250/3), a=1/e, b=0.
//
// cover="core" works in (M,U,v), M=L^2.
// cover="square" substitutes M=L^2 and works in (L,U,v), resolving in v.
// cover="squareL" uses order (U,v,L), resolving in L.
// The contact-open product is M*v*(U^2-4v^2), respectively
// L*v*(U^2-4v^2).
//
// Usage:
//   timeout 300s magma -b p:=11 cover:=core code/..._modular.m
//   timeout 300s magma -b p:=11 cover:=square code/..._modular.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);
if not assigned p then p:=11; elif Type(p) eq MonStgElt then p:=StringToInteger(p); end if;
if not assigned cover then cover:="core"; end if;
if not assigned genus_max_degree then genus_max_degree:=0;
elif Type(genus_max_degree) eq MonStgElt then
    genus_max_degree:=StringToInteger(genus_max_degree);
end if;
if not assigned square_genus_max_degree then square_genus_max_degree:=0;
elif Type(square_genus_max_degree) eq MonStgElt then
    square_genus_max_degree:=StringToInteger(square_genus_max_degree);
end if;
if not assigned do_branch then do_branch:=false;
elif Type(do_branch) eq MonStgElt then
    do_branch:=do_branch in {"true","True","1","yes","Yes"};
end if;
require IsPrime(p) and p notin {2,3,5}: "use a prime away from 2,3,5";

k:=GF(p);
K<u>:=FunctionField(k);
t:=4*(u^2+u-6)/(u^2+6);
e:=-(K!25/3)*t^2/(t^4-25*t^2+K!1250/3);
a:=1/e;

print "CONTACT6_M612_WEIGHTED_R3_P8_CORE_MODULAR";
print "PRIME",p,"COVER",cover;
print "PARAMETER_DEGREES_t_e",
      <Degree(Numerator(t)),Degree(Denominator(t))>,
      <Degree(Numerator(e)),Degree(Denominator(e))>;

if cover eq "core" then
    P<z,M,U,v>:=PolynomialRing(K,4,"grevlex");
    MM:=M;
    open:=M*v*(U^2-4*v^2);
elif cover eq "coreM" then
    P<z,U,v,M>:=PolynomialRing(K,4,"grevlex");
    MM:=M;
    open:=M*v*(U^2-4*v^2);
elif cover eq "square" then
    P<z,L,U,v>:=PolynomialRing(K,4,"grevlex");
    MM:=L^2;
    open:=L*v*(U^2-4*v^2);
elif cover eq "squareL" then
    P<z,U,v,L>:=PolynomialRing(K,4,"grevlex");
    MM:=L^2;
    open:=L*v*(U^2-4*v^2);
else
    error "cover must be core, coreM, square, or squareL";
end if;

c1:=2*a+6; c2:=a^2-15; c3:=K!22; c4:=2*a-15; c5:=K!6;
B3:=c5*MM+3*U;
Delta3:=4*c4*MM+12*(U^2+v^2)-B3^2;
F3:=B3*Delta3+16*v^3-8*c3*MM-8*U^3-48*U*v^2;
F2:=Delta3^2+64*B3*v^3-64*c2*MM-192*(U^2*v^2+v^4);
F1:=Delta3*v^3-4*c1*MM-12*U*v^4;
I:=ideal<P|F3,F2,F1,z*open-1>;

print "GREVLEX_START";
time G0:=GroebnerBasis(I);
dim,pars:=Dimension(I);
deg:=Dimension(quo<P|I>);
print "GREVLEX_BASIS_SIZE",#G0;
print "DIMENSION",dim,"PARAMETERS",pars,"AFFINE_DEGREE",deg;
print "GREVLEX_SHAPES",
      [<TotalDegree(g),[Degree(g,j):j in [1..4]],#Terms(g)>:g in G0];

print "LEX_START";
time IL:=ChangeOrder(I,"lex");
time GL:=GroebnerBasis(IL);
print "LEX_BASIS_SIZE",#GL;
univ:=[g:g in GL|&and[Degree(g,j) eq 0:j in [1..3]] and Degree(g,4) gt 0];
resolve_label:=cover eq "squareL" select "L" else
               (cover eq "coreM" select "M" else "v");
print "UNIVARIATE_COUNT",resolve_label,#univ;
if #univ gt 0 then
    R:=UnivariatePolynomial(univ[#univ]);
    R/:=LeadingCoefficient(R);
    print "RESOLVENT",resolve_label,"DEGREE",Degree(R);
    time fac:=Factorization(R);
    print "RESOLVENT_FACTOR_DEGREES",resolve_label,
          [<Degree(q[1]),q[2]>:q in fac];
    for qi in [1..#fac] do
        q:=fac[qi][1];
        if do_branch then
            discq:=Discriminant(q);
            print "FACTOR",qi,"DEGREE",Degree(q),
                  "DISCRIMINANT_NUM_DEGREE",Degree(Numerator(discq)),
                  "DISCRIMINANT_DEN_DEGREE",Degree(Denominator(discq));
            print "FACTOR",qi,"DISCRIMINANT_NUM_FACTOR_DEGREES",
                  [<Degree(ff[1]),ff[2]>:ff in Factorization(Numerator(discq))];
            print "FACTOR",qi,"DISCRIMINANT_DEN_FACTOR_DEGREES",
                  [<Degree(ff[1]),ff[2]>:ff in Factorization(Denominator(discq))];
        end if;
        if genus_max_degree gt 0 and Degree(q) le genus_max_degree then
            print "FACTOR",qi,"FUNCTION_FIELD_GENUS_START";
            time Fq:=FunctionField(q:Check:=true);
            time gq:=Genus(Fq);
            print "FACTOR",qi,"FUNCTION_FIELD_GENUS",gq;
        end if;
        if square_genus_max_degree gt 0 and Degree(q) le square_genus_max_degree then
            if not assigned Fq then Fq:=FunctionField(q:Check:=true); end if;
            is_sq_M:=IsSquare(Fq.1);
            print "FACTOR",qi,"M_IS_SQUARE_IN_SUPPORT_FIELD",is_sq_M;
            if not is_sq_M then
                PFq<ellq>:=PolynomialRing(Fq);
                print "FACTOR",qi,"SQUARE_EXTENSION_GENUS_START";
                time EFq:=ext<Fq|ellq^2-Fq.1>;
                time egq:=Genus(EFq);
                print "FACTOR",qi,"SQUARE_EXTENSION_GENUS",egq;
            end if;
        end if;
    end for;
    if cover eq "coreM" then
        KL<ell>:=PolynomialRing(K);
        RL:=Evaluate(R,ell^2);
        print "SQUARE_PULLBACK_DEGREE",Degree(RL);
        time facL:=Factorization(RL);
        print "SQUARE_PULLBACK_FACTOR_DEGREES",
              [<Degree(q[1]),q[2]>:q in facL];
    end if;
end if;
print "CONTACT6_M612_WEIGHTED_R3_P8_CORE_MODULAR_DONE";
quit;
