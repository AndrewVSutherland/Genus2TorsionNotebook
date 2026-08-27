//////////////////////////////////////////////////////////////////////
// First local analysis of the ten p=17 parameter disks retained by the
// P8 extra-3 sieve.
//
// Disk types:
//   e=0 endpoint: tau=2,-3;
//   a=-2 common-root collision: tau=1,3,11,15;
//   DC=0 double-root collision: tau=0,8,12,infinity.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q:=Rationals(); Z:=Integers(); p:=17; k:=GF(p);
QPs<sQ>:=PolynomialRing(Q); kPs<s>:=PolynomialRing(k);
FK:=FieldOfFractions(kPs);

function ReducePolynomial(f)
    ans:=kPs!0;
    for i in [0..Degree(f)] do
        c:=Coefficient(f,i);
        assert Denominator(c) mod p ne 0;
        ans+:=(k!Numerator(c)/k!Denominator(c))*s^i;
    end for;
    return ans;
end function;

function EDataAffine(u)
    A:=u^2+u-k!6; B:=u^2+k!6;
    D:=k!768*A^4-k!1200*A^2*B^2+k!1250*B^4;
    if B eq 0 then return false,k!0,k!0,"P8_parameter_pole"; end if;
    if D eq 0 then return false,k!0,k!0,"e_pole"; end if;
    e:=-k!400*A^2*B^2/D;
    if e eq 0 then return false,e,k!0,"e_zero"; end if;
    a:=1/e;
    if a+2 eq 0 then return false,e,a,"a_plus_2_common_root"; end if;
    if a+3 eq 0 then return false,e,a,"a_plus_3_degree_boundary"; end if;
    if a^2-6*a-15 eq 0 then return false,e,a,"DB_double_root"; end if;
    if -8*a-15 eq 0 then return false,e,a,"DC_double_root"; end if;
    return true,e,a,"good";
end function;

function EMapAffine(u)
    A:=u^2+u-k!6; B:=u^2+k!6;
    D:=k!768*A^4-k!1200*A^2*B^2+k!1250*B^4;
    return -k!400*A^2*B^2/D;
end function;

badkeys:=[0,1,2,3,8,11,12,14,15];
classification:=[];
for key in badkeys do
    ok,e,a,reason:=EDataAffine(k!key);
    Append(~classification,<key,Z!e,Z!a,reason>);
end for;
einf:=-k!200/k!409; ainf:=1/einf;
assert -8*ainf-15 eq 0;
Append(~classification,<17,Z!einf,Z!ainf,"DC_double_root">);

// Derivatives of e(tau), including reciprocal coordinate w=1/tau at
// infinity.  Nonzero values certify that every collision disk is etale
// over the corresponding e-disk.
function EDerivativeAffine(u)
    A:=u^2+u-k!6; B:=u^2+k!6;
    Ap:=2*u+1; Bp:=2*u;
    Num:=-k!400*A^2*B^2;
    Den:=k!768*A^4-k!1200*A^2*B^2+k!1250*B^4;
    Np:=-k!400*(2*A*Ap*B^2+2*A^2*B*Bp);
    Dp:=k!768*4*A^3*Ap
        -k!1200*(2*A*Ap*B^2+2*A^2*B*Bp)
        +k!1250*4*B^3*Bp;
    return (Np*Den-Num*Dp)/Den^2;
end function;
derivatives:=[<key,Z!EDerivativeAffine(k!key)>:
              key in [0,1,3,8,11,12,15]];
// At infinity use w=1/tau: Abar=1+w-6w^2, Bbar=1+6w^2.
A0:=k!1; B0:=k!1; Ap0:=k!1; Bp0:=k!0;
N0:=-k!400*A0^2*B0^2;
D0:=k!768*A0^4-k!1200*A0^2*B0^2+k!1250*B0^4;
Np0:=-k!400*(2*A0*Ap0*B0^2+2*A0^2*B0*Bp0);
Dp0:=k!768*4*A0^3*Ap0
     -k!1200*(2*A0*Ap0*B0^2+2*A0^2*B0*Bp0)
     +k!1250*4*B0^3*Bp0;
Append(~derivatives,<17,Z!((Np0*D0-N0*Dp0)/D0^2)>);
assert &and[rec[2] ne 0:rec in derivatives];

// Endpoint e=0 expansions and E3 initial system.
Ebars:=AssociativeArray();
for rec in [<"ZERO_2",Q!2>,<"ZERO_MINUS3",Q!-3>] do
    tau:=QPs!(rec[2]+p*sQ);
    AA:=tau^2+tau-6; BB:=tau^2+6;
    enum:=-400*AA^2*BB^2;
    eden:=768*AA^4-1200*AA^2*BB^2+1250*BB^4;
    enump2:=QPs!(enum/p^2); assert p^2*enump2 eq enum;
    Ebars[rec[1]]:=FK!ReducePolynomial(enump2)/FK!ReducePolynomial(eden);
end for;

// Efficient enumeration of E3 using E,L,U,V and the first two equations
// to recover N,R.
e3_by_E:=AssociativeArray();
for EE in k do
    sols:=[];
    if EE ne 0 then
        for LL in k do
            if LL eq 0 then continue; end if;
            for UU in k do for VV in k do
                if VV eq 0 or UU^2-4*VV^2 eq 0 then continue; end if;
                NN:=(3*UU+6*LL^2)/2;
                RR:=(2*LL^2/EE-NN^2+3*UU^2+3*VV^2)/2;
                g3:=2*VV^3+2*NN*RR-UU^3-6*UU*VV^2;
                g2:=RR^2+2*NN*VV^3-3*UU^2*VV^2-3*VV^4;
                g1:=2*RR*VV^3-3*UU*VV^4;
                if g1 ne 0 or g2 ne 0 or g3 ne 0 then continue; end if;
                // The exact E3 fixed-variable Jacobian has rank five.
                S<E,L,U,V,N,R>:=PolynomialRing(k,6);
                geqs:=[2*N-3*U-6*L^2,
                    E*(N^2+2*R-3*U^2-3*V^2)-2*L^2,
                    2*V^3+2*N*R-U^3-6*U*V^2,
                    R^2+2*N*V^3-3*U^2*V^2-3*V^4,
                    2*R*V^3-3*U*V^4];
                pt:=[EE,LL,UU,VV,NN,RR];
                Jfix:=Matrix(k,5,5,
                    &cat[[Evaluate(Derivative(geqs[i],j),pt):j in [2..6]]
                          :i in [1..5]]);
                Append(~sols,<Z!LL,Z!UU,Z!VV,Z!NN,Z!RR,Rank(Jfix)>);
            end for; end for;
        end for;
    end if;
    e3_by_E[Z!EE]:=sols;
end for;

// Collision fibers.  Enumerate signed contact classes and retain only
// support-disjoint points.  Rank five in the five contact variables gives
// a fixed-parameter Hensel branch in every punctured disk.
function CollisionData(e0)
    P<x>:=PolynomialRing(k); a:=1/e0;
    h6:=1+a*x+x^3; f:=h6^2-(x-1)^6;
    S<L,U,V,N,R>:=PolynomialRing(k,5);
    c1:=2*a+6; c2:=a^2-15; c3:=k!22; c4:=2*a-15; c5:=k!6;
    eqs:=[2*N-3*U-c5*L^2,
          N^2+2*R-3*U^2-3*V^2-c4*L^2,
          2*V^3+2*N*R-U^3-6*U*V^2-c3*L^2,
          R^2+2*N*V^3-3*U^2*V^2-3*V^4-c2*L^2,
          2*R*V^3-3*U*V^4-c1*L^2];
    derivs:=[[Derivative(eqs[i],j):j in [1..5]]:i in [1..5]];
    uP:=(x-1)^2; sols:=[];
    raw:=0; gcd_bad:=0; pairing_den_bad:=0;
    for LL in k do
        if LL eq 0 then continue; end if;
        M:=LL^2;
        for UU in k do for VV in k do
            if VV eq 0 or UU^2-4*VV^2 eq 0 then continue; end if;
            NN:=(3*UU+c5*M)/2;
            RR:=(3*UU^2+3*VV^2+c4*M-NN^2)/2;
            pt:=[LL,UU,VV,NN,RR];
            if not &and[Evaluate(g,pt) eq 0:g in eqs] then continue; end if;
            raw+:=1;
            q:=x^2+UU*x+VV^2;
            if Degree(GCD(q,f)) gt 0 then gcd_bad+:=1; continue; end if;
            H:=x^3+NN*x^2+RR*x+VV^3;
            den:=LL^2*Resultant(uP,LL*h6-H);
            if den eq 0 then pairing_den_bad+:=1; continue; end if;
            pairing:=Resultant(q,H-LL*h6)/den;
            Jfix:=Matrix(k,5,5,
                &cat[[Evaluate(derivs[i][j],pt):j in [1..5]]:
                      i in [1..5]]);
            Append(~sols,<Z!LL,Z!UU,Z!VV,Z!NN,Z!RR,
                          Rank(Jfix),Z!pairing>);
        end for; end for;
    end for;
    return f,sols,<raw,gcd_bad,pairing_den_bad>;
end function;

fcommon,common_sols,common_stats:=CollisionData(k!8); // a=-2
fdc,dc_sols,dc_stats:=CollisionData(k!4);             // a=13, DC=0

assert Ebars["ZERO_2"] eq FK!(4*s^2);
assert Ebars["ZERO_MINUS3"] eq FK!(15*s^2);
nonempty_E:=[e0:e0 in [1..16]|#e3_by_E[e0] gt 0];
nonsquare_E:=[e0:e0 in [1..16]|not IsSquare(k!e0)];
assert nonempty_E eq nonsquare_E;
assert &and[&and[q[6] eq 5:q in e3_by_E[e0]]:e0 in nonempty_E];
assert common_stats eq <0,0,0> and #common_sols eq 0;
assert dc_stats eq <2,0,0> and #dc_sols eq 2;
assert &and[q[6] eq 5 and q[7] eq 1:q in dc_sols];

print "CONTACT6_M612_P8_P17_BAD_DISKS";
print "CLASSIFICATION_key_e_a_reason",classification;
print "ETALENESS_DERIVATIVES_key_de",derivatives;
print "ENDPOINT_EBARS",Ebars["ZERO_2"],Ebars["ZERO_MINUS3"];
print "E3_NONEMPTY_E_COUNTS",
      [<e0,#e3_by_E[e0],Sort([sol[6]:sol in e3_by_E[e0]])>:
       e0 in nonempty_E];
print "ZERO_2_ROWS_s_E_count",
      [<Z!s0,Z!Evaluate(Numerator(Ebars["ZERO_2"]),s0)/
                    Evaluate(Denominator(Ebars["ZERO_2"]),s0),
          #e3_by_E[Z!(Evaluate(Numerator(Ebars["ZERO_2"]),s0)/
                    Evaluate(Denominator(Ebars["ZERO_2"]),s0))]>:
       s0 in k|s0 ne 0];
print "ZERO_MINUS3_ROWS_s_E_count",
      [<Z!s0,Z!Evaluate(Numerator(Ebars["ZERO_MINUS3"]),s0)/
                    Evaluate(Denominator(Ebars["ZERO_MINUS3"]),s0),
          #e3_by_E[Z!(Evaluate(Numerator(Ebars["ZERO_MINUS3"]),s0)/
                    Evaluate(Denominator(Ebars["ZERO_MINUS3"]),s0))]>:
       s0 in k|s0 ne 0];
print "COMMON_ROOT_FIBER_FACTOR_DEGREES",
      [<Degree(fe[1]),fe[2]>:fe in Factorization(fcommon)],
      "RAW_GCD_BAD_PAIRING_DEN_BAD",common_stats,
      "SIGNED_OPEN_SOLUTIONS",#common_sols,
      "RANK_COUNTS",Sort([<r,#[q:q in common_sols|q[6] eq r]>:
                           r in Seqset([q[6]:q in common_sols])]),
      "PAIRING_COUNTS",Sort([<r,#[q:q in common_sols|q[7] eq r]>:
                              r in Seqset([q[7]:q in common_sols])]);
print "COMMON_ROOT_SAMPLES",common_sols[1..Min(#common_sols,20)];
print "DC_FIBER_FACTOR_DEGREES",
      [<Degree(fe[1]),fe[2]>:fe in Factorization(fdc)],
      "RAW_GCD_BAD_PAIRING_DEN_BAD",dc_stats,
      "SIGNED_OPEN_SOLUTIONS",#dc_sols,
      "RANK_COUNTS",Sort([<r,#[q:q in dc_sols|q[6] eq r]>:
                           r in Seqset([q[6]:q in dc_sols])]),
      "PAIRING_COUNTS",Sort([<r,#[q:q in dc_sols|q[7] eq r]>:
                              r in Seqset([q[7]:q in dc_sols])]);
print "DC_SAMPLES",dc_sols[1..Min(#dc_sols,20)];

print "CAVEAT",
      "rank-5 DC points certify punctured-disk branches; weighted common-root poles and non-E3 endpoint signatures are unresolved";
print "CONTACT6_M612_P8_P17_BAD_DISKS_DONE";
quit;
