//////////////////////////////////////////////////////////////////////
// Exact affine p=5,7 boundary congruences on the compact T_B cover.
//
// On aden!=0, a=anum/aden.  A rational [6,12] target reducing to finite
// (s,r,K) at p=5 or 7 must lie on the displayed curve-discriminant boundary,
// because the good-open local target mask is empty at both primes.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q:=Rationals(); Z:=Integers();
R<s,r,K>:=PolynomialRing(Q,3,"grevlex");

aden:=s*(r^2-K)-r;
anum:=-(s*(-3*r^2+4*r-3*K)-r*(-3+4*s^2*r-2*K));

// Rebuild the irreducible compact cover H.
RA<a0,s0,r0,K0>:=PolynomialRing(Q,4,"grevlex");
A3:=a0-3+4*s0^2*r0-2*K0;
H1:=s0*((a0-3)*r0^2+4*r0-K0*(a0+3))-r0*A3;
H2:=8*s0^2*(2*s0^2*r0^2+2*(a0-3)*r0+2
             -K0*(2*s0^2-6))-A3^2-32*s0^3*r0;
dA:=Coefficient(H1,a0,1);
nA:=-Evaluate(H1,[RA!0,s0,r0,K0]);
degA:=Degree(H2,a0);
Hred:=&+[Coefficient(H2,a0,i)*nA^i*dA^(degA-i):i in [0..degA]];
assert IsDivisibleBy(Hred,s0^2*K0);
rawH:=ExactQuotient(Hred,s0^2*K0);
toR:=hom<RA -> R | R!0,s,r,K>;
H:=R!toR(rawH);

function Primitive(g)
    if g eq 0 then return g; end if;
    den:=LCM([Denominator(c):c in Coefficients(g)]);
    vals:=[Z!(den*c):c in Coefficients(g)];
    cont:=GCD(vals); if cont eq 0 then cont:=1; end if;
    return Parent(g)!((Q!den/Q!Abs(cont))*g);
end function;
H:=Primitive(H);

KF:=FieldOfFractions(R);
PX<X>:=PolynomialRing(KF);
aa:=KF!anum/KF!aden;
bb:=2*(KF!s)^2-3;
q1:=(bb+3)*X^2+(aa-3)*X+2;
q2:=2*X^2+(bb-3)*X+(aa+3);
f:=X*q1*q2;
discNum:=Primitive(R!Numerator(Discriminant(f)));

print "CONTACT6_M612_TB_BOUNDARY_MODP";
print "H_SHAPE",<TotalDegree(H),Degree(H,K),#Terms(H)>;
print "DISCNUM_SHAPE",
      <TotalDegree(discNum),Degree(discNum,s),Degree(discNum,r),
       Degree(discNum,K),#Terms(discNum)>;
print "DISCNUM_Q_FACTORIZATION",
      [<TotalDegree(fe[1]),Degree(fe[1],s),Degree(fe[1],r),
        Degree(fe[1],K),#Terms(fe[1]),fe[2]>
       : fe in Factorization(discNum)];

for p in [5,7] do
    k:=GF(p);
    Rp<sp,rp,kp>:=PolynomialRing(k,3,"grevlex");
    mp:=hom<R -> Rp | sp,rp,kp>;
    Hp:=Rp!mp(H);
    Dp:=Rp!mp(discNum);
    facD:=Factorization(Dp);
    print "PRIME",p;
    print " DISC_FACTORIZATION_SHAPES",
          [<TotalDegree(fe[1]),Degree(fe[1],sp),Degree(fe[1],rp),
            Degree(fe[1],kp),#Terms(fe[1]),fe[2]>:fe in facD];
    for fe in facD do
        if #Terms(fe[1]) le 40 then print " DISC_FACTOR",fe; end if;
    end for;
    print " H_DISC_GCD",Factorization(GCD(Hp,Dp));

    boundaryTriples:=[];
    recoveryTriples:=[];
    squareBoundaryTriples:=[];
    goodAllowed:=0;
    squares:={z^2:z in k|z ne 0};
    for sv in k do for rv in k do for kv in k do
        vals:=[sv,rv,kv];
        if Evaluate(Hp,vals) ne 0 then continue; end if;
        dv:=Evaluate(mp(aden),vals);
        if dv eq 0 then
            Append(~recoveryTriples,<Z!sv,Z!rv,Z!kv>);
            continue;
        end if;
        if Evaluate(Dp,vals) eq 0 then
            triple:=<Z!sv,Z!rv,Z!kv>;
            Append(~boundaryTriples,triple);
            if kv in squares then Append(~squareBoundaryTriples,triple); end if;
            continue;
        end if;
        if kv notin squares then continue; end if;
        av:=Evaluate(mp(anum),vals)/dv;
        bv:=2*sv^2-3;
        Pk<xx>:=PolynomialRing(k);
        ff:=xx*((bv+3)*xx^2+(av-3)*xx+2)
              *(2*xx^2+(bv-3)*xx+(av+3));
        if Degree(ff) ne 5 or Discriminant(ff) eq 0 then
            continue;
        end if;
        AG,phi:=AbelianGroup(Jacobian(HyperellipticCurve(ff)));
        inv:=Invariants(AG);
        has612:=#[n:n in inv|(Z!n) mod 6 eq 0] ge 2
             and #[n:n in inv|(Z!n) mod 12 eq 0] ge 1;
        if has612 then goodAllowed+:=1; end if;
    end for; end for; end for;
    print " AFFINE_BOUNDARY_TRIPLES",boundaryTriples;
    print " SQUARE_K_BOUNDARY_TRIPLES",squareBoundaryTriples;
    print " RECOVERY_TRIPLES",recoveryTriples;
    print " GOOD_ALLOWED_612",goodAllowed;
end for;
print "CONTACT6_M612_TB_BOUNDARY_MODP_DONE";
quit;
