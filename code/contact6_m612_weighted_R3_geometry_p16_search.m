//////////////////////////////////////////////////////////////////////
// Rational-point reconnaissance on the P16 component through its
// elliptic quotient z=mu^2.
//
// The quotient is E: y^2=x^3-12x of rank one.  Pull small Mordell-Weil
// combinations back to (e,z) and retain those with z a rational square.
// Also solve the exact H1,H0 fibers above visible z=0 plane points.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MWBound then MWBound:=20; end if;
if Type(MWBound) eq MonStgElt then MWBound:=StringToInteger(MWBound); end if;

Q:=Rationals();
A<e,mu,nu>:=PolynomialRing(Q,3,"grevlex");
PX<X>:=PolynomialRing(A);
D:=3+5*e;
A1:=(-2*e-3*e^2)*X^2+(6*e+10*e^2)*X+(1-3*e^2);
A2:=2*e*X^2-(1+3*e);
R3:=2-3*X^2;
S:=R3*(mu*X+nu)^2-D*A1*A2;
s4:=Coefficient(S,4); s3:=Coefficient(S,3);
s2:=Coefficient(S,2); s1:=Coefficient(S,1); s0:=Coefficient(S,0);
H1:=8*s4^2*s1-s3*(4*s4*s2-s3^2);
H0:=64*s4^3*s0-(4*s4*s2-s3^2)^2;

B<ee,mm>:=PolynomialRing(Q,2,"grevlex");
toB:=hom<A -> B | ee,mm,B!0>;
fac:=Factorization(B!toB(Resultant(H1,H0,nu)));
P16:=[f[1]:f in fac|Degree(f[1],ee) eq 16 and Degree(f[1],mm) eq 8][1];
T<et,z>:=PolynomialRing(Q,2,"grevlex");
coeffToT:=hom<B -> T | et,T!0>;
Pz:=&+[coeffToT(B!Coefficient(P16,mm,2*j))*z^j:j in [0..4]];
Cz:=ProjectiveClosure(Curve(AffineSpace(T),Pz));
q0:=Cz![-2/5,4/75,1];
assert not IsSingular(q0);
Eraw,mp:=EllipticCurve(Cz,q0);
E,mmin:=MinimalModel(Eraw);
assert aInvariants(E) eq [0,0,0,-12,0];

print "CONTACT6_M612_WEIGHTED_R3_GEOMETRY_P16_SEARCH";
print "ELLIPTIC_QUOTIENT",E,"RANK_BOUNDS",RankBounds(E);
print "MINIMAL_MAP_DOMAIN",Domain(mmin),"CODOMAIN",Codomain(mmin);
print "CURVE_MAP_HAS_INVERSE",HasKnownInverse(mp);

// Exact fibers over rational z=0 points seen on the plane quotient.
Un<n>:=PolynomialRing(Q);
for ev in [Q!0,Q!-3/5,Q!-6/17] do
    h1:=Un!Evaluate(H1,[ev,Q!0,n]);
    h0:=Un!Evaluate(H0,[ev,Q!0,n]);
    gg:=GCD(h1,h0);
    print "ZERO_MU_FIBER","e",ev,"s4",Evaluate(s4,[ev,0,0]),
          "H1",h1,"H0",h0,"GCD",gg;
    if gg ne 0 then print " ZERO_MU_ROOTS",Roots(gg); end if;
end for;

if HasKnownInverse(mp) then
try
    rawToMin:=mmin;
    minToRaw:=Inverse(mmin);
    curveToRaw:=mp;
    rawToCurve:=Inverse(mp);
    print "INVERSE_MAPS_READY";
    gens:=Generators(E);
    print "MW_GENERATORS",gens;
    TG,tmap:=TorsionSubgroup(E);
    torspts:=[tmap(t):t in TG];
    hits:=[];
    checked:=0;
    for nmult in [-MWBound..MWBound] do
        for tp in torspts do
            P:=tp;
            if #gens gt 0 then P+:=nmult*gens[1]; end if;
            try
                Praw:=minToRaw(P);
                q:=rawToCurve(Praw);
                cc:=Coordinates(q);
                if cc[3] eq 0 then continue; end if;
                ev:=cc[1]/cc[3]; zv:=cc[2]/cc[3];
                checked+:=1;
                ok,mv:=IsSquare(zv);
                if ok then
                    Append(~hits,<nmult,tp,ev,zv,mv,q>);
                    print "SQUARE_Z_HIT",hits[#hits];
                end if;
            catch err
                continue;
            end try;
        end for;
    end for;
    print "MW_SEARCH_BOUND",MWBound,"PULLBACKS_CHECKED",checked,
          "SQUARE_Z_HIT_COUNT",#hits;
catch err
    print "MW_PULLBACK_SEARCH_FAILED",err`Object;
end try;
else
    print "MW_PULLBACK_SEARCH_SKIPPED","curve-to-elliptic map has no stored inverse";
end if;
print "CONTACT6_M612_WEIGHTED_R3_GEOMETRY_P16_SEARCH_DONE";
quit;
