//////////////////////////////////////////////////////////////////////
// First weighted layer in the four forced p=7 parameter disks of P8.
//
// The endpoint parameter e has a double zero in every disk.  For the
// primitive layer tau=tau0+7*s, write e=7^2*E.  The only contact-open
// endpoint signature of this valuation is E3:
//   (v(e),v(L),v(U),v(v),v(N),v(R))=(2,-1,-2,-2,-2,-4).
// This script derives the leading unit E(s), enumerates the exact E3
// initial system, and checks fixed-parameter Jacobian ranks.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q:=Rationals(); Z:=Integers(); k:=GF(7);
QPs<sQ>:=PolynomialRing(Q); FQ:=FieldOfFractions(QPs);
kPs<s>:=PolynomialRing(k); FK:=FieldOfFractions(kPs);

function ReducePolynomial(f)
    ans:=kPs!0;
    for i in [0..Degree(f)] do
        c:=Coefficient(f,i);
        assert Denominator(c) mod 7 ne 0;
        ans+:=(k!Numerator(c)/k!Denominator(c))*s^i;
    end for;
    return ans;
end function;

function ReduceRationalFunction(q)
    return FK!ReducePolynomial(Numerator(q))/FK!ReducePolynomial(Denominator(q));
end function;

centers:=[<"ZERO_2",Q!2>,<"ZERO_MINUS3",Q!-3>,
         <"POLE_1",Q!1>,<"POLE_MINUS1",Q!-1>];
Ebars:=AssociativeArray();
for rec in centers do
    label:=rec[1]; tau0:=rec[2];
    tau:=QPs!(tau0+7*sQ);
    AA:=tau^2+tau-6; BB:=tau^2+6;
    enum:=-400*AA^2*BB^2;
    eden:=768*AA^4-1200*AA^2*BB^2+1250*BB^4;
    enum49:=QPs!(enum/49);
    assert 49*enum49 eq enum;
    Ebar:=FK!ReducePolynomial(enum49)/FK!ReducePolynomial(eden);
    Ebars[label]:=Ebar;
end for;

// E3 initial equations after e=7^2 E, L=l/7,
// U=u/7^2, v=w/7^2, N=n/7^2, R=r/7^4.
S<E,L,U,V,N,R>:=PolynomialRing(k,6,"grevlex");
eqs:=[
    2*N-3*U-6*L^2,
    E*(N^2+2*R-3*U^2-3*V^2)-2*L^2,
    2*V^3+2*N*R-U^3-6*U*V^2,
    R^2+2*N*V^3-3*U^2*V^2-3*V^4,
    2*R*V^3-3*U*V^4
];

// Derive these forms from the exact five coefficient equations, rather
// than trusting a hand-entered initial system.  Multiplication by E^2 is
// harmless on the E != 0 chart and clears all E-denominators.
TQ<Eq,Lq,Uq,Vq,Nq,Rq>:=PolynomialRing(Q,6,"grevlex");
KTQ:=FieldOfFractions(TQ);
eex:=49*Eq; lex:=Lq/7; uex:=Uq/49; vex:=Vq/49;
nex:=Nq/49; rex:=Rq/2401;
h5:=2*nex-3*uex-6*lex^2;
h4:=nex^2+2*rex-3*uex^2-3*vex^2-(2/eex-15)*lex^2;
h3:=2*vex^3+2*nex*rex-uex^3-6*uex*vex^2-22*lex^2;
h2:=rex^2+2*nex*vex^3-3*uex^2*vex^2-3*vex^4
    -(1/eex^2-15)*lex^2;
h1:=2*rex*vex^3-3*uex*vex^4-(2/eex+6)*lex^2;
cleared:=[];
for rec in [<h5,2>,<h4,4>,<h3,6>,<h2,8>,<h1,10>] do
    g:=KTQ!(Eq^2*7^rec[2]*rec[1]);
    assert Denominator(g) eq 1;
    Append(~cleared,TQ!Numerator(g));
end for;

function ReduceMultivariate(g)
    ans:=S!0;
    for mon in Monomials(g) do
        c:=MonomialCoefficient(g,mon); ex:=Exponents(mon);
        assert Denominator(c) mod 7 ne 0;
        smon:=S!1;
        for i in [1..6] do smon*:=S.i^ex[i]; end for;
        ans+:=(k!Numerator(c)/k!Denominator(c))*smon;
    end for;
    return ans;
end function;

reduced_exact:=[ReduceMultivariate(g):g in cleared];
assert reduced_exact eq
       [E^2*eqs[1],E*eqs[2],E^2*eqs[3],E^2*eqs[4],E^2*eqs[5]];

// Weil pairing with the marked order-3 line.  The formula is validated
// independently in contact6_m612_p8_weil_pairing.m.  Compute its E3
// weighted initial value, so the degree-12 orthogonal supports can be
// separated from the degree-27 nonorthogonal supports.
PX<x>:=PolynomialRing(KTQ);
h6ex:=1+x/eex+x^3;
qex:=x^2+uex*x+vex^2;
Hex:=x^3+nex*x^2+rex*x+vex^3;
uPex:=(x-1)^2;
pairing:=Resultant(qex,Hex-lex*h6ex)/
         (lex^2*Resultant(uPex,lex*h6ex-Hex));

function CoeffVal7(c)
    return Valuation(Numerator(c),7)-Valuation(Denominator(c),7);
end function;

function NormalizePolynomial7(g)
    vals:=[CoeffVal7(MonomialCoefficient(g,mon)):mon in Monomials(g)];
    mn:=Min(vals);
    return TQ!(g/(Q!7)^mn),mn;
end function;

pnum,pnumval:=NormalizePolynomial7(TQ!Numerator(pairing));
pden,pdenval:=NormalizePolynomial7(TQ!Denominator(pairing));
assert pnumval eq pdenval;
pnumbar:=ReduceMultivariate(pnum);
pdenbar:=ReduceMultivariate(pden);
derivs_fixed:=[[Derivative(eqs[i],j):j in [2..6]]:i in [1..5]];

solutions_by_E:=AssociativeArray();
for e0 in k do
    sols:=[];
    if e0 ne 0 then
        for l0 in k do for u0 in k do for v0 in k do
        for n0 in k do for r0 in k do
            if l0 eq 0 or v0 eq 0 or u0^2-4*v0^2 eq 0 then continue; end if;
            pt:=[e0,l0,u0,v0,n0,r0];
            if not &and[Evaluate(g,pt) eq 0:g in eqs] then continue; end if;
            Jfix:=Matrix(k,5,5,
                &cat[[Evaluate(derivs_fixed[i][j-1],pt):j in [2..6]]
                      :i in [1..5]]);
            pn0:=Evaluate(pnumbar,pt); pd0:=Evaluate(pdenbar,pt);
            pair0:=pd0 eq 0 select -1 else Z!(pn0/pd0);
            Append(~sols,<Z!l0,Z!u0,Z!v0,Z!n0,Z!r0,
                          Rank(Jfix),pair0>);
        end for; end for; end for; end for; end for;
    end if;
    solutions_by_E[Z!e0]:=sols;
end for;

print "CONTACT6_M612_P8_P7_BAD_DISKS";
print "FIRST_LAYER_EBARS";
for rec in centers do print rec[1],Ebars[rec[1]]; end for;
assert Ebars["ZERO_2"] eq FK!(3*s^2);
assert Ebars["ZERO_MINUS3"] eq FK!(6*s^2);
assert Ebars["POLE_1"] eq FK!(2*(1+2*s)^2);
assert Ebars["POLE_MINUS1"] eq FK!(4*(1-2*s)^2);

print "E3_SOLUTIONS_BY_E";
for e0 in [1..6] do
    print e0,solutions_by_E[e0];
end for;
assert [e0:e0 in [1..6]|#solutions_by_E[e0] gt 0] eq [1,2,4];
assert &and[&and[sol[6] eq 5:sol in solutions_by_E[e0]]:
            e0 in [1,2,4]];
print "E3_PAIRING_INITIAL_NUMERATOR",pnumbar;
print "E3_PAIRING_INITIAL_DENOMINATOR",pdenbar;
print "E3_ORTHOGONAL_COUNTS_BY_E",
      [<e0,#[sol:sol in solutions_by_E[e0]|sol[7] eq 1],
             [sol[7]:sol in solutions_by_E[e0]]>:e0 in [1,2,4]];

for rec in centers do
    label:=rec[1]; Ebar:=Ebars[label];
    rows:=[]; deep:=[];
    for s0 in k do
        den0:=Evaluate(Denominator(Ebar),s0);
        if den0 eq 0 then Append(~deep,Z!s0); continue; end if;
        e0:=Evaluate(Numerator(Ebar),s0)/den0;
        if e0 eq 0 then Append(~deep,Z!s0); continue; end if;
        Append(~rows,<Z!s0,Z!e0,#solutions_by_E[Z!e0],
                       [<sol[6],sol[7]>:sol in solutions_by_E[Z!e0]]>);
    end for;
    print "DISK",label,"FIRST_LAYER_ROWS_s_E_count_fixedrank",rows,
          "DEEP_SUBDISKS",deep;
end for;

print "CONCLUSION_ZERO_DISKS",
      "no orthogonal E3 point for any primitive first-layer s; the same nonsquare unit excludes every scaled E3 depth";
print "CONCLUSION_POLE_DISKS",
      "two signed orthogonal fixed-parameter full-rank E3 lifts for every primitive first-layer s";
print "CAVEAT",
      "zero-disk non-E3 weighted signatures beyond the existing bounded endpoint fan remain unresolved";
print "CONTACT6_M612_P8_P7_BAD_DISKS_DONE";
quit;
