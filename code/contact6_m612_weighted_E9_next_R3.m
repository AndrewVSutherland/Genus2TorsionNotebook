//////////////////////////////////////////////////////////////////////
// R3 compatibility for the continued p=5 E9 core residues.
//
// The core scaling is e=25E.  The mod-25 core residues that continue
// to the regular blowup have
//
//   E=1,6,11,16,21 or E=4,9,14,19,24 (mod 25).
//
// For each resulting e modulo 625, uniquely Hensel-lift both smooth
// endpoint R3 solutions (mu,nu)=(1,0),(4,0) from mod 5 through mod 625.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Z:=Integers();
Q:=Rationals();
p:=5;
k:=GF(p);

A<e,mu,nu>:=PolynomialRing(Q,3,"grevlex");
PX<X>:=PolynomialRing(A);

D:=3+5*e;
A1:=(-2*e-3*e^2)*X^2+(6*e+10*e^2)*X+(1-3*e^2);
A2:=2*e*X^2-(1+3*e);
R3:=2-3*X^2;
S:=R3*(mu*X+nu)^2-D*A1*A2;
s4:=Coefficient(S,4);
s3:=Coefficient(S,3);
s2:=Coefficient(S,2);
s1:=Coefficient(S,1);
s0:=Coefficient(S,0);
H1:=8*s4^2*s1-s3*(4*s4*s2-s3^2);
H0:=64*s4^3*s0-(4*s4*s2-s3^2)^2;
Hs:=[H1,H0];

function EvalZ(g,pt)
    z:=Evaluate(g,[Q!a:a in pt]);
    assert Denominator(z) eq 1;
    return Z!z;
end function;

function LiftR3(evalue,pt,modulus)
    nextmod:=p*modulus;
    lifts:=[];
    for dm in [0..p-1] do for dn in [0..p-1] do
        q:=[(pt[1]+modulus*dm) mod nextmod,
           (pt[2]+modulus*dn) mod nextmod];
        full:=[evalue,q[1],q[2]];
        if &and[(EvalZ(h,full) mod nextmod) eq 0:h in Hs] then
            Append(~lifts,q);
        end if;
    end for; end for;
    assert #lifts eq 1;
    return lifts[1];
end function;

Eresidues:=[1,6,11,16,21,4,9,14,19,24];
print "CONTACT6_M612_WEIGHTED_E9_NEXT_R3";
print "CORE_E_RESIDUES_MOD25",Eresidues;
for E0 in Eresidues do
    e0:=(25*E0) mod 625;
    branches:=[];
    for root in [[1,0],[4,0]] do
        q:=root;
        q:=LiftR3(e0,q,5);
        assert q in [[4,0],[21,0]];
        q:=LiftR3(e0,q,25);
        q:=LiftR3(e0,q,125);
        assert &and[(EvalZ(h,[e0,q[1],q[2]]) mod 625) eq 0:h in Hs];
        assert (EvalZ(s4,[e0,q[1],q[2]]) mod 5) ne 0;
        Append(~branches,q);
    end for;
    assert branches[1] ne branches[2];
    print "E_MOD25",E0,"e_MOD625",e0,
          "R3_BRANCHES_mu_nu_MOD625",branches;
end for;

print "R3_VERDICT",
      "both smooth R3 branches extend through mod 625 for every continued core E residue";
print "CONTACT6_M612_WEIGHTED_E9_NEXT_R3_DONE";
quit;
