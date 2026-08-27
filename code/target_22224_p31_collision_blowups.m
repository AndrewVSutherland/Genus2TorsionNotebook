//////////////////////////////////////////////////////////////////////
// Exact checks supporting target_22224_p31_collision_blowups.py.
//
// Run from the repository root:
//   magma -b code/target_22224_p31_collision_blowups.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q:=Rationals();

if not assigned log_file then
    log_file:="results/target_22224_p31_collision_blowups_exact.log";
end if;
SetLogFile(log_file : Overwrite:=true);

print "TARGET_22224_P31_COLLISION_BLOWUPS_EXACT_START";

//////////////////////////////////////////////////////////////////////
// The full 2-descent coordinates of D4.
//////////////////////////////////////////////////////////////////////

R<a,b,c,d>:=PolynomialRing(Q,4);
RX<X>:=PolynomialRing(R);
g4:=(X-a*b)*(X-c*d)-X*(a+b)*(c+d);
vals:=[Evaluate(g4,z):z in [0,-a^2,-b^2,-c^2,-d^2]];
print "G4_EVALUATIONS_FACTORED";
for z in vals do print Factorization(z); end for;
assert vals[1] eq a*b*c*d;
assert vals[2] eq a*(a+b)*(a+c)*(a+d);
assert vals[3] eq b*(a+b)*(b+c)*(b+d);
assert vals[4] eq c*(a+c)*(b+c)*(c+d);
assert vals[5] eq d*(a+d)*(b+d)*(c+d);

//////////////////////////////////////////////////////////////////////
// The known record is a negative exact Kummer control.
//////////////////////////////////////////////////////////////////////

A:=Q!13/187; B:=Q!77/13; C:=Q!17/7;
rho:=Q!104880/17017; sigma:=Q!243120/17017; D:=sigma/rho;
assert A*B*C eq 1;
W:=[A,B,C,D];
K:=[];
for i in [1..4] do
    ki:=W[i];
    for j in [1..4] do if i ne j then ki*:=W[i]+W[j]; end if; end for;
    Append(~K,ki);
end for;
print "RECORD_NORMALIZED_W",W;
print "RECORD_KUMMER_PRODUCTS",K;
print "RECORD_KUMMER_SQUARES",[IsSquare(k):k in K],"D_square",IsSquare(D);
assert not &and[IsSquare(k):k in K];

//////////////////////////////////////////////////////////////////////
// Exact internal collision: x=y.
//
// If r^2=2x+1 and k^2=x(x+2), Y=2k satisfies
// Y^2=r^4+2r^2-3.  Its elliptic model has rank zero and torsion Z/4;
// all rational quartic points have r=+-1 (hence x=0) or are at infinity.
//////////////////////////////////////////////////////////////////////

Pr<rr>:=PolynomialRing(Q);
Cint:=HyperellipticCurve(rr^4+2*rr^2-3);
Pint:=Cint![1,0];
Eint,phiint:=EllipticCurve(Cint,Pint);
Eintmin,mint:=MinimalModel(Eint);
rlo,rhi:=RankBounds(Eintmin);
Tint,mapTint:=TorsionSubgroup(Eintmin);
print "INTERNAL_QUARTIC",Cint;
print "INTERNAL_ELLIPTIC_MODEL",Eint,"minimal",Eintmin;
print "INTERNAL_RANK_BOUNDS",rlo,rhi,"torsion",Invariants(Tint);
assert rlo eq 0 and rhi eq 0 and Invariants(Tint) eq [4];
invphiint:=Inverse(phiint); invmint:=Inverse(mint);
intpts:=[invphiint(invmint(mapTint(t))):t in Tint];
print "INTERNAL_ALL_RATIONAL_POINTS",intpts;
print "INTERNAL_AFFINE_R_VALUES",[p[1]/p[3]:p in intpts|p[3] ne 0];
assert {p[1]/p[3]:p in intpts|p[3] ne 0} eq {Q!-1,Q!1};

//////////////////////////////////////////////////////////////////////
// Exact same-sign external collision A=D=1.
//
// The Kummer square conditions are B=2*t^2 and Y^2=4*t^4+1.
// This quartic also has rank zero.  Its affine rational points have t=0,
// so B=0; the remaining points are at infinity.  Thus the exact boundary
// has no nonzero rational target point, although it has live Q_31 disks.
//////////////////////////////////////////////////////////////////////

Pt<tt>:=PolynomialRing(Q);
Cext:=HyperellipticCurve(4*tt^4+1);
Pext:=Cext![0,1];
Eext,phiext:=EllipticCurve(Cext,Pext);
Eextmin,mext:=MinimalModel(Eext);
elo,ehi:=RankBounds(Eextmin);
Text,mapText:=TorsionSubgroup(Eextmin);
print "EXTERNAL_QUARTIC",Cext;
print "EXTERNAL_ELLIPTIC_MODEL",Eext,"minimal",Eextmin;
print "EXTERNAL_RANK_BOUNDS",elo,ehi,"torsion",Invariants(Text);
assert elo eq 0 and ehi eq 0 and Invariants(Text) eq [2,2];
invphiext:=Inverse(phiext); invmext:=Inverse(mext);
extpts:=[invphiext(invmext(mapText(t))):t in Text];
print "EXTERNAL_ALL_RATIONAL_POINTS",extpts;
print "EXTERNAL_AFFINE_T_VALUES",[p[1]/p[3]:p in extpts|p[3] ne 0];
assert {p[1]/p[3]:p in extpts|p[3] ne 0} eq {Q!0};

print "TARGET_22224_P31_COLLISION_BLOWUPS_EXACT_DONE";
UnsetLogFile();
quit;
