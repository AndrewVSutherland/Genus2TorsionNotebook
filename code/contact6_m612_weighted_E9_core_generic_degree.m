//////////////////////////////////////////////////////////////////////
// Bounded generic-degree diagnostic for the b=0 cubic-contact core.
//
// Work over Q(e), quotient the harmless sign L -> -L by M=L^2, eliminate
// the two linear variables N,R, and remove the obvious M=0 component.
// The remaining zero-dimensional generic fiber has length 40.  This is a
// warning about the ambient core cover, not an irreducibility assertion for
// the particular E9 component and not a full contact-open saturation.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q:=Rationals(); K<e>:=FunctionField(Q);
R<M,U,v>:=PolynomialRing(K,3,"grevlex");

N:=(3*U+6*M)/2;
Rc:=(3*U^2+3*v^2+(2/e-15)*M-N^2)/2;
F3:=2*v^3+2*N*Rc-U^3-6*U*v^2-22*M;
F2:=Rc^2+2*N*v^3-3*U^2*v^2-3*v^4-(1/e^2-15)*M;
F1:=2*Rc*v^3-3*U*v^4-(2/e+6)*M;
eqs:=[R!Numerator(f):f in [F3,F2,F1]];

// The raw positive-dimensional contaminant is visible without elimination.
assert &and[Evaluate(f,[K!0,2*v,v]) eq 0:f in eqs];
Iraw:=ideal<R|eqs>;
dimraw:=Dimension(Iraw);
assert dimraw eq 1;

t0:=Cputime();
I:=Saturation(Iraw,ideal<R|M>);
sat_seconds:=Cputime(t0);
dim:=Dimension(I);
assert dim eq 0;

LMs:=[LeadingMonomial(g):g in Basis(I)];
exps:=[Exponents(mm):mm in LMs];
bounds:=[];
for j in [1..3] do
    vals:=[a[j]:a in exps|a[j] gt 0 and
          &and[a[h] eq 0:h in [1..3]|h ne j]];
    assert #vals gt 0;
    Append(~bounds,Min(vals));
end for;
degree:=0;
for a in [0..bounds[1]-1] do for b in [0..bounds[2]-1] do
for c in [0..bounds[3]-1] do
    q:=[a,b,c];
    if &or[&and[q[j] ge lm[j]:j in [1..3]]:lm in exps] then continue; end if;
    degree+:=1;
end for; end for; end for;
assert degree eq 40;

print "CONTACT6_M612_WEIGHTED_E9_CORE_GENERIC_DEGREE";
print "VARIABLES",["M=L^2","U","v"],"BASE","Q(e)";
print "ELIMINATED_N",N;
print "ELIMINATED_R",Rc;
print "EQUATION_SHAPES",[<TotalDegree(f),#Terms(f)>:f in eqs];
print "RAW_DIMENSION",dimraw,
      "OBVIOUS_COMPONENT","M=0,U=2v";
print "M_SATURATION_SECONDS",sat_seconds,"SATURATED_DIMENSION",dim;
print "LEADING_MONOMIALS",LMs;
print "PURE_POWER_BOUNDS",bounds,"GENERIC_FIBER_LENGTH",degree;
print "WARNING",
      "length 40 is before full contact-open decomposition and does not rule out a degree-one E9 component";
print "CONTACT6_M612_WEIGHTED_E9_CORE_GENERIC_DEGREE_DONE";
quit;
