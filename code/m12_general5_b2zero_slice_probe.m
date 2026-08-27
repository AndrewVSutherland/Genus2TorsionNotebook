//////////////////////////////////////////////////////////////////////
// Memory-bounded rational hyperplane probes of the generic linear-B
// M(12)+5 sign quotient.
//
// This deliberately works in the five-variable chart (b,w,c,d,e), whose
// equations are much smaller than the eliminated four-variable equations.
// It intersects the expected quotient curve with one rational hyperplane,
// computes a zero-dimensional Groebner basis when possible, and reports
// rational points only after substitution into all four exact equations.
//
// Examples:
//   magma -b var:="c" num:=1 den:=1 MemGB:=2 \
//       code/m12_general5_b2zero_slice_probe.m
//   magma -b var:="d" num:=-1 den:=1 MemGB:=2 \
//       code/m12_general5_b2zero_slice_probe.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned var then var := "c"; end if;
if var notin {"b","w","c","d","e"} then
    error "var must be b, w, c, d, or e";
end if;
if not assigned num then num := 0;
elif Type(num) eq MonStgElt then num := StringToInteger(num); end if;
if not assigned den then den := 1;
elif Type(den) eq MonStgElt then den := StringToInteger(den); end if;
require den ne 0: "den must be nonzero";
if not assigned MemGB then MemGB := 2;
elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();
R<b,w,c,d,e> := PolynomialRing(Q,5,"grevlex");
P<Z> := PolynomialRing(R);
x := d*Z-c;
L := b+(2*b-1)*x;
H := x+w*(1+b*x);
G := P!(L*(L*H^2+4*b*(1+x)^2*(w*L-x^2)));
g := [R!Coefficient(G,i) : i in [0..5]];
g0:=g[1]; g1:=g[2]; g2:=g[3]; g3:=g[4]; g4:=g[5]; g5:=g[6];

ell3 := (5*e/2)*g0-g1;
ell4 := (5/2+15*e^2/8)*g0-g2;
ell5 := (5/2+15*e^2/8)*g0-g3;
ell6 := (5*e/2)*g0-g4;
ell7 := g0-g5;
delta := (2-e)^3;
C1 := ell3-ell7;
C2 := (32*e^2-33*e+58)*ell3-20*ell5;
C3 := (19*e-6)*ell3-8*ell6;
C4 := 4*(19*e-6)*ell3^2-32*ell4*ell3+5*delta*g0^2;
eqs := [C1,C2,C3,C4];

idx := Index(["b","w","c","d","e"],var);
S<t1,t2,t3,t4> := PolynomialRing(Q,4,"grevlex");
remaining := [i:i in [1..5]|i ne idx];
vals := [S!0:i in [1..5]];
vals[idx] := Q!num/den;
for i in [1..#remaining] do vals[remaining[i]] := S.i; end for;
phi := hom<R -> S | vals>;
eqs_slice := [phi(f):f in eqs];
I := ideal<S|eqs_slice>;
print "SLICE_BEGIN",var,Rationals()!num/den;
print "REMAINING_VARIABLES",[["b","w","c","d","e"][i]:i in remaining];
print "EQUATION_SHAPES",[<TotalDegree(f),#Terms(f)> : f in eqs_slice];
t0 := Cputime();
B := Basis(I);
print "RAW_BASIS_END","cpu",Cputime(t0),"len",#B,
      "degrees",[TotalDegree(f):f in B],"terms",[#Terms(f):f in B];
if &or[f ne 0 and TotalDegree(f) eq 0:f in B] then
    print "SLICE_EMPTY_RAW";
    quit;
end if;

// Only cheap chart factors are saturated here.  Large discriminants and
// resultants are checked on any recovered point instead of being inserted
// into the initial ideal.
cheap_names := ["b","w","b-1","2b-1","d","e-2","e+2","ell3"];
cheap := [phi(f):f in [b,w,b-1,2*b-1,d,e-2,e+2,ell3]];
J := I;
for i in [1..#cheap] do
    if cheap[i] eq 0 then
        print "CHEAP_FACTOR_ZERO_ON_SLICE",cheap_names[i];
        print "SLICE_EMPTY_CHEAP";
        quit;
    end if;
    nf := NormalForm(cheap[i],Basis(J));
    if nf eq 0 then
        print "CHEAP_FACTOR_CONTAINS_SLICE",cheap_names[i];
    end if;
    t0 := Cputime();
    J := Saturation(J,ideal<S|nf>);
    BJ := Basis(J);
    print "SAT_END",cheap_names[i],"cpu",Cputime(t0),"len",#BJ,
          "unit",&or[f ne 0 and TotalDegree(f) eq 0:f in BJ];
    if &or[f ne 0 and TotalDegree(f) eq 0:f in BJ] then
        print "SLICE_EMPTY_CHEAP";
        quit;
    end if;
end for;

BJ := Basis(J);
try
    dim,ds := Dimension(J);
    print "CHEAP_DIMENSION",dim,"COMPONENT_DEGREES",ds;
catch err
    print "CHEAP_DIMENSION_FAILED",err`Object;
end try;

// ChangeOrder uses the already-computed zero-dimensional grevlex quotient
// (FGLM where applicable); recreating a lex ideal from its generators makes
// Magma start a much more expensive lex Buchberger computation.
t0 := Cputime();
JL := ChangeOrder(J,"lex");
BL := GroebnerBasis(JL);
print "LEX_BASIS_END","cpu",Cputime(t0),"len",#BL,
      "max_degree",Maximum([TotalDegree(f):f in BL]),
      "max_terms",Maximum([#Terms(f):f in BL]);
for f in BL do
    if #Terms(f) le 30 then print "LEX_POLY",f; end if;
end for;
RL := Generic(JL);
for f in BL do
    support_vars := [j:j in [1..4]|Degree(f,RL.j) gt 0];
    if #support_vars le 1 and TotalDegree(f) gt 1 then
        fac := Factorization(f);
        print "LEX_UNIVARIATE","variable_indices",support_vars,
              "variable_name",["b","w","c","d","e"][remaining[support_vars[1]]],
              "degree",TotalDegree(f),"terms",#Terms(f);
        print "LEX_UNIVARIATE_FACTOR_SHAPES",
              [<TotalDegree(fe[1]),fe[2]>:fe in fac];
        for fe in fac do
            if TotalDegree(fe[1]) eq 1 then
                print "LEX_RATIONAL_FACTOR",fe;
            end if;
        end for;
    end if;
end for;
print "SLICE_DONE";
