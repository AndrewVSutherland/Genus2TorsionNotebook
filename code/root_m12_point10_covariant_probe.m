//////////////////////////////////////////////////////////////////////
// Low-memory point-order-10 probe on the compact M(12) chart.
//
// Move the visible Weierstrass root L=b+(2b-1)x to z=0 and multiply
// the quintic by the square (2b-1)^4.  For P=(r,c), put E=P-infinity
// and T=(0,0)-infinity.  The exact relation 5E=T is equivalent to
//
//   F(z)-H(z)^2 = kappa*z*(z-r)^5,  H(0)=0,
//
// for a cubic H.  Matching the Taylor expansion of H^2 to F through
// order four gives two covariants in (b,w,r); F(r) must additionally
// be a nonzero rational square.
//
// Modes:
//   summary : factor/report the two characteristic-zero covariants;
//   local   : exhaustive finite-field census and exact group checks.
//
// Examples:
//   magma -b mode:=summary MemGB:=2 code/root_m12_point10_covariant_probe.m
//   magma -b mode:=local p:=7 MemGB:=2 code/root_m12_point10_covariant_probe.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned mode then mode := "summary"; end if;
if mode notin {"summary","local"} then
    error "mode must be summary or local";
end if;
if not assigned p then p := 7;
elif Type(p) eq MonStgElt then p := StringToInteger(p); end if;
if not assigned MemGB then MemGB := 2;
elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
if not assigned sample_limit then sample_limit := 12;
elif Type(sample_limit) eq MonStgElt then
    sample_limit := StringToInteger(sample_limit);
end if;
SetMemoryLimit(MemGB*10^9);

function TransformedQuintic(P,b,w)
    z := P.1;
    d := 2*b-1;
    // x=(z-b)/d, L=z, and d^4 is a square scaling of the y-model.
    Hn := z-b+w*(d+b*(z-b));
    return P!(z*(d^2*z*Hn^2
           +4*b*(z+b-1)^2*(w*z*d^2-(z-b)^2)));
end function;

function Covariants(P,F,r)
    A0 := Evaluate(F,r);
    A1 := Evaluate(Derivative(F),r);
    A2 := Evaluate(Derivative(F,2),r)/2;
    A3 := Evaluate(Derivative(F,3),r)/6;
    A4 := Evaluate(Derivative(F,4),r)/24;
    Q2 := 4*A0*A2-A1^2;
    E3 := 8*A0^2*A3-A1*Q2;
    // Existence of a cubic Taylor square root through order four.
    E4 := 64*A0^3*A4-Q2^2-4*A1*E3;
    // Its value at the Weierstrass root z=0 is zero.
    W0 := 16*A0^3-8*A1*A0^2*r+2*Q2*A0*r^2-E3*r^3;
    return E4,W0,A0,A1,A2,A3,A4,Q2,E3;
end function;

if mode eq "summary" then
    Q := Rationals();
    R<b,w,r> := PolynomialRing(Q,3,"grevlex");
    P<z> := PolynomialRing(R);
    F := TransformedQuintic(P,b,w);
    E4,W0,A0,A1,A2,A3,A4,Q2,E3 := Covariants(P,F,r);
    print "M12_POINT10_COVARIANT_SUMMARY";
    print "F_degree",Degree(F),"F_terms",#Terms(F),
          "F_lc",LeadingCoefficient(F),"F_const",Coefficient(F,0);
    for item in [<"E4",E4>,<"W0",W0>] do
        g := item[2];
        print item[1],"degree",TotalDegree(g),"degree_r",Degree(g,r),
              "terms",#Terms(g);
        fac := Factorization(g);
        print item[1],"factors",
              [<TotalDegree(fe[1]),Degree(fe[1],r),#Terms(fe[1]),fe[2]>:
               fe in fac];
        for fe in fac do
            if #Terms(fe[1]) le 80 then print item[1],"FACTOR",fe; end if;
        end for;
    end for;
    G := GCD(E4,W0);
    print "GCD","degree",TotalDegree(G),"terms",#Terms(G),G;
    assert Evaluate(F,0) eq 0;
    print "M12_POINT10_COVARIANT_SUMMARY_DONE";
    quit;
end if;

require IsPrime(p) and p notin {2,3,5}: "use a prime other than 2,3,5";
k := GF(p);
P<z> := PolynomialRing(k);
elts := [a:a in k];
raw := 0; open := 0; square_lifts := 0; exact10 := 0;
cyclic60_finite := 0; samples := 0;
for b in elts do for w in elts do
    d := 2*b-1;
    if b*w*(b-1)*d eq 0 then continue; end if;
    F := TransformedQuintic(P,b,w);
    if Degree(F) ne 5 or Discriminant(F) eq 0 then continue; end if;
    C := HyperellipticCurve(F); J := Jacobian(C);
    T := J![z,0];
    if Order(T) ne 2 then continue; end if;
    z12 := 1-b;
    eta12 := d^2*(1-b)*(w*(1-b)-1);
    if Evaluate(F,z12) ne eta12^2 then continue; end if;
    D12 := J![z-z12,eta12];
    if Order(D12) ne 12 then continue; end if;
    for r in elts do
        E4,W0,A0,A1,A2,A3,A4,Q2,E3 := Covariants(P,F,r);
        if E4 ne 0 or W0 ne 0 then continue; end if;
        raw +:= 1;
        if r eq 0 or A0 eq 0 then continue; end if;
        open +:= 1;
        if not IsSquare(A0) then continue; end if;
        square_lifts +:= 1;
        c := Sqrt(A0);
        E := J![z-r,c];
        if 5*E ne T or Order(E) ne 10 then
            print "EXACT_RELATION_FAILURE",p,b,w,r,c,
                  "ordE",Order(E),"fiveE",5*E,"T",T;
            continue;
        end if;
        exact10 +:= 1;
        D60 := D12+E;
        if Order(D60) eq 60 then cyclic60_finite +:= 1; end if;
        if samples lt sample_limit then
            samples +:= 1;
            print "LOCAL_SAMPLE","p",p,"b",b,"w",w,"r",r,"c",c,
                  "ordE",Order(E),"ordD12",Order(D12),
                  "ordSum",Order(D60),"Jorder",#J;
        end if;
    end for;
end for; end for;
print "LOCAL_COUNT","p",p,"raw",raw,"open",open,
      "square_lifts",square_lifts,"exact10",exact10,
      "cyclic60_finite",cyclic60_finite,"samples",samples;
print "M12_POINT10_COVARIANT_LOCAL_DONE";
quit;
