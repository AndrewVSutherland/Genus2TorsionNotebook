//////////////////////////////////////////////////////////////////////
// Plane elimination and rationality test for the genus-zero rational
// 3-torsion cover on the minus elliptic quotient of (1,r,s,rs).
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals(); Z := Integers();
B<r,y> := PolynomialRing(Q,2);
K := FieldOfFractions(B); PU<u> := PolynomialRing(K);

function PrimitiveNumerator(g)
    h := B!Numerator(g);
    cc := [c:c in Coefficients(h)|c ne 0];
    den := #cc eq 0 select 1 else LCM([Denominator(c):c in cc]);
    h *:= den;
    nums := [Z!c:c in Coefficients(h)|c ne 0];
    content := #nums eq 0 select 1 else GCD([Abs(n):n in nums]);
    return content gt 1 select h/content else h;
end function;

for family in [5,6] do
    AA := family eq 5 select 16 else 144;
    BB := family eq 5 select 9 else 25;
    ss := K!(BB-AA*r)/(AA-BB*r); kk := K!r*ss;
    A := -2*kk; C1 := 1+kk^2; C2 := r^2+ss^2;
    a2 := A+C1+C2; a4 := A*C1+A*C2+C1*C2; a6 := A*C1*C2;
    psi := 3*u^4+4*a2*u^3+6*a4*u^2+12*a6*u+(4*a2*a6-a4^2);
    cubic := (u+A)*(u+C1)*(u+C2);
    elim := PrimitiveNumerator(Resultant(psi,u^0*y^2-cubic));
    fac := Factorization(elim);
    print "PRODUCT_PLANE_ELIM","family",family,
          "degree_r",Degree(elim,1),"degree_y",Degree(elim,2),
          "total_degree",TotalDegree(elim),"terms",#Terms(elim),
          "factors",[<Degree(z[1],1),Degree(z[1],2),z[2]>:z in fac];
    assert #fac eq 1;
    F := fac[1][1];
    A2 := AffineSpace(Q,2);
    Caff := Curve(A2,F);
    Cproj := ProjectiveClosure(Caff);
    FC := FunctionField(Cproj);
    AF,toAF := AlgorithmicFunctionField(FC);
    print "PRODUCT_PLANE_FUNCTION_FIELD","family",family,"genus",Genus(AF);
    pls := Places(AF,1);
    print "PRODUCT_PLANE_DEGREE1_PLACES","family",family,#pls;
    if #pls gt 0 then
        z,images := Parametrization(AF,Divisor(pls[1]));
        print "PRODUCT_PLANE_PARAMETRIZATION","family",family,
              "parameter",z,"images",images;
    end if;
end for;

quit;
