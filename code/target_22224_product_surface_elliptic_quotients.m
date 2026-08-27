//////////////////////////////////////////////////////////////////////
// Elliptic-quotient geometry of the product surface
//             (a,b,c,d)=(1,r,s,r*s).
//
// Put k=r*s.  The genus-two curve has the rational involution
//
//   iota: (x,y) |-> (k^2/x, k^3*y/x^3).
//
// With u=x+k^2/x, the two elliptic quotients are
//
//   E_+: v^2=(u+2k)(u+1+k^2)(u+r^2+s^2),
//   E_-: v^2=(u-2k)(u+1+k^2)(u+r^2+s^2).
//
// The (2,2)-isogeny J -> E_+ x E_- is prime to 3, so rational
// 3-torsion on J is equivalent to rational 3-torsion on one quotient.
// This script pulls the two 3-division quartics back to the fixed
// auxiliary fibers t=1/2 and t=2/3 and factors them over Q[r,u].
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
R<r,u> := PolynomialRing(Q,2);
K := FieldOfFractions(R);
PX<X> := PolynomialRing(K);

function PrimitiveNumerator(g)
    h := R!Numerator(g);
    raw := [c:c in Coefficients(h)|c ne 0];
    den := #raw eq 0 select 1 else LCM([Denominator(c):c in raw]);
    h *:= den;
    nums := [Integers()!c:c in Coefficients(h)|c ne 0];
    content := #nums eq 0 select 1 else GCD([Abs(z):z in nums]);
    return content gt 1 select h/content else h;
end function;

function Psi3ForRoots(A,B,C)
    a2 := A+B+C;
    a4 := A*B+A*C+B*C;
    a6 := A*B*C;
    return 3*u^4+4*a2*u^3+6*a4*u^2+12*a6*u+(4*a2*a6-a4^2);
end function;

function DivisionCoverGenera(AA,BB,sign)
    Kt<rr> := RationalFunctionField(Q);
    PU<uu> := PolynomialRing(Kt);
    ss := (BB-AA*rr)/(AA-BB*rr);
    kk := rr*ss;
    A := sign*2*kk; B := 1+kk^2; C := rr^2+ss^2;
    a2 := A+B+C; a4 := A*B+A*C+B*C; a6 := A*B*C;
    psi := 3*uu^4+4*a2*uu^3+6*a4*uu^2+12*a6*uu
           +(4*a2*a6-a4^2);
    assert IsIrreducible(psi);
    FF<ux> := FunctionField(psi);
    gx := Genus(FF);
    cubic := (ux+FF!A)*(ux+FF!B)*(ux+FF!C);
    if IsSquare(cubic) then return gx,gx,true; end if;
    PY<Y> := PolynomialRing(FF);
    EE<vy> := ext<FF|Y^2-cubic>;
    return gx,Genus(EE),false;
end function;

out := Open("results/target_22224_product_surface_elliptic_quotients.txt","w");
for family in [5,6] do
    AA := family eq 5 select 16 else 144;
    BB := family eq 5 select 9 else 25;
    ss := K!(BB-AA*r)/(AA-BB*r);
    kk := K!r*ss;
    common1 := 1+kk^2;
    common2 := r^2+ss^2;
    print "PRODUCT_FIBER",family,"A",AA,"B",BB,"s",ss;
    fprintf out,"FAMILY %o A=%o B=%o s=%o\n",family,AA,BB,ss;
    for sign in [1,-1] do
        psi := Psi3ForRoots(sign*2*kk,common1,common2);
        F := PrimitiveNumerator(psi);
        fac := Factorization(F);
        gx,gy,ysquare := DivisionCoverGenera(AA,BB,sign);
        print "PSI3_PULLBACK","family",family,"sign",sign,
              "total_degree",TotalDegree(F),"degree_r",Degree(F,1),
              "degree_u",Degree(F,2),"terms",#Terms(F),
              "factor_degrees",[<Degree(z[1],1),Degree(z[1],2),z[2]>:z in fac],
              "x_coordinate_cover_genus",gx,
              "rational_point_cover_genus",gy,"y_already_square",ysquare;
        fprintf out,"sign=%o psi3=%o\n",sign,F;
        fprintf out,"sign=%o factorization=%o\n",sign,fac;
        fprintf out,"sign=%o x_coordinate_cover_genus=%o\n",sign,gx;
        fprintf out,"sign=%o rational_point_cover_genus=%o y_already_square=%o\n",
            sign,gy,ysquare;
    end for;
    fprintf out,"\n";
end for;
delete out;
print "PRODUCT_ELLIPTIC_QUOTIENT_FILE results/target_22224_product_surface_elliptic_quotients.txt";
quit;
