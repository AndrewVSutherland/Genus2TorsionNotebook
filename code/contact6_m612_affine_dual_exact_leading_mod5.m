//////////////////////////////////////////////////////////////////////
// Exact dual-halving covariants on the surviving affine mod-5 cones.
//
// The corrected affine blowup audit leaves
//
//   DB, U+2v : 8 square-cover-smooth cones,
//   DC, U+2v : 8 square-cover-smooth cones,
//
// all over (a,b)=(0,0).  R3 is the relevant dual class on DB=W^2;
// R2 is the relevant dual class on DC=W^2.  This file substitutes the
// common special point into the exact square-quartic halving covariants
// and audits the valid leading-coefficient chart.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
k:=GF(5);
A<m,n>:=PolynomialRing(k,2,"grevlex");
PX<X>:=PolynomialRing(A);

function SquareQuarticEquations(S)
    s4:=Coefficient(S,4); s3:=Coefficient(S,3);
    s2:=Coefficient(S,2); s1:=Coefficient(S,1);
    s0:=Coefficient(S,0);
    E1:=8*s4^2*s1-s3*(4*s4*s2-s3^2);
    E0:=64*s4^3*s0-(4*s4*s2-s3^2)^2;
    return A!E1,A!E0,s4;
end function;

// At a=b=0 modulo 5, Delta=0 and all three dual quadratics coincide.
Delta:=A!0;
R1:=2*(X^2+1);
R2:=2*(X^2+1);
R3:=2*(X^2+1);

procedure Audit(label,u,v,w,number_of_core_cones)
    S:=u*(m*X+n)^2-Delta*v*w;
    E1,E0,s4:=SquareQuarticEquations(S);

    assert E1 eq 4*m^5*n;
    assert E0 eq -m^6*(m^2+n^2);
    assert s4 eq 2*m^2;

    all:=[]; valid:=[]; degenerate:=[]; degenerate_ranks:=[];
    dE:=[[Derivative(E1,j):j in [1,2]],
         [Derivative(E0,j):j in [1,2]]];
    for m0 in k do for n0 in k do
        pt:=[m0,n0];
        if Evaluate(E1,pt) ne 0 or Evaluate(E0,pt) ne 0 then
            continue;
        end if;
        Append(~all,<Integers()!m0,Integers()!n0>);
        if Evaluate(s4,pt) ne 0 then
            Append(~valid,<Integers()!m0,Integers()!n0>);
        else
            Append(~degenerate,<Integers()!m0,Integers()!n0>);
            J:=Matrix(k,2,2,
                [Evaluate(dE[i][j],pt):i in [1,2],j in [1,2]]);
            Append(~degenerate_ranks,Rank(J));
        end if;
    end for; end for;

    print "CLASS",label,"CORE_SQUARE_CONES",number_of_core_cones;
    print "SPECIALIZED_R1_R2_R3",R1;
    print "E1",E1,"E0",E0,"LEADING_COEFFICIENT_s4",s4;
    print "COVARIANT_POINTS",all;
    print "VALID_QUARTIC_CHART_POINTS_s4_nonzero",valid;
    print "DEGENERATE_m0_POINTS",degenerate,
          "COVARIANT_JACOBIAN_RANKS",degenerate_ranks;
    print "CERTIFIED_SMOOTH_EXACT_HALVING_CONES",0;
    print "UNRESOLVED_DEGENERATE_HALVING_BOUNDARY_PER_CORE_CONE",#degenerate;
end procedure;

print "CONTACT6_M612_AFFINE_DUAL_EXACT_LEADING_MOD5";
Audit("R3_on_DB_square",R3,R1,R2,8);
Audit("R2_on_DC_square",R2,R1,R3,8);
print "CONCLUSION valid exact square-quartic chart empty for all 16 cones";
print "CAUTION m=0 has s4=0 and rank 0; it is unresolved, not a half";
print "CONTACT6_M612_AFFINE_DUAL_EXACT_LEADING_MOD5_DONE";
quit;
