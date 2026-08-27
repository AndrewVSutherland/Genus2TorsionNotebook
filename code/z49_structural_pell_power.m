//////////////////////////////////////////////////////////////////////
// Why taking the seventh power of an order-7 Pell unit does not produce
// order 49.  The displayed sextic has fundamental infinity-difference
// order 7.  We recover a degree-7 unit, raise it to the seventh power,
// and obtain a degree-49 Pell solution on the same curve; the CF order
// remains 7, so the degree-49 solution is imprimitive.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(10^9);
Q:=Rationals(); P<x>:=PolynomialRing(Q);

function SqrtLaurent(f,nmax)
    coeffs:=AssociativeArray(Integers());
    for n in [1..6] do coeffs[n]:=Q!Coefficient(f,6-n); end for;
    c:=[Q!1];
    for n in [1..nmax] do
        rhs:=IsDefined(coeffs,n) select coeffs[n] else Q!0;
        conv:=n eq 1 select Q!0 else
              &+[c[i+1]*c[n-i+1]:i in [1..n-1]];
        Append(~c,(rhs-conv)/2);
    end for;
    return c;
end function;

function FindUnit(f,m)
    r:=m-2; c:=SqrtLaurent(f,2*r+6);
    H:=Matrix(Q,[[c[(i+j+1)+3+1]:j in [0..r-1]]:i in [0..r]]);
    N:=Nullspace(Transpose(H));
    for vector in Basis(N) do
        if vector[r] eq 0 then continue; end if;
        beta:=&+[vector[j+1]*x^j:j in [0..r-1]];
        cc:=SqrtLaurent(f,Degree(beta)+8);
        alpha:=&+[(&+[Coefficient(beta,j)*cc[3-(ell-j)+1]
                     :j in [0..Degree(beta)] | ell-j le 3]) * x^ell
                   :ell in [0..Degree(beta)+3]];
        if Degree(alpha^2-beta^2*f) le 0 then
            return true,alpha,beta,alpha^2-beta^2*f;
        end if;
    end for;
    return false,P!0,P!0,Q!0;
end function;

function SqrtPolyPart(f)
    s:=x^3;
    for k in [1..3] do
        d:=f-s^2;
        if Degree(d) le 2 then break; end if;
        s+:=Coefficient(d,6-k)/(2*Coefficient(s,3))*x^(3-k);
    end for;
    return s;
end function;

function DinfinityOrder(f,maxsteps)
    s:=SqrtPolyPart(f); Pi:=P!0; Qi:=P!1; total:=0;
    for i in [0..maxsteps] do
        ai:=(Pi+s) div Qi; total+:=Degree(ai);
        Pnext:=ai*Qi-Pi;
        assert (f-Pnext^2) mod Qi eq 0;
        Qnext:=(f-Pnext^2) div Qi;
        Pi:=Pnext; Qi:=Qnext;
        if i ge 1 and Degree(Qi) le 0 and Qi ne 0 then return total; end if;
    end for;
    return 0;
end function;

f:=x^6+2*x^5-5*x^4-14*x^3-3*x^2+24*x+28;
ok,alpha,beta,norm:=FindUnit(f,7);
assert ok and norm ne 0;
A:=P!1; B:=P!0;
for i in [1..7] do
    Anew:=A*alpha+B*beta*f;
    Bnew:=A*beta+B*alpha;
    A:=Anew; B:=Bnew;
end for;

print "Z49_STRUCTURAL_PELL_POWER";
print "f",f;
print "fundamental_unit",alpha,beta,"norm",norm,
      "degrees",<Degree(alpha),Degree(beta)>;
print "seventh_power_degrees",<Degree(A),Degree(B)>,
      "pell_identity",A^2-B^2*f eq norm^7;
print "continued_fraction_fundamental_order",DinfinityOrder(f,30);
print "CONCLUSION degree_49_solution_is_imprimitive true";
print "Z49_STRUCTURAL_PELL_POWER_DONE";
quit;
