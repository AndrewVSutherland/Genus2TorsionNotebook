//////////////////////////////////////////////////////////////////////
// Necessary finite-group masks on the transverse HLP deformation
//   F_t=F_HLP+t*(1+x).
// A rational [6,12] subgroup injects at every good p>3, so the finite
// Jacobian must contain two invariant factors divisible by 6 and one by 12.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Z:=Integers();

function Has612(inv)
    return #[n:n in inv|(Z!n) mod 6 eq 0] ge 2 and
           #[n:n in inv|(Z!n) mod 12 eq 0] ge 1;
end function;

print "M612_HLP_TRANSVERSE_FINITE";
for p in [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89] do
    k:=GF(p); P<x>:=PolynomialRing(k);
    fseed:=187392*x^6-118767*x^4-118767*x^2+187392;
    good:=[]; allowed:=[]; bad:=[]; records:=[];
    for t in k do
        f:=fseed+t*(1+x);
        if Degree(f) notin {5,6} or Discriminant(f) eq 0 then
            Append(~bad,Z!t); continue;
        end if;
        Append(~good,Z!t);
        try
            A,mp:=AbelianGroup(Jacobian(HyperellipticCurve(f)));
            inv:=Invariants(A);
            if Has612(inv) then
                Append(~allowed,Z!t); Append(~records,<Z!t,inv>);
            end if;
        catch e
            Append(~bad,Z!t);
        end try;
    end for;
    print "p",p,"good",#good,"bad",Sort(bad),
          "allowed",Sort(allowed),"records",records;
end for;
quit;
