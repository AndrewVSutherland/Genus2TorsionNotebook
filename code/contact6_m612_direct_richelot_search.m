//////////////////////////////////////////////////////////////////////
// Boundary-aware direct search for [6,12] in the contact-6 (a,b) plane.
//
// Two targets are tested simultaneously:
//   (S) the source Jacobian itself has exact torsion [6,12];
//   (D) the source contains [6,6] and its distinguished rational
//       Richelot codomain has exact torsion [6,12].
//
// At each finite prime we retain bad/boundary reductions conservatively.
// On good reductions we use exact finite Jacobian invariant factors.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

if not assigned height then height := 30;
elif Type(height) eq MonStgElt then height := StringToInteger(height); end if;
if not assigned prime_bound then prime_bound := 31;
elif Type(prime_bound) eq MonStgElt then prime_bound := StringToInteger(prime_bound); end if;
if not assigned max_exact then max_exact := 200;
elif Type(max_exact) eq MonStgElt then max_exact := StringToInteger(max_exact); end if;
if not assigned progress_interval then progress_interval := 100000;
elif Type(progress_interval) eq MonStgElt then progress_interval := StringToInteger(progress_interval); end if;

function HeightRationals(B)
    vals := []; seen := {};
    for d in [1..B] do
        for n in [-B..B] do
            if GCD(n,d) ne 1 then continue; end if;
            q := Q!n/d; k := Sprint(q);
            if k notin seen then Include(~seen,k); Append(~vals,q); end if;
        end for;
    end for;
    return vals;
end function;

function Has66(inv)
    return #[n : n in inv | (Z!n) mod 6 eq 0] ge 2;
end function;

function Has612(inv)
    return #[n : n in inv | (Z!n) mod 6 eq 0] ge 2 and
           #[n : n in inv | (Z!n) mod 12 eq 0] ge 1;
end function;

function ContactFactors(K,a,b)
    PK<X> := PolynomialRing(K);
    A := X;
    B := (b+3)*X^2+(a-3)*X+2;
    C := 2*X^2+(b-3)*X+(a+3);
    return A,B,C,A*B*C;
end function;

function Good(f)
    return Degree(f) in {5,6} and Discriminant(f) ne 0;
end function;

function Bracket(A,B)
    return Derivative(A)*B-A*Derivative(B);
end function;

function CoeffDet(A,B,C)
    K := BaseRing(Parent(A));
    M := Matrix(K,3,3,
        [Coefficient(A,i) : i in [0..2]] cat
        [Coefficient(B,i) : i in [0..2]] cat
        [Coefficient(C,i) : i in [0..2]]);
    return Determinant(M);
end function;

function DualFinite(A,B,C,source_order)
    Delta := CoeffDet(A,B,C);
    if Delta eq 0 then return false,Parent(A)!0; end if;
    g0 := Delta*Bracket(B,C)*Bracket(C,A)*Bracket(A,B);
    for sgn in [1,-1] do
        g := sgn*g0;
        if not Good(g) then continue; end if;
        try
            if #Jacobian(HyperellipticCurve(g)) eq source_order then
                return true,g;
            end if;
        catch e
            continue;
        end try;
    end for;
    return false,Parent(A)!0;
end function;

function ResidueTables(p)
    F := GF(p);
    goodS := {Z|}; allowS := {Z|};
    goodD := {Z|}; allowD := {Z|};
    for aa in [0..p-1] do
        for bb in [0..p-1] do
            a := F!aa; b := F!bb; key := aa+p*bb;
            A,B,C,f := ContactFactors(F,a,b);
            if not Good(f) then continue; end if;
            Include(~goodS,key);
            Js := Jacobian(HyperellipticCurve(f));
            Gs,mp := AbelianGroup(Js); invS := Invariants(Gs);
            if Has612(invS) then Include(~allowS,key); end if;
            if not Has66(invS) then
                // The dual cannot have rational (3,3) if the good source does not.
                Include(~goodD,key);
                continue;
            end if;
            ok,g := DualFinite(A,B,C,#Js);
            if not ok then
                // Degenerate dual reduction is a boundary and is retained.
                continue;
            end if;
            Include(~goodD,key);
            Gd,mpd := AbelianGroup(Jacobian(HyperellipticCurve(g)));
            if Has612(Invariants(Gd)) then Include(~allowD,key); end if;
        end for;
    end for;
    return goodS,allowS,goodD,allowD;
end function;

function ReduceKey(a,b,p)
    F := GF(p);
    try
        aa := F!a; bb := F!b;
    catch e
        return false,0;
    end try;
    return true,Z!aa+p*(Z!bb);
end function;

function IntegralScale(f)
    d := LCM([Denominator(Coefficient(f,i)) : i in [0..Degree(f)]]);
    return P!(d^2*f);
end function;

function TorsionData(C)
    f,h := HyperellipticPolynomials(C);
    if Degree(h) ge 0 then return false,[],P!0; end if;
    fI := IntegralScale(P!f);
    CI := HyperellipticCurve(fI); JI := Jacobian(CI);
    G,mp := TorsionSubgroup(JI);
    return true,Invariants(G),fI;
end function;

function SimpleCert(f)
    C := HyperellipticCurve(f);
    for p in [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97] do
        try
            fp := ChangeRing(f,GF(p));
            if Degree(fp) ne Degree(f) or Discriminant(fp) eq 0 then continue; end if;
            Lp := LPolynomial(ChangeRing(C,GF(p))); fac := Factorization(Lp);
            if #fac eq 1 and fac[1][2] eq 1 and Degree(fac[1][1]) eq 4 then
                return true,p,Lp;
            end if;
        catch e
            continue;
        end try;
    end for;
    return false,0,P!0;
end function;

params := HeightRationals(height);
primes := [p : p in PrimesUpTo(prime_bound) | p notin {2,3}];
tables := [];

print "CONTACT6_M612_DIRECT_RICHELOT_SEARCH";
print "height",height,"parameter_count",#params,"prime_bound",prime_bound,
      "max_exact",max_exact;
print "PRECOMPUTE";
for p in primes do
    goodS,allowS,goodD,allowD := ResidueTables(p);
    Append(~tables,<p,goodS,allowS,goodD,allowD>);
    print "p",p,"good_source",#goodS,"allow_source612",#allowS,
          "good_dual",#goodD,"allow_dual612",#allowD;
end for;

checked := 0; mask_survivors := 0; open122 := 0; exact := 0; hits := 0;
for a in params do
    for b in params do
        checked +:= 1;
        if progress_interval gt 0 and checked mod progress_interval eq 0 then
            print "PROGRESS",checked,"mask",mask_survivors,"open122",open122,
                  "exact",exact,"hits",hits;
        end if;
        possibleS := true; possibleD := true;
        for T in tables do
            p := T[1]; ok,key := ReduceKey(a,b,p);
            if not ok then continue; end if;
            if key in T[2] then
                if key notin T[3] then possibleS := false; end if;
                if key in T[4] and key notin T[5] then possibleD := false; end if;
            end if;
            if not possibleS and not possibleD then break; end if;
        end for;
        if not possibleS and not possibleD then continue; end if;
        mask_survivors +:= 1;

        A,B,C,f0 := ContactFactors(Q,a,b); f := P!f0;
        if not Good(f) then continue; end if;
        if Sort([Degree(z[1]) : z in Factorization(f)]) ne [1,2,2] then continue; end if;
        open122 +:= 1;
        simple,pcert,Lp := SimpleCert(f);
        if not simple then continue; end if;
        if exact ge max_exact then continue; end if;

        CS := HyperellipticCurve(f); JS := Jacobian(CS);
        okS,invS,fSI := TorsionData(CS); exact +:= 1;
        print "EXACT_SOURCE",a,b,"possibleS",possibleS,"possibleD",possibleD,
              "torsion",invS,"pcert",pcert;
        if possibleS and invS eq [6,12] then
            hits +:= 1; print "HIT_SOURCE_6_12",a,b,fSI;
        end if;
        if possibleD and Has66(invS) then
            try
                Rs := RichelotIsogenousSurfaces(JS);
                for i in [1..#Rs] do
                    if Type(Rs[i]) ne JacHyp then continue; end if;
                    okD,invD,fDI := TorsionData(Curve(Rs[i]));
                    print "EXACT_DUAL",a,b,i,"torsion",invD;
                    if invD eq [6,12] then
                        hits +:= 1; print "HIT_DUAL_6_12",a,b,i,fDI;
                    end if;
                end for;
            catch e
                print "DUAL_ERROR",a,b,e`Object;
            end try;
        end if;
    end for;
end for;

print "DONE","checked",checked,"mask_survivors",mask_survivors,
      "open122",open122,"exact",exact,"hits",hits;
quit;
