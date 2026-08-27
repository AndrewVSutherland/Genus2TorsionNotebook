//////////////////////////////////////////////////////////////////////
// Class-by-class finite halving diagnostic for the contact-6 [6,6]
// source and its distinguished Richelot dual.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Z := Integers();

if assigned primes then
    if Type(primes) eq MonStgElt then
        prime_list := [StringToInteger(s) : s in Split(primes, ",") | #s gt 0];
    else
        prime_list := primes;
    end if;
else
    prime_list := [5,7,11,13,17,19,23,29,31];
end if;

function Has66(inv)
    return #[n : n in inv | (Z!n) mod 6 eq 0] ge 2;
end function;

function Good(f)
    return Degree(f) in {5,6} and Discriminant(f) ne 0;
end function;

function Bracket(A,B)
    return Derivative(A)*B-A*Derivative(B);
end function;

function CoefficientMatrixDeterminant(A,B,C)
    F := BaseRing(Parent(A));
    M := Matrix(F,3,3,
        [Coefficient(A,i) : i in [0..2]] cat
        [Coefficient(B,i) : i in [0..2]] cat
        [Coefficient(C,i) : i in [0..2]]);
    return Determinant(M);
end function;

function DualPolynomial(A,B,C,source_order)
    Delta := CoefficientMatrixDeterminant(A,B,C);
    if Delta eq 0 then
        return false, Parent(A)!0;
    end if;
    g0 := Delta*Bracket(B,C)*Bracket(C,A)*Bracket(A,B);
    for sign in [1,-1] do
        g := sign*g0;
        if Good(g) and #Jacobian(HyperellipticCurve(g)) eq source_order then
            return true,g;
        end if;
    end for;
    return false,Parent(A)!0;
end function;

procedure AddMask(~counts,mask)
    if IsDefined(counts,mask) then
        counts[mask] +:= 1;
    else
        counts[mask] := 1;
    end if;
end procedure;

function IsDivisibleBy2Finite(D,G,phi)
    a := D @@ phi;
    coords := Eltseq(a);
    invs := Invariants(G);
    for i in [1..#coords] do
        if GCD(2,invs[i]) ne 1 and coords[i] mod 2 ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;

function HalfMask(J,polys,G,phi)
    mask := 0;
    bit := 1;
    for q in polys do
        T := J![q/LeadingCoefficient(q),0];
        if IsDivisibleBy2Finite(T,G,phi) then
            mask +:= bit;
        end if;
        bit *:= 2;
    end for;
    return mask;
end function;

print "CONTACT6_M612_CLASS_FINITE";
print "mask bits source=(T0,TB,TC), dual=(R1,R2,R3)";

for p in prime_list do
    F := GF(p);
    P<x> := PolynomialRing(F);
    src_counts := AssociativeArray();
    dual_counts := AssociativeArray();
    src66 := 0;
    dual_tested := 0;

    for a in F do
        for b in F do
            A := x;
            B := (b+3)*x^2+(a-3)*x+2;
            C := 2*x^2+(b-3)*x+(a+3);
            f := A*B*C;
            if not Good(f) then continue; end if;
            Js := Jacobian(HyperellipticCurve(f));
            Gs, phis := AbelianGroup(Js);
            if not Has66(Invariants(Gs)) then continue; end if;
            src66 +:= 1;
            AddMask(~src_counts,HalfMask(Js,[A,B,C],Gs,phis));

            ok,g := DualPolynomial(A,B,C,#Js);
            if not ok then continue; end if;
            Jd := Jacobian(HyperellipticCurve(g));
            Gd, phid := AbelianGroup(Jd);
            R1 := Bracket(B,C);
            R2 := Bracket(C,A);
            R3 := Bracket(A,B);
            dual_tested +:= 1;
            AddMask(~dual_counts,HalfMask(Jd,[R1,R2,R3],Gd,phid));
        end for;
    end for;

    print "p",p,"source66",src66,
          "source_masks",Sort([<k,src_counts[k]> : k in Keys(src_counts)]),
          "dual_tested",dual_tested,
          "dual_masks",Sort([<k,dual_counts[k]> : k in Keys(dual_counts)]);
end for;

quit;
