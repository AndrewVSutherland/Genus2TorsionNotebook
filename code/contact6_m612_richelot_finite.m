//////////////////////////////////////////////////////////////////////
// Finite-field diagnostic for [6,12] via the distinguished Richelot
// quotient of the contact-6 [1,2,2] factorization
//
//   f = A*B*C,
//   A = x,
//   B = (b+3)x^2 + (a-3)x + 2,
//   C = 2x^2 + (b-3)x + (a+3).
//
// We retain residues for which the source contains [6,6] and the dual
// contains [6,12].  This tests local viability of the exact proposed move.
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
prime_list := [Z!p : p in prime_list | p notin {2,3}];

function Has66(inv)
    return #[n : n in inv | (Z!n) mod 6 eq 0] ge 2;
end function;

function Has612(inv)
    return #[n : n in inv | (Z!n) mod 6 eq 0] ge 2 and
           #[n : n in inv | (Z!n) mod 12 eq 0] ge 1;
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
        return false, Parent(A)!0, "delta0";
    end if;
    g0 := Delta*Bracket(B,C)*Bracket(C,A)*Bracket(A,B);
    for sign in [1,-1] do
        g := sign*g0;
        if not Good(g) then
            continue;
        end if;
        try
            if #Jacobian(HyperellipticCurve(g)) eq source_order then
                return true,g,sign eq 1 select "+" else "-";
            end if;
        catch e
            continue;
        end try;
    end for;
    return false,Parent(A)!0,"mismatch";
end function;

print "CONTACT6_M612_RICHELOT_FINITE";
print "primes", prime_list;

for p in prime_list do
    F := GF(p);
    P<x> := PolynomialRing(F);
    good := 0;
    source66 := 0;
    source612 := 0;
    dual_good := 0;
    dual612 := 0;
    failures := 0;
    residues := [];
    samples := [];

    for a in F do
        for b in F do
            A := x;
            B := (b+3)*x^2+(a-3)*x+2;
            C := 2*x^2+(b-3)*x+(a+3);
            f := A*B*C;
            if not Good(f) then
                continue;
            end if;
            good +:= 1;
            Js := Jacobian(HyperellipticCurve(f));
            Gs, phis := AbelianGroup(Js);
            invS := Invariants(Gs);
            if not Has66(invS) then
                continue;
            end if;
            source66 +:= 1;
            if Has612(invS) then
                source612 +:= 1;
            end if;
            ok,g,sgn := DualPolynomial(A,B,C,#Js);
            if not ok then
                failures +:= 1;
                continue;
            end if;
            dual_good +:= 1;
            Gd, phid := AbelianGroup(Jacobian(HyperellipticCurve(g)));
            invD := Invariants(Gd);
            if Has612(invD) then
                dual612 +:= 1;
                Append(~residues,<Z!a,Z!b>);
                if #samples lt 8 then
                    Append(~samples,<Z!a,Z!b,invS,invD,sgn>);
                end if;
            end if;
        end for;
    end for;
    print "p",p,"good",good,"source66",source66,"source612",source612,
          "dual_good",dual_good,"dual612",dual612,
          "formula_fail",failures,"residues",residues,"samples",samples;
end for;

quit;
