//////////////////////////////////////////////////////////////////////
// Larger-prime local scan for the R = -8 full second-stage conditions.
//
// This is the fixed-fiber version of agent_m18_416_fiber_local_scan.m:
// enumerate m residue classes modulo p^2 (plus valuation +/-1 classes)
// and test the simultaneous S_A and S_B p-adic square conditions.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();

if not assigned MaxPrime then MaxPrime := 101; end if;
if Type(MaxPrime) eq MonStgElt then MaxPrime := StringToInteger(MaxPrime); end if;
if not assigned Prec then Prec := 80; end if;
if Type(Prec) eq MonStgElt then Prec := StringToInteger(Prec); end if;

K2<R,m> := RationalFunctionField(Q, 2);
PX<x> := PolynomialRing(K2);
K := -2*R*(R^2-1);
w := (m^2 + K)/(m^2 - K);
t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
c4 := R + 2 + 4*t;
Apol := x^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*x + R^4;
Bpol := c4*x^2 + (R^2 + 4*R + 1 + 8*t)*x + (2*R^2 + R + 4*t);
XR := -c4*R;
At := PX![K2!co : co in Coefficients(c4^2*Evaluate(Apol, x/c4))];
Bt := PX![K2!co : co in Coefficients(c4*Evaluate(Bpol, x/c4))];
G := 2*(R^2-1)*(R*(2*R+1) - w^2*(R+2));
alphaA := XR + Coefficient(At,1)/2;
alphaB := XR + Coefficient(Bt,1)/2;
okA, hA := IsSquare(Evaluate(At, XR)/G); assert okA;
okB, hB := IsSquare(Evaluate(Bt, XR)/G); assert okB;

funcs := [G, alphaA, hA, alphaB, hB];
fnum := [Numerator(f) : f in funcs];
fden := [Denominator(f) : f in funcs];
Rv := Q!-8;

function EvalAll(mv)
    vals := [];
    for i in [1..#funcs] do
        dd := Evaluate(fden[i], [Rv, mv]);
        if dd eq 0 then
            return false, [];
        end if;
        Append(~vals, Evaluate(fnum[i], [Rv, mv])/dd);
    end for;
    return true, vals;
end function;

function PointPasses(mv, Qp)
    ok, v := EvalAll(mv);
    if not ok then
        return false;
    end if;
    Gv := v[1]; aA := v[2]; hAv := v[3]; aB := v[4]; hBv := v[5];
    if Gv eq 0 or hAv eq 0 or hBv eq 0 then
        return false;
    end if;
    okG, sg := IsSquare(Qp!Gv);
    if not okG then
        return false;
    end if;
    passA := false;
    for sgn in [1,-1] do
        VA := Qp!(2*aA) + sgn*(Qp!(2*hAv))*sg;
        if VA ne 0 and IsSquare(VA) then
            passA := true;
            break;
        end if;
    end for;
    if not passA then
        return false;
    end if;
    for sgn in [1,-1] do
        VB := Qp!(2*aB) + sgn*(Qp!(2*hBv))*sg;
        if VB ne 0 and IsSquare(VB) then
            return true;
        end if;
    end for;
    return false;
end function;

function RepsForPrime(p)
    reps := [];
    if p eq 2 then
        for a in [1..31 by 2] do Append(~reps, Q!a); end for;
        for a in [1..15] do Append(~reps, Q!(2*a)); end for;
        for a in [1..31 by 2] do Append(~reps, Q!a/2); end for;
        for a in [1..31 by 2] do Append(~reps, Q!a/4); end for;
        for a in [1..15 by 2] do Append(~reps, Q!(4*a)); end for;
    else
        for a in [1..p-1] do
            for b in [0..p-1] do
                Append(~reps, Q!(a + b*p));
            end for;
        end for;
        for c in [1..p-1] do Append(~reps, Q!c/p); end for;
        for c in [1..p-1] do Append(~reps, Q!(c*p)); end for;
    end if;
    return reps;
end function;

printf "R=-8 local second-stage scan, MaxPrime=%o Prec=%o\n", MaxPrime, Prec;
for p in PrimesUpTo(MaxPrime) do
    Qp := pAdicField(p, Prec);
    reps := RepsForPrime(p);
    found := false;
    witness := Q!0;
    tested := 0;
    for mv in reps do
        tested +:= 1;
        if PointPasses(mv, Qp) then
            found := true;
            witness := mv;
            break;
        end if;
    end for;
    if found then
        printf "p=%o soluble, witness m=%o, tested_until=%o/%o\n", p, witness, tested, #reps;
    else
        printf "LOCAL OBSTRUCTION at p=%o after testing %o classes\n", p, #reps;
        break;
    end if;
end for;
print "DONE";
quit;
