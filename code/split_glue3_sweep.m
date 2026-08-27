// BHLS (3,3)-family direct sweep: C_{a,b,c,d,t}: t y^2 = (x^3+3ax+2b)(2dx^3+3cx^2+1),
// 12ac+16bd = 1; elliptic quotients E1: t y^2 = f1, E2: t y^2 = f2 explicit.
// Funnel: torsion of E1,E2 (fast) -> genus-2 exact torsion when product big.
// Params: Sh, NSh shard on a-index; Hn, Hd heights; MinProd escalation.
SetColumns(0);
SetSeed(1);
if not assigned Sh then Sh := 0; elif Type(Sh) eq MonStgElt then Sh := StringToInteger(Sh); end if;
if not assigned NSh then NSh := 3; elif Type(NSh) eq MonStgElt then NSh := StringToInteger(NSh); end if;
if not assigned Hn then Hn := 8; elif Type(Hn) eq MonStgElt then Hn := StringToInteger(Hn); end if;
if not assigned Hd then Hd := 3; elif Type(Hd) eq MonStgElt then Hd := StringToInteger(Hd); end if;
if not assigned MinProd then MinProd := 25; elif Type(MinProd) eq MonStgElt then MinProd := StringToInteger(MinProd); end if;
SetMemoryLimit(3*10^9);
Q := Rationals();
P<x> := PolynomialRing(Q);
vals := [];
for den in [1..Hd] do for num in [-Hn..Hn] do
    if num eq 0 or GCD(Abs(num), den) ne 1 then continue; end if;
    Append(~vals, Q!num/den);
end for; end for;
Append(~vals, Q!0);
twists := [1,-1,2,-2,3,-3,5,-5,6,-6];
printf "GRID %o values, twists %o, shard %o/%o\n", #vals, #twists, Sh, NSh;
nproc := 0; nesc := 0;
seen := {};
for ia in [1..#vals] do
    if (ia mod NSh) ne Sh then continue; end if;
    a := vals[ia];
    for b in vals, c in vals do
        if b eq 0 then continue; end if;   // solve d from the relation
        d := (1 - 12*a*c)/(16*b);
        D1 := a^3 + b^2; D2 := c^3 + d^2;
        if D1 eq 0 or D2 eq 0 then continue; end if;
        f  := (x^3 + 3*a*x + 2*b)*(2*d*x^3 + 3*c*x^2 + 1);
        if not IsSquarefree(f) then continue; end if;
        f1 := x^3 + 12*(2*a^2*d - b*c)*x^2 + 12*(16*a*d^2 + 3*c^2)*D1*x + 512*D1^2*d^3;
        f2 := x^3 + 12*(2*b*c^2 - a*d)*x^2 + 12*(16*b^2*c + 3*a^2)*D2*x + 512*D2^2*b^3;
        nproc +:= 1;
        for t in twists do
            okE := true;
            try
                E1 := QuadraticTwist(EllipticCurve(f1), t);
                E2 := QuadraticTwist(EllipticCurve(f2), t);
            catch e okE := false; end try;
            if not okE then continue; end if;
            t1 := #TorsionSubgroup(E1);
            if t1 lt 3 then continue; end if;
            t2 := #TorsionSubgroup(E2);
            if t1*t2 lt MinProd then continue; end if;
            // escalate
            C := HyperellipticCurve(t*f);
            ii := G2Invariants(C);
            if ii in seen then continue; end if;
            Include(~seen, ii);
            nesc +:= 1;
            Cm := C;
            try Cm := ReducedMinimalWeierstrassModel(C); catch e ; end try;
            T := TorsionSubgroup(Jacobian(SimplifiedModel(Cm)));
            printf "G3HIT a=%o b=%o c=%o t=%o E1tors=%o E2tors=%o JTORSION %o order %o\n",
                a, b, c, t, t1, t2, Invariants(T), #T;
        end for;
    end for;
    printf "PROGRESS ia=%o/%o processed %o escalated %o\n", ia, #vals, nproc, nesc;
end for;
printf "GLUE3_DONE shard %o/%o\n", Sh, NSh;
quit;
