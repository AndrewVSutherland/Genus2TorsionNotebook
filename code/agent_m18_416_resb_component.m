//////////////////////////////////////////////////////////////////////
//  Projected branch components for the [4,16] second-halving surface
//  Sigma' : E1core = E0core = 0 in A^4_{R,w,a,b}.
//
//  Continuation of agent_m18_416_branch_discriminant.m.  That script
//  found, modulo p=101 and p=103,
//
//      Res_b(E1core,E0core)
//        = (w-1)^3 (w+1)^3 (R+1)^8 (degree 5)^8 (degree 52)
//
//  in A^3_{R,w,a}.  The purpose here is to isolate the repeated degree-5
//  projected component, preferably over Q, and print enough data to decide
//  whether it is a genuine rational component rather than a mod-p accident.
//
//  Usage:
//    magma -b do_rational:=true primes:="101,103,107" \
//        code/agent_m18_416_resb_component.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned do_rational then
    do_rational := true;
elif Type(do_rational) eq MonStgElt then
    do_rational := do_rational in {"true", "True", "1", "yes"};
end if;
if not assigned primes then
    prime_list := [101,103,107,109];
elif Type(primes) eq MonStgElt then
    prime_list := [StringToInteger(s) : s in Split(primes, ",") | #s gt 0];
else
    prime_list := primes;
end if;
if not assigned print_large_factor then
    print_large_factor := false;
elif Type(print_large_factor) eq MonStgElt then
    print_large_factor := print_large_factor in {"true", "True", "1", "yes"};
end if;

Q := Rationals();
Z := Integers();
R4<R,w,a,b> := PolynomialRing(Q, 4, "grevlex");
K4 := FieldOfFractions(R4);
PX<x> := PolynomialRing(K4);

t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
A := x^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*x + R^4;
B := (R + 2 + 4*t)*x^2 + (R^2 + 4*R + 1 + 8*t)*x + (2*R^2 + R + 4*t);
f := x*A*B;
c4 := R + 2 + 4*t;
q := x^2 + a*x + b;
D := f - c4*(x+R)*q^2;
d4 := Coefficient(D,4); d3 := Coefficient(D,3); d2 := Coefficient(D,2);
d1 := Coefficient(D,1); d0 := Coefficient(D,0);
E1 := 8*d4^2*d1 - d3*(4*d4*d2 - d3^2);
E0 := 64*d4^3*d0 - (4*d4*d2 - d3^2)^2;

function ClearPrim4(g)
    gn := Numerator(g);
    den := LCM([Denominator(co) : co in Coefficients(gn)]);
    gn := den*gn;
    gg := GCD([Z!co : co in Coefficients(gn)]);
    return R4!(gn/gg);
end function;

E1n := ClearPrim4(E1);
E0n := ClearPrim4(E0);
d4n := ClearPrim4(d4);
while IsDivisibleBy(E1n, R-1) do E1n := ExactQuotient(E1n, R-1); end while;
while IsDivisibleBy(E0n, R-1) do E0n := ExactQuotient(E0n, R-1); end while;
d4core4 := d4n;
while IsDivisibleBy(d4core4, R-1) do d4core4 := ExactQuotient(d4core4, R-1); end while;

printf "E1core deg %o terms %o | E0core deg %o terms %o\n",
    TotalDegree(E1n), #Terms(E1n), TotalDegree(E0n), #Terms(E0n);
printf "d4 cleared factorization:\n";
for fr in Factorization(d4n) do
    printf "  <deg %o, terms %o, mult %o>  %o\n",
        TotalDegree(fr[1]), #Terms(fr[1]), fr[2], fr[1];
end for;
printf "d4core after stripping R-1: deg %o terms %o  %o\n",
    TotalDegree(d4core4), #Terms(d4core4), d4core4;

R3<R_,w_,a_> := PolynomialRing(Q, 3, "grevlex");
Pb := PolynomialRing(R3);
h34 := hom<R4 -> Pb | [Pb!R_, Pb!w_, Pb!a_, Pb.1]>;
E1b := h34(E1n);
E0b := h34(E0n);
d4core := R3!Coefficient(h34(d4core4), 0);
printf "deg_b E1core=%o deg_b E0core=%o\n", Degree(E1b), Degree(E0b);

function Primitive3(g)
    if g eq 0 then return g; end if;
    den := LCM([Denominator(co) : co in Coefficients(g)]);
    gn := R3!(den*g);
    gg := GCD([Z!co : co in Coefficients(gn)]);
    if gg ne 0 then gn := R3!(gn/gg); end if;
    return gn;
end function;

procedure PrintFactorProfile(fac, label)
    printf "%o factor profile:\n", label;
    for fr in fac do
        g := fr[1]; m := fr[2];
        printf "  <deg %o, terms %o, mult %o>", TotalDegree(g), #Terms(g), m;
        if #Terms(g) le 80 or print_large_factor then
            printf "  %o", g;
        end if;
        printf "\n";
    end for;
end procedure;

if do_rational then
    tim := Cputime();
    ResQ := Primitive3(Resultant(E1b, E0b));
    printf "Q: Res_b degree %o terms %o (%.1o s)\n",
        TotalDegree(ResQ), #Terms(ResQ), Cputime(tim);
    tim := Cputime();
    facQ := Factorization(ResQ);
    printf "Q: factorization completed in %.1o s\n", Cputime(tim);
    PrintFactorProfile(facQ, "Q Res_b");
    printf "Q: d4core^8 divides Res_b? %o\n", IsDivisibleBy(ResQ, d4core^8);
    if IsDivisibleBy(ResQ, d4core^8) then
        sat := ExactQuotient(ResQ, d4core^8);
        printf "Q: Res_b / d4core^8 degree %o terms %o\n",
            TotalDegree(sat), #Terms(sat);
        PrintFactorProfile(Factorization(sat), "Q d4-saturated Res_b");
        geom := sat;
        for bd in [w_-1, w_+1, R_+1] do
            while IsDivisibleBy(geom, bd) do geom := ExactQuotient(geom, bd); end while;
        end for;
        printf "Q: after also stripping w=+-1 and R=-1: degree %o terms %o\n",
            TotalDegree(geom), #Terms(geom);
        facGeom := Factorization(geom);
        PrintFactorProfile(facGeom, "Q geometric residual projection");
        for fr in facGeom do
            g := fr[1];
            printf "Q_GEOM_FACTOR deg %o mult %o terms %o degrees(R,w,a)=(%o,%o,%o)\n",
                TotalDegree(g), fr[2], #Terms(g), Degree(g,1), Degree(g,2), Degree(g,3);
        end for;
    end if;

    repeated_low := [fr[1] : fr in facQ | TotalDegree(fr[1]) le 8 and fr[2] ge 2];
    printf "Q: repeated low-degree factors count %o\n", #repeated_low;
    for g in repeated_low do
        printf "Q_REPEATED_FACTOR deg %o terms %o\n%o\n",
            TotalDegree(g), #Terms(g), g;
    end for;

    core := ResQ;
    for bd in [w_-1, w_+1, R_+1] do
        while IsDivisibleBy(core, bd) do core := ExactQuotient(core, bd); end while;
    end for;
    printf "Q: after stripping (w-1),(w+1),(R+1): degree %o terms %o\n",
        TotalDegree(core), #Terms(core);
    facCore := Factorization(core);
    PrintFactorProfile(facCore, "Q stripped core");
end if;

for pp in prime_list do
    Fp := GF(pp);
    R3p<Rp,wp,ap> := PolynomialRing(Fp, 3, "grevlex");
    Pbp := PolynomialRing(R3p);
    hp := hom<R4 -> Pbp | [Pbp!Rp, Pbp!wp, Pbp!ap, Pbp.1]>;
    tim := Cputime();
    Resp := Resultant(hp(E1n), hp(E0n));
    printf "p=%o: Res_b degree %o terms %o (%.1o s)\n",
        pp, TotalDegree(Resp), #Terms(Resp), Cputime(tim);
    tim := Cputime();
    facp := Factorization(Resp);
    printf "p=%o: factorization completed in %.1o s\n", pp, Cputime(tim);
    PrintFactorProfile(facp, Sprintf("p=%o Res_b", pp));
    deg5 := [fr : fr in facp | TotalDegree(fr[1]) eq 5 and fr[2] ge 2];
    printf "p=%o: repeated degree-5 factors %o\n", pp, #deg5;
    for fr in deg5 do
        printf "P%o_DEG5_MULT%o terms %o\n%o\n", pp, fr[2], #Terms(fr[1]), fr[1];
    end for;
end for;

print "DONE";
quit;
