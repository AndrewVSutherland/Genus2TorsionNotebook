
//////////////////////////////////////////////////////////////////////
//  FAST genuine-Z/24 finder on the d=0 slice (pure polynomial test).
//
//  On d=0 (p=r(r+t)/2), split-q3-through-0 3-torsion exists over Q iff
//  Eq4(m), Eq3(m) (deg 3,5 in m; coeffs rational in r,t) have a common
//  rational root m0 with kappa=f6-m0^2 <> 0 and beta<>0.  This is a fast
//  gcd + rational-roots test -- no Jacobian, no point counting.
//
//  Purpose: is the genuine cover X={(r,t,m):Eq4=Eq3=0,kap<>0} rational
//  (=> genuine Z/24 abundant at low height) or high genus (=> sparse)?
//  Reports every genuine (r,t,m0), so their density/height answers this.
//
//  Usage: magma -b H:=80 NParts:=3 Part:=0 agent_a2_24_d0_fast.m
//////////////////////////////////////////////////////////////////////
SetColumns(0);
if not assigned MemGB then MemGB := 3; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q := Rationals(); Z := Integers();
Pm<m> := PolynomialRing(Q);

if not assigned H then H := 60; elif Type(H) eq MonStgElt then H := StringToInteger(H); end if;
if not assigned NParts then NParts := 1; elif Type(NParts) eq MonStgElt then NParts := StringToInteger(NParts); end if;
if not assigned Part then Part := 0; elif Type(Part) eq MonStgElt then Part := StringToInteger(Part); end if;
if not assigned progress then progress := 1000000; elif Type(progress) eq MonStgElt then progress := StringToInteger(progress); end if;

function HeightRationals(HH)
    vals := [];
    for den in [1..HH] do for num in [-HH..HH] do
        if GCD(num, den) eq 1 then Append(~vals, Q!num/den); end if;
    end for; end for;
    return Sort(Setseq(Seqset(vals)));
end function;

// genuine split-through-0 3-torsion roots m0 for (r,t) on d=0; [] if none
function GenuineRoots(rv, tv)
    if rv eq 0 or tv eq 0 or rv*tv eq 1 then return []; end if;
    pv := rv*(rv+tv)/2;
    e := tv^2 - 2*pv*tv/rv; lambda := rv/tv;
    A := rv^2 - lambda;
    if A eq 0 then return []; end if;
    B := 2*rv*pv - 2*lambda*(pv + rv*tv) + 2*rv*lambda;
    C := pv^2 + 2*pv*rv^2 - rv^4 - rv^3*tv - rv*pv^2/tv
         - lambda*(rv^2 + e) + 2*lambda*(rv*pv + rv^2*tv - 3*pv*tv + rv*tv^2);
    if C eq 0 then return []; end if;
    // f = (A x^2+Bx+C)(x^4 + A x^2+Bx+C); need f coeffs f1..f6
    // compute directly
    q0:=C; q1:=B; q2:=A;               // q = q2 x^2+q1 x+q0
    // W = x^4 + q ; f = q * W
    // f coeffs:
    f0 := q0*q0;
    f1 := q0*q1 + q1*q0;
    f2 := q0*q2 + q1*q1 + q2*q0;
    f3 := q1*q2 + q2*q1;               // from q*(A x^2+Bx+C) parts + q*x^4 shift
    // careful: do it by polynomial multiply
    PX := PolynomialRing(Q); X := PX.1;
    qq := q2*X^2 + q1*X + q0; ff := qq*(X^4 + qq);
    f0:=Coefficient(ff,0); f1:=Coefficient(ff,1); f2:=Coefficient(ff,2);
    f3:=Coefficient(ff,3); f4:=Coefficient(ff,4); f5:=Coefficient(ff,5); f6:=Coefficient(ff,6);
    l := C; j := f1/(2*l); n := (f2 - j^2)/(2*l);
    kap := f6 - m^2; bN := 2*m*n - f5;
    Eq4 := 3*kap*(f4 - n^2 - 2*m*j) - bN^2;
    Eq3 := 27*kap^2*(f3 - 2*m*l - 2*n*j) + bN^3;
    if Eq4 eq 0 or Eq3 eq 0 then return []; end if;
    g := GCD(Eq4, Eq3);
    if Degree(g) lt 1 then return []; end if;
    good := [];
    for rt in Roots(g) do
        m0 := rt[1];
        kv := f6 - m0^2;
        if kv eq 0 then continue; end if;             // spurious kap=0
        bv := 2*m0*n - f5; beta := bv/(3*kv);
        if beta eq 0 then continue; end if;           // q3 degenerate
        Append(~good, <m0, beta, kv>);
    end for;
    return good;
end function;

printf "d=0 FAST genuine-Z/24 finder H=%o Part=%o/%o\n", H, Part, NParts;
vals := HeightRationals(H);
tested := 0; hits := 0;
ridx := 0;
for rv in vals do
    ridx +:= 1;
    if (ridx mod NParts) ne Part then continue; end if;
    if rv eq 0 then continue; end if;
    for tv in vals do
        if tv eq 0 or rv*tv eq 1 then continue; end if;
        tested +:= 1;
        if tested mod progress eq 0 then
            printf "PROGRESS tested=%o hits=%o\n", tested, hits;
        end if;
        gr := GenuineRoots(rv, tv);
        if #gr eq 0 then continue; end if;
        hits +:= 1;
        for rec in gr do
            printf "GENUINE r=%o t=%o m0=%o beta=%o kap=%o\n", rv, tv, rec[1], rec[2], rec[3];
        end for;
    end for;
end for;
printf "SEARCH_DONE tested=%o genuine_hits=%o\n", tested, hits;
print "DONE";
quit;
