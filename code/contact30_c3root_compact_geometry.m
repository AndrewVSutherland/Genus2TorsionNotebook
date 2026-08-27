//////////////////////////////////////////////////////////////////////
//  Compact geometry of the C3-rational-root cover of the simultaneous
//  contact-5/contact-6 order-30 family.
//
//  The old eliminated equation in (R,rho) has bidegree (20,3).  Before
//  substituting the R-parametrization, use the genus-zero core variables
//  (u,s), where q=A-e is recovered as q=N/D.  The cubic factor is
//
//      C3 = 2*x^3 + (A-3)*x^2 + (B+3)*x + (C-1),
//
//  A=(s+q)/2, B=(15-s*q)/2, C=(u^2+1)/(2*u).  Thus C3(rho)=0 is linear
//  in q and, on rho*(rho-s) != 0,
//
//      q = -[4*rho^3+(s-6)*rho^2+21*rho+(u-1)^2/u]
//            /[rho*(rho-s)].
//
//  After clearing u, the compact cover is H(u,s)=E(u,s,rho)=0 below.
//  This script factors the recovery identities and the s-resultant,
//  asks Magma for the genera of projected plane components when feasible,
//  and performs a small direct rational-height scan in u.
//
//  Typical diagnostic:
//      magma -b Height:=30 PointBound:=100 \
//          code/contact30_c3root_compact_geometry.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(6*10^9);

if not assigned Height then
    Height := 30;
elif Type(Height) eq MonStgElt then
    Height := StringToInteger(Height);
end if;
if not assigned PointBound then
    PointBound := 100;
elif Type(PointBound) eq MonStgElt then
    PointBound := StringToInteger(PointBound);
end if;
if not assigned MaxGenusDegree then
    // The genuine degree-32 component did not normalize within the short
    // diagnostic window.  Opt in explicitly if a long genus job is wanted.
    MaxGenusDegree := 12;
elif Type(MaxGenusDegree) eq MonStgElt then
    MaxGenusDegree := StringToInteger(MaxGenusDegree);
end if;

Q := Rationals();
R<u,s,rho> := PolynomialRing(Q, 3, "grevlex");
K := FieldOfFractions(R);

// Irreducible genus-zero core remaining after eliminating q.
H :=
    u^10 + 20*u^9 + 6*u^8*s + 223*u^8 + 30*u^7*s
    - 21*u^6*s^2 - 2*u^5*s^3 + 1380*u^7 - 438*u^6*s
    - 60*u^5*s^2 + 34*u^4*s^3 - 6*u^3*s^4 + 4005*u^6
    - 3525*u^5*s + 1557*u^4*s^2 - 488*u^3*s^3
    + 105*u^2*s^4 - 15*u*s^5 + s^6 + 2796*u^5
    - 2256*u^4*s + 780*u^3*s^2 - 150*u^2*s^3
    + 12*u*s^4 + 767*u^4 - 420*u^3*s + 75*u^2*s^2
    - 4*u*s^3 + 70*u^3 - 12*u^2*s - u^2;

D := u^6 + 6*u^4*s - 2*u^4 + 15*u^3*s - u*s^3 + u^2;
N := 15*u^5 + 90*u^4 + 20*u^3*s - 6*u^2*s^2 + 231*u^3
   + 2*u^2*s - 15*u*s^2 + 90*u^2 - 20*u*s + 15*u - 2*s;
qrat := K!N/K!D;

// Original two contact equations after substituting q=N/D.
F2rat := u^6*qrat - u^3*s^2*qrat^2 + 6*u^4*s*qrat - 15*u^5
       + u^4*s - 2*u^4*qrat + 30*u^3*s*qrat - 90*u^4
       + 6*u^2*s*qrat - 231*u^3 - 2*u^2*s + u^2*qrat
       - 90*u^2 - 15*u + s;
F3rat := -u^3*s*qrat^2 + u^4 + 15*u^3*qrat + u*s^2*qrat
       + 20*u^3 - 6*u^2*s + 6*u^2*qrat - 15*u*s - 20*u - 1;
G2 := R!Numerator(K!F2rat);
G3 := R!Numerator(K!F3rat);

// Cleared C3(rho)=0 after equating q=N/D to the linear root formula.
E := u*rho*(rho-s)*N
   + D*(u*(4*rho^3+(s-6)*rho^2+21*rho)+(u-1)^2);

print "COMPACT C3-root cover";
print "H total_degree", TotalDegree(H), "terms", #Terms(H);
print "D total_degree", TotalDegree(D), "terms", #Terms(D);
print "N total_degree", TotalDegree(N), "terms", #Terms(N);
print "E total_degree", TotalDegree(E), "terms", #Terms(E),
      "degree_rho", Degree(E, rho), "degree_s", Degree(E, s);

print "FACTOR H";
for fe in Factorization(H) do
    print fe[2], "degree", TotalDegree(fe[1]), "terms", #Terms(fe[1]), fe[1];
end for;
print "FACTOR recovery numerator G2";
for fe in Factorization(G2) do
    print fe[2], "degree", TotalDegree(fe[1]), "terms", #Terms(fe[1]), fe[1];
end for;
print "FACTOR recovery numerator G3";
for fe in Factorization(G3) do
    print fe[2], "degree", TotalDegree(fe[1]), "terms", #Terms(fe[1]), fe[1];
end for;
print "FACTOR E";
for fe in Factorization(E) do
    print fe[2], "degree", TotalDegree(fe[1]), "terms", #Terms(fe[1]), fe[1];
end for;

print "SPECIAL ROOT CHARTS";
print "rho=0 factorization";
for fe in Factorization(Evaluate(E, [u,s,R!0])) do print fe; end for;
print "rho=s factorization";
for fe in Factorization(Evaluate(E, [u,s,s])) do print fe; end for;

print "RESULTANT eliminating s";
time ResS := Resultant(H, E, s);
print "ResS total_degree", TotalDegree(ResS), "terms", #Terms(ResS);
resfac := Factorization(ResS);
print "ResS factors", #resfac;
for fe in resfac do
    print "RES_FACTOR", "multiplicity", fe[2],
          "total_degree", TotalDegree(fe[1]), "terms", #Terms(fe[1]);
    if #Terms(fe[1]) le 80 then
        print fe[1];
    end if;
end for;

// Plane geometry of nontrivial projected factors in (u,rho).
R2<U,Rho> := PolynomialRing(Q, 2);
toR2 := hom<R -> R2 | U, R2!0, Rho>;
print "PROJECTED COMPONENT GEOMETRY";
for fe in resfac do
    g := R2!toR2(fe[1]);
    if TotalDegree(g) le 1 then
        continue;
    end if;
    print "COMPONENT", "multiplicity", fe[2], "degree", TotalDegree(g),
          "terms", #Terms(g);
    if TotalDegree(g) gt MaxGenusDegree then
        print " geometry_skipped_degree_limit", MaxGenusDegree;
        continue;
    end if;
    try
        Caff := Curve(AffineSpace(R2), g);
        Cp := ProjectiveClosure(Caff);
        print " dimension", Dimension(Cp), "projective_degree", Degree(Cp),
              "normalized_genus", Genus(Cp);
        try
            pts := Points(Cp : Bound := PointBound);
            print " points_bound", PointBound, "count", #pts;
            if #pts le 30 then print " points", pts; end if;
        catch ep
            print " points_failed", ep`Object;
        end try;
    catch eg
        print " geometry_failed", eg`Object;
    end try;
end for;

//////////////////////////////////////////////////////////////////////
// Direct rational-height scan on the compact core.
//////////////////////////////////////////////////////////////////////

function HeightRationals(B)
    vals := [];
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num, den) eq 1 then
                Append(~vals, Q!num/Q!den);
            end if;
        end for;
    end for;
    return vals;
end function;

QS<Sv> := PolynomialRing(Q);
QR<Rv> := PolynomialRing(Q);

function CoreAtU(uv)
    return
        uv^10 + 20*uv^9 + 6*uv^8*Sv + 223*uv^8 + 30*uv^7*Sv
        - 21*uv^6*Sv^2 - 2*uv^5*Sv^3 + 1380*uv^7
        - 438*uv^6*Sv - 60*uv^5*Sv^2 + 34*uv^4*Sv^3
        - 6*uv^3*Sv^4 + 4005*uv^6 - 3525*uv^5*Sv
        + 1557*uv^4*Sv^2 - 488*uv^3*Sv^3 + 105*uv^2*Sv^4
        - 15*uv*Sv^5 + Sv^6 + 2796*uv^5 - 2256*uv^4*Sv
        + 780*uv^3*Sv^2 - 150*uv^2*Sv^3 + 12*uv*Sv^4
        + 767*uv^4 - 420*uv^3*Sv + 75*uv^2*Sv^2
        - 4*uv*Sv^3 + 70*uv^3 - 12*uv^2*Sv - uv^2;
end function;

function RecoveryAtUS(uv, sv)
    Dv := uv^6 + 6*uv^4*sv - 2*uv^4 + 15*uv^3*sv - uv*sv^3 + uv^2;
    Nv := 15*uv^5 + 90*uv^4 + 20*uv^3*sv - 6*uv^2*sv^2
        + 231*uv^3 + 2*uv^2*sv - 15*uv*sv^2 + 90*uv^2
        - 20*uv*sv + 15*uv - 2*sv;
    return Dv, Nv;
end function;

function C3AtUS(uv, sv, qv)
    A := (sv+qv)/2;
    B := (15-sv*qv)/2;
    C := (uv^2+1)/(2*uv);
    return 2*Rv^3 + (A-3)*Rv^2 + (B+3)*Rv + (C-1);
end function;

params := HeightRationals(Height);
u_checked := 0;
core_points := 0;
open_points := 0;
c3root_points := 0;
seen_core := {};
hits := [];

// Seed the four known smooth core points so the diagnostic validates its
// normalization even when Height is smaller than 125.
seed_points := [
    <Q!125,Q!5415>, <Q!125,Q!2715>,
    <Q!1/125,Q!831/3125>, <Q!1/125,Q!-69/3125>
];

for seed in seed_points do
    uv := seed[1]; sv := seed[2];
    key := Sprint(<uv,sv>);
    Include(~seen_core, key);
    core_points +:= 1;
    Dv, Nv := RecoveryAtUS(uv, sv);
    if uv ne 0 and uv^2 ne 1 and Dv ne 0 then
        open_points +:= 1;
        qv := Nv/Dv;
        roots := Roots(C3AtUS(uv,sv,qv));
        print "SEED", uv, sv, "C3_roots", roots;
        for rt in roots do
            c3root_points +:= 1;
            Append(~hits, <uv,sv,rt[1],"seed">);
        end for;
    end if;
end for;

for uv in params do
    u_checked +:= 1;
    if uv eq 0 or uv^2 eq 1 then
        continue;
    end if;
    hspec := CoreAtU(uv);
    if hspec eq 0 then
        continue;
    end if;
    for sr in Roots(hspec) do
        sv := Q!sr[1];
        key := Sprint(<uv,sv>);
        if key in seen_core then
            continue;
        end if;
        Include(~seen_core, key);
        core_points +:= 1;
        Dv, Nv := RecoveryAtUS(uv, sv);
        if Dv eq 0 then
            continue;
        end if;
        open_points +:= 1;
        qv := Nv/Dv;
        roots := Roots(C3AtUS(uv,sv,qv));
        if #roots gt 0 then
            print "C3ROOT", "u", uv, "s", sv, "q", qv, "roots", roots;
        end if;
        for rt in roots do
            c3root_points +:= 1;
            Append(~hits, <uv,sv,rt[1],"height">);
        end for;
    end for;
end for;

print "LOW_HEIGHT_SUMMARY", "Height", Height,
      "u_checked", u_checked,
      "core_points", core_points,
      "open_points", open_points,
      "c3root_points", c3root_points,
      "hits", hits;

quit;
