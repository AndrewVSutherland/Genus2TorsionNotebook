//////////////////////////////////////////////////////////////////////
//  Genus diagnostic for HPL C-square + F1 line slices.
//
//  For each small integer direction (ds,dm,dq), form the HPL slice from
//  code/m2248_hpl_c_f1_slice_search.m.  The first remaining geometric
//  square condition is D=d^2.  This script computes the squarefree
//  polynomial defining y^2 = D(u) and reports the lowest-degree slices.
//
//  A degree 3 or 4 squarefree part would give a genus-one slice worth
//  studying as an elliptic curve.  High degree means the naive line slice
//  is unlikely to be a useful fibration.
//
//  Typical run:
//      magma -b dir_bound:=4 top:=20 code/m2248_hpl_slice_genus_diagnostic.m
//////////////////////////////////////////////////////////////////////

if not assigned dir_bound then
    dir_bound := 3;
elif Type(dir_bound) eq MonStgElt then
    dir_bound := StringToInteger(dir_bound);
end if;

if not assigned top then
    top := 20;
elif Type(top) eq MonStgElt then
    top := StringToInteger(top);
end if;

if not assigned show_polys then
    show_polys := false;
elif Type(show_polys) eq MonStgElt then
    show_polys := show_polys eq "true" or show_polys eq "True" or show_polys eq "1";
end if;

Q := Rationals();
R<u> := PolynomialRing(Q);
K<U> := FunctionField(Q);

function SquarefreePartPolynomial(f)
    if f eq 0 then
        return R!0;
    end if;
    fac := Factorization(f);
    out := R!1;
    for ff in fac do
        if ff[2] mod 2 eq 1 then
            out *:= ff[1];
        end if;
    end for;
    return out;
end function;

function ConicParam(s, m)
    den := 1 + (s^2 - 1)*m^2;
    r := s*((s^2 - 1)*m^2 - 1)/den;
    c := -2*s*m/den;
    return r, c;
end function;

function ToPolynomial(f)
    return R![ Q!Coefficient(f, i) : i in [0..Degree(f)] ];
end function;

rho0 := Q!58466134224 / Q!53109477625;
sigma0 := Q!719363573659505664 / Q!749082246897952705;
tau0 := Q!307598400 / Q!352612321;
F10 := (1 + rho0)*(1 + sigma0)*(1 + tau0);
ok10, q10 := IsSquare(F10);
assert ok10;
m20 := (1 + rho0/sigma0)/((sigma0^2 - 1)*(1 - rho0/sigma0));
okm, m0 := IsSquare(m20);
assert okm;

rows := [];

for ds in [-dir_bound..dir_bound] do
    for dm in [-dir_bound..dir_bound] do
        for dq in [-dir_bound..dir_bound] do
            if ds eq 0 and dm eq 0 and dq eq 0 then
                continue;
            end if;

            sigma := K!sigma0*(1 + (K!ds)*U);
            m := K!m0*(1 + (K!dm)*U);
            q1 := K!q10*(1 + (K!dq)*U);
            rho, c := ConicParam(sigma, m);
            tau := q1^2/((1 + rho)*(1 + sigma)) - 1;
            Dval := (tau^2 - rho^2)/(tau^2 - 1);

            num := ToPolynomial(Numerator(Dval));
            den := ToPolynomial(Denominator(Dval));
            sf := SquarefreePartPolynomial(num*den);
            deg := Degree(sf);
            genus := deg le 0 select -1 else Floor((deg - 1)/2);
            Append(~rows, <deg, genus, [ds,dm,dq], sf>);
        end for;
    end for;
end for;

Sort(~rows);
print "HPL C-square + F1 slice D=d^2 genus diagnostic";
print "dir_bound", dir_bound, "directions", #rows;
for i in [1..Minimum(top,#rows)] do
    row := rows[i];
    if show_polys then
        print "rank", i, "degree", row[1], "genus", row[2], "dir", row[3], "poly", row[4];
    else
        print "rank", i, "degree", row[1], "genus", row[2], "dir", row[3];
    end if;
end for;
