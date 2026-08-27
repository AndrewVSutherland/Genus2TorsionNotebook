//////////////////////////////////////////////////////////////////////
// Exact monodromy/genus certificate for the marked auxiliary covers on
//
//   F_t = F_HLP + t*(2+x-x^2+x^3+x^4+x^5+x^6).
//
// The support covers are finite etale over the smooth-fiber locus:
//   * nonzero J[3] modulo +/- has degree (3^4-1)/2 = 40;
//   * primitive J[4] modulo +/- has degree (4^4-2^4)/2 = 120.
//
// A simple nodal fiber acts by one symplectic transvection.  This file
// checks that the sextic pencil has ten simple finite nodal fibers, no
// degeneration at infinity, and connected root monodromy.  It then
// enumerates the exact transvection cycle types on the two finite sets
// and applies Riemann--Hurwitz.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Q := Rationals();
Qt<t> := PolynomialRing(Q);
P<x> := PolynomialRing(Qt);

F0 := 187392*x^6-118767*x^4-118767*x^2+187392;
GA := 2+x-x^2+x^3+x^4+x^5+x^6;
F := F0+t*GA;

disc := Discriminant(F);
disc_sf := SquarefreePart(disc);
print "M612_HLP_GA_MONODROMY_GEOMETRY";
print "discriminant_degree",Degree(disc);
print "discriminant_squarefree",Degree(disc_sf) eq Degree(disc);
print "discriminant_factor_degrees_Q",
      [<Degree(a[1]),a[2]> : a in Factorization(disc)];
print "direction_sextic_discriminant_nonzero",Discriminant(GA) ne 0;

K<T> := FunctionField(Q);
PK<X> := PolynomialRing(K);
FK := PK![ K!Coefficient(F,i) : i in [0..Degree(F)] ];
facF := Factorization(FK);
print "sextic_factor_degrees_Qt",[<Degree(a[1]),a[2]> : a in facF];
assert #facF eq 1 and Degree(facF[1][1]) eq 6 and facF[1][2] eq 1;
assert Degree(disc) eq 10 and Degree(disc_sf) eq 10;
assert Discriminant(GA) ne 0;

function CanonicalSign(v,n)
    nv := [(-a) mod n : a in v];
    return v lt nv select v else nv;
end function;

function PrimitiveSignClasses(n)
    tuples := CartesianPower({0..n-1},4);
    out := {};
    for tup in tuples do
        v := [Integers()!a : a in tup];
        primitive := n eq 3 select (&or[a ne 0 : a in v])
                                  else (&or[IsOdd(a) : a in v]);
        if primitive then Include(~out,CanonicalSign(v,n)); end if;
    end for;
    return out;
end function;

function Transvection(v,n)
    // Symplectic coordinates (e1,e2,f1,f2), vanishing vector e1.
    c := (-v[3]) mod n;
    w := [(v[1]+c) mod n,v[2],v[3],v[4]];
    return CanonicalSign(w,n);
end function;

function CycleData(n)
    S := PrimitiveSignClasses(n);
    seen := {};
    lengths := [];
    for v in S do
        if v in seen then continue; end if;
        w := v;
        ell := 0;
        repeat
            Include(~seen,w);
            ell +:= 1;
            w := Transvection(w,n);
        until w in seen;
        Append(~lengths,ell);
    end for;
    Sort(~lengths);
    return #S,lengths,#S-#lengths;
end function;

d3,cyc3,ram3 := CycleData(3);
d4,cyc4,ram4 := CycleData(4);
print "J3_mod_sign_degree",d3;
print "J3_transvection_cycle_counts",
      [<ell,#[a : a in cyc3 | a eq ell]> : ell in Seqset(cyc3)];
print "J3_ramification_per_node",ram3;
print "J4_primitive_mod_sign_degree",d4;
print "J4_transvection_cycle_counts",
      [<ell,#[a : a in cyc4 | a eq ell]> : ell in Seqset(cyc4)];
print "J4_ramification_per_node",ram4;

nodes := Degree(disc);
g3 := 1+(-2*d3+nodes*ram3) div 2;
g4 := 1+(-2*d4+nodes*ram4) div 2;
print "contact_cover_genus",g3;
print "halving_cover_genus",g4;
quit;
