//////////////////////////////////////////////////////////////////////
// Exact function-field geometry of the C3-root cover in the
// simultaneous contact-5/contact-6 order-30 family.
//
// Default run (including the three descended Richelot-halving covers):
//   magma -b code/contact30_c3root_trigonal_geometry.m
//
// The two halving covers on the original genus-12 curve can be recomputed
// with ComputeOriginalHalvingCovers:=true.  They are omitted by default
// because the quotient computations are both stronger and much faster.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(6*10^9);

if not assigned ComputeHalvingCovers then
    ComputeHalvingCovers := true;
elif Type(ComputeHalvingCovers) eq MonStgElt then
    ComputeHalvingCovers := StringToLower(ComputeHalvingCovers) in
        { "true", "t", "yes", "y", "1" };
end if;
if not assigned ComputeOriginalHalvingCovers then
    ComputeOriginalHalvingCovers := false;
elif Type(ComputeOriginalHalvingCovers) eq MonStgElt then
    ComputeOriginalHalvingCovers := StringToLower(
        ComputeOriginalHalvingCovers) in { "true", "t", "yes", "y", "1" };
end if;

Q := Rationals();
K<R> := FunctionField(Q);
P<rho> := PolynomialRing(K);

function RootData(eps)
    t := (5*R^2-20*R+19)/(R^2-5);
    Y := -2*(5*R^2-22*R+25)/(R^2-5);
    u := t^3;
    s := t^5+t^4+(K!5/2)*t^3+(K!1/2)*t
       + eps*t*(t-K!1/2)*(t+1)*Y;
    C := (u^2+1)/(2*u);
    den := u^6+6*u^4*s-2*u^4+15*u^3*s-u*s^3+u^2;
    num := 15*u^5+90*u^4+20*u^3*s-6*u^2*s^2+231*u^3
       +2*u^2*s-15*u*s^2+90*u^2-20*u*s+15*u-2*s;
    q := num/den;
    A := (s+q)/2;
    B := (15-s*q)/2;
    f := 2*rho^3+(A-3)*rho^2+(B+3)*rho+(C-1);
    return f,A,B,C;
end function;

procedure PrintDifferent(D)
    print "different_degree", Degree(D);
    for pl in Support(D) do
        print "  branch_place", Minimum(pl),
              "degree", Degree(pl),
              "e", RamificationIndex(pl),
              "f", InertiaDegree(pl),
              "different_exponent", Valuation(D,pl);
    end for;
end procedure;

psi := (5*R-7)/(3*R-5);
fminus,Aminus,Bminus,Cminus := RootData(-1);
fplus,Aplus,Bplus,Cplus := RootData(1);

print "ORIGINAL TRIGONAL ROOT COVER";
for eps in [-1,1] do
    f,A,B,C := RootData(eps);
    assert IsIrreducible(f);
    assert &and[ Evaluate(Coefficient(f,i),psi) eq Coefficient(f,i)
                 : i in [0..3] ];
    disc := Discriminant(f);
    print "branch",eps,"irreducible",true;
    print "discriminant_numerator_factorization",
          Factorization(Numerator(disc));
    print "discriminant_denominator_factorization",
          Factorization(Denominator(disc));
    F<a> := FunctionField(f);
    assert Genus(F) eq 12;
    D := DifferentDivisor(F);
    assert Degree(D) eq 28;
    print "normalized_genus",Genus(F);
    PrintDifferent(D);
end for;

// The two signs parametrize the same curve.  This base transformation
// fixes t and reverses Y, hence changes eps.
phi := (11*R-25)/(5*R-11);
assert &and[ Evaluate(Coefficient(fplus,i),phi) eq Coefficient(fminus,i)
             : i in [0..3] ];
print "branch_swap_R",phi;
print "root_cover_involution_R",psi,"rho_fixed",true;

//////////////////////////////////////////////////////////////////////
// Quotient by psi.  Its invariant is
//
//     z=(3R^2-7)/(3R-5),
//
// and R satisfies R^2-zR+(5z-7)/3=0.
//////////////////////////////////////////////////////////////////////

L<z> := FunctionField(Q);
PT<T> := PolynomialRing(L);
M<r> := ext<L | T^2-z*T+(5*z-7)/3>;

function QuotientData(eps)
    RR := r;
    t := (5*RR^2-20*RR+19)/(RR^2-5);
    Y := -2*(5*RR^2-22*RR+25)/(RR^2-5);
    u := t^3;
    s := t^5+t^4+(M!5/2)*t^3+(M!1/2)*t
       + eps*t*(t-M!1/2)*(t+1)*Y;
    C := (u^2+1)/(2*u);
    den := u^6+6*u^4*s-2*u^4+15*u^3*s-u*s^3+u^2;
    num := 15*u^5+90*u^4+20*u^3*s-6*u^2*s^2+231*u^3
       +2*u^2*s-15*u*s^2+90*u^2-20*u*s+15*u-2*s;
    q := num/den;
    A := (s+q)/2;
    B := (15-s*q)/2;
    // Exact invariance check in the quadratic basis [1,r].
    assert Eltseq(A)[2] eq 0 and Eltseq(B)[2] eq 0 and Eltseq(C)[2] eq 0;
    aa := L!Eltseq(A)[1];
    bb := L!Eltseq(B)[1];
    cc := L!Eltseq(C)[1];
    PL<x> := PolynomialRing(L);
    return 2*x^3+(aa-3)*x^2+(bb+3)*x+(cc-1),aa,bb,cc;
end function;

fm,Am,Bm,Cm := QuotientData(-1);
fp,Ap,Bp,Cp := QuotientData(1);

// The sign-changing map phi induces this map on z.
zswap := ((L!47/10)*z-22)/(z-L!47/10);
assert &and[ Evaluate(Coefficient(fp,i),zswap) eq Coefficient(fm,i)
             : i in [0..3] ];
print "QUOTIENT TRIGONAL COVER";
print "invariant_z", "(3*R^2-7)/(3*R-5)";
print "branch_swap_z",zswap;

for eps in [-1,1] do
    f,A,B,C := QuotientData(eps);
    assert IsIrreducible(f);
    disc := Discriminant(f);
    print "quotient_branch",eps;
    print "discriminant_numerator_factorization",
          Factorization(Numerator(disc));
    print "discriminant_denominator_factorization",
          Factorization(Denominator(disc));
    F<a> := FunctionField(f);
    assert Genus(F) eq 6;
    D := DifferentDivisor(F);
    assert Degree(D) eq 16;
    print "normalized_genus",Genus(F);
    PrintDifferent(D);
end for;

// For eps=-1, a Q-automorphism of the genus-6 quotient must preserve its
// unique trigonal pencil.  It fixes the unique rational e=3 branch z=14/3
// and either fixes or swaps the rational e=2 branches z=2,5.  The only
// nonidentity base candidate is the following involution.  It fails to
// preserve the unique degree-2 branch place, so it cannot lift.
candidate := (106*z-532)/(21*z-106);
assert Evaluate(candidate,L!14/3) eq L!14/3;
assert Evaluate(candidate,L!2) eq L!5;
assert Evaluate(candidate,L!5) eq L!2;
facqm := Factorization(Numerator(Discriminant(fm)));
quad_places := [ fe[1] : fe in facqm | Degree(fe[1]) eq 2
                 and IsOdd(fe[2]) ];
assert #quad_places eq 1;
q2 := quad_places[1];
ratio := Evaluate(q2,candidate)*Denominator(candidate)^2/q2;
isconstant,constant := IsCoercible(Q,ratio);
assert not isconstant;
print "only_nonidentity_Q_base_candidate",candidate;
print "candidate_preserves_degree2_branch",false;
print "Aut_Q_genus12_root_cover", "C2 generated by R -> psi";
print "Aut_Q_genus6_quotient", "trivial";

//////////////////////////////////////////////////////////////////////
// Richelot-halving squareclasses.  With A0=Q2 and
// B0=C3/(x-rho), their discriminants on the root cover are
//
// DA=(B-3)^2-4(A+3)(C+1),
// DB=(A-3)^2-4(A-3)rho-12rho^2-8(B+3).
//
// Both are psi-invariant and descend to the genus-6 quotient.
//////////////////////////////////////////////////////////////////////

if ComputeHalvingCovers then
    Fq<a> := FunctionField(fm);
    assert Genus(Fq) eq 6;
    DA := Fq!((Bm-3)^2-4*(Am+3)*(Cm+1));
    DB := (Fq!(Am-3))^2-4*(Fq!(Am-3))*a-12*a^2-8*(Fq!(Bm+3));
    PFq<w> := PolynomialRing(Fq);
    names := [ "DA", "DB", "DA*DB" ];
    vals := [ DA, DB, DA*DB ];
    expected_genera := [ 17,18,24 ];
    expected_ramification := [ 12,14,26 ];
    print "DESCENDED RICHELOT-HALVING DOUBLE COVERS";
    for i in [1..3] do
        G<v> := FunctionField(w^2-vals[i]);
        g := Genus(G);
        relram := 2*g-2-2*(2*Genus(Fq)-2);
        assert g eq expected_genera[i];
        assert relram eq expected_ramification[i];
        print names[i],"genus",g,
              "relative_ramification_degree",relram;
    end for;
    // The third Richelot gate W requires DA and DB to be squares
    // simultaneously, not merely DA*DB to be a square.  Its V4 cover has
    // the three double covers above as quotients.  The standard V4 genus
    // relation is g(V4)=g(DA)+g(DB)+g(DA*DB)-2*g(base).
    v4genus := &+expected_genera - 2*Genus(Fq);
    assert v4genus eq 47;
    print "W_BOTH_DA_AND_DB_V4", "genus",v4genus;
end if;

if ComputeOriginalHalvingCovers then
    Fo<a> := FunctionField(fminus);
    DA := Fo!((Bminus-3)^2-4*(Aminus+3)*(Cminus+1));
    DB := (Fo!(Aminus-3))^2-4*(Fo!(Aminus-3))*a-12*a^2
        -8*(Fo!(Bminus+3));
    PFo<w> := PolynomialRing(Fo);
    print "ORIGINAL GENUS-12 HALVING COVERS";
    for pair in [ <"DA",DA,35>, <"DB",DB,36> ] do
        G<v> := FunctionField(w^2-pair[2]);
        assert Genus(G) eq pair[3];
        print pair[1],"genus",Genus(G);
    end for;
end if;

quit;
