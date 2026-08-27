//////////////////////////////////////////////////////////////////////
//  Bounded Magma decomposition diagnostics for the nonautomorphic
//  [2,6,6] contact-6 slice eps=+1, r=4.
//
//  Typical calls (always wrap expensive modes in an external timeout):
//    magma -b p:=13 mode:="base" this_file
//    magma -b p:=13 mode:="primary" this_file
//    magma -b p:=13 mode:="lex" this_file
//    magma -b p:=0  mode:="base" this_file       // exact Q
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned p then p := 13;
elif Type(p) eq MonStgElt then p := StringToInteger(p); end if;
if not assigned mode then mode := "base"; end if;

K := p eq 0 select Rationals() else GF(p);
R<b,M,U,v> := PolynomialRing(K, 4, "grevlex");
Px<x> := PolynomialRing(R);

r := K!4;
a := 3 - (b+3)*r - 2/r;
h := 1 + a*x + b*x^2 + x^3;
f := h^2 - (x-1)^6;
c := [Coefficient(f,i) : i in [0..5]];
B := c[6]*M + 3*U;
Delta := 4*c[5]*M + 12*(U^2+v^2) - B^2;
F1 := Delta*v^3 - 4*c[2]*M - 12*U*v^4;
F2 := Delta^2 + 64*B*v^3 - 64*c[3]*M
      - 192*(U^2*v^2+v^4);
F3 := B*Delta + 16*v^3 - 8*c[4]*M - 8*U^3 - 48*U*v^2;
q := x^2 + U*x + v^2;

AutoA := [
    2*b^3*r^4 + 30*b^2*r^4 - 36*b^2*r^3 + 8*b^2*r^2
      + 126*b*r^4 - 296*b*r^3 + 216*b*r^2 - 48*b*r
      + 162*r^4 - 612*r^3 + 816*r^2 - 464*r + 96,
    -b^4*r^4 - 12*b^3*r^4 + 12*b^3*r^3 - 2*b^3*r^2
      - 54*b^2*r^4 + 108*b^2*r^3 - 66*b^2*r^2 + 12*b^2*r
      - 108*b*r^4 + 324*b*r^3 - 342*b*r^2 + 152*b*r - 24*b
      - 81*r^4 + 324*r^3 - 470*r^2 + 300*r - 72,
    b^4*r^3 + 6*b^3*r^3 - 6*b^3*r^2 - 14*b^2*r^2
      + 12*b^2*r - 54*b*r^3 + 102*b*r^2 - 48*b*r
      - 81*r^3 + 270*r^2 - 284*r + 96
];
AutoB := [
    12*b^3*r^4 - 12*b^3*r^3 + 4*b^3*r^2 + 108*b^2*r^4
      - 188*b^2*r^3 + 120*b^2*r^2 - 36*b^2*r + 4*b^2
      + 324*b*r^4 - 804*b*r^3 + 732*b*r^2 - 296*b*r + 48*b
      + 324*r^4 - 1044*r^3 + 1224*r^2 - 612*r + 108,
    -6*b^3*r^4 + 6*b^3*r^3 - 2*b^3*r^2 - 54*b^2*r^4
      + 94*b^2*r^3 - 66*b^2*r^2 + 24*b^2*r - 4*b^2
      - 162*b*r^4 + 378*b*r^3 - 342*b*r^2 + 144*b*r - 24*b
      - 162*r^4 + 450*r^3 - 470*r^2 + 216*r - 36,
    6*b^2*r^2 - 6*b^2*r + 2*b^2 + 24*b*r^3 - 24*b*r^2
      + 4*b*r + 72*r^3 - 142*r^2 + 90*r - 18
];

structural := (b+3)*M*v*(U^2-4*v^2)*(a+b+2);
smooth := Discriminant(f);
coprime := Resultant(q,f);

function EvalAt(g, pt)
    return Evaluate(g, <pt[1],pt[2],pt[3],pt[4]>);
end function;

function InIdealAt(I, pt)
    return &and[EvalAt(g,pt) eq 0 : g in Basis(I)];
end function;

procedure Summary(label, I)
    dim, degs := Dimension(I);
    print label, "dimension", dim, "component_degrees", degs;
    try print label, "degree", Degree(Homogenization(I)); catch e print label, "degree_error", e`Object; end try;
    print label, "basis_length", #Basis(I), "basis_summary",
          [<TotalDegree(g),#Terms(g)> : g in Basis(I)];
end procedure;

print "contact6_m36_266_r4_magma_decompose";
print "characteristic", p, "mode", mode;
print "equation_summary", [<TotalDegree(g),#Terms(g)> : g in [F1,F2,F3]];
print "boundary_degrees", TotalDegree(structural), TotalDegree(smooth), TotalDegree(coprime);

Iraw := ideal<R | F1,F2,F3>;
Summary("raw", Iraw);
Istruct := Saturation(Iraw, ideal<R | structural>);
Summary("structural", Istruct);
Ismooth := Saturation(Istruct, ideal<R | smooth>);
Iopen := Saturation(Ismooth, ideal<R | coprime>);
InoA := Saturation(Iopen, ideal<R | AutoA>);
InoAB := Saturation(InoA, ideal<R | AutoB>);
print "saturation_equalities",
      Istruct eq Ismooth, Ismooth eq Iopen, Iopen eq InoA, InoA eq InoAB;
Summary("nonautomorphic_open", InoAB);

if p eq 13 then known := <K!5,K!1,K!4,K!6>;
elif p eq 19 then known := <K!5,K!6,K!9,K!3>;
else known := <K!0,K!0,K!0,K!0>; end if;
if p in {13,19} then
    print "known_point", known,
          "on_raw", InIdealAt(Iraw,known),
          "on_open", InIdealAt(InoAB,known),
          "structural", EvalAt(structural,known),
          "smooth", EvalAt(smooth,known),
          "coprime", EvalAt(coprime,known),
          "autoA", [EvalAt(g,known) : g in AutoA],
          "autoB", [EvalAt(g,known) : g in AutoB];
end if;

if mode eq "radical" then
    Irad := Radical(InoAB);
    Summary("radical", Irad);
    print "already_radical", Irad eq InoAB;
    print "radical_is_prime", IsPrime(Irad);
elif mode eq "primary" then
    comps := PrimaryDecomposition(InoAB);
    print "primary_count", #comps;
    for i in [1..#comps] do
        C := comps[i];
        Summary(Sprintf("primary_%o",i), C);
        print "primary", i, "is_prime", IsPrime(C);
        if p in {13,19} then print "primary", i, "contains_known", InIdealAt(C,known); end if;
    end for;
elif mode in {"lex", "split"} then
    L<bb,MM,UU,vv> := PolynomialRing(K,4,"lex");
    phi := hom<R -> L | bb,MM,UU,vv>;
    IL := ideal<L | [phi(g) : g in Basis(InoAB)]>;
    GL := GroebnerBasis(IL);
    print "lex_basis_length", #GL;
    print "lex_basis_summary", [<TotalDegree(g),#Terms(g)> : g in GL];
    for inds in [[4],[3,4],[2,3,4]] do
        candidates := [g : g in GL |
            &and[Degree(g,j) eq 0 : j in [1..4] | not j in inds]];
        print "lex_elimination_subset", inds, "count", #candidates;
        for g in candidates do
            fac := Factorization(g);
            print " elim", <TotalDegree(g),#Terms(g)>,
                  "factor_degrees", [<TotalDegree(t[1]),t[2]> : t in fac];
            if p in {13,19} then
                print " factor_values_at_known",
                      [Evaluate(t[1], known) : t in fac];
            end if;
        end for;
    end for;
    if mode eq "split" then
        plane_candidates := [g : g in GL |
            Degree(g,1) eq 0 and Degree(g,2) eq 0];
        require #plane_candidates eq 1:
            "expected one elimination polynomial in U,v";
        plane_factors := Factorization(plane_candidates[1]);
        print "split_plane_factor_degrees",
              [<TotalDegree(t[1]),t[2]> : t in plane_factors];
        psi := hom<L -> R | b,M,U,v>;
        A2<u,w> := AffineSpace(K,2);
        for i in [1..#plane_factors] do
            plane_factor := plane_factors[i][1];
            Csplit := InoAB + ideal<R | psi(plane_factor)>;
            dim_split, degs_split := Dimension(Csplit);
            print "split", i, "dimension", dim_split,
                  "component_degrees", degs_split;
            print "split", i, "is_prime", IsPrime(Csplit);
            if p in {13,19} then
                print "split", i, "contains_known", InIdealAt(Csplit,known);
            end if;
            cp := Evaluate(plane_factor, <K!0,K!0,u,w>);
            PC := ProjectiveClosure(Curve(A2,cp));
            try
                print "split", i, "plane_genus", Genus(PC);
            catch e
                print "split", i, "plane_genus_error", e`Object;
            end try;
            try
                print "split", i, "plane_geometric_genus", GeometricGenus(PC);
            catch e
                print "split", i, "plane_geometric_genus_error", e`Object;
            end try;
        end for;
    end if;
end if;

print "DONE";
