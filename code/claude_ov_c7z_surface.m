// ===========================================================================
// Lane 2 (claude_ov_c7z) -- the FULL order-112 locus of the contact-7 chart,
// beyond the decided c=2 slice.
//
// PART A.  Universal relations (re-derived in claude_ov_c7z_xdecide.m):
//    S = { (z1..z5) : e1 = 6, e5 = 2, e3 = 2 e4 - 2 }  subset A^5.
//    Its rational points ARE exactly the order-112 configurations (all z_i in Q).
//    Projective closure: a complete intersection of a hyperplane and two
//    hypersurfaces of degrees 4 and 5 in P^5, i.e. a (4,5) complete-intersection
//    SURFACE in P^4, so K_S = O_S(4+5-5) = O_S(4)  ->  general type.
//
// PART B.  Slice z_1 = c.  Cab_c is the locus where the residual quartic splits
//    into two rational quadratics; the three-root curve is the double cover
//    u^2 = d1, the order-112 curve X_c the (Z/2)^2-cover u^2=d1, v^2=d2.
//    We compute genus(Cab_c), genus(three-root curve) and genus(X_c) for many c.
// ===========================================================================
SetColumns(0);
if not assigned MemGB then MemGB := 16; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q := Rationals();

printf "======== PART A: the order-112 surface S ========\n";
P4<z1,z2,z3,z4,tt> := ProjectiveSpace(Q,4);
z5 := 6*tt - z1 - z2 - z3 - z4;
zz := [z1,z2,z3,z4,z5];
function Esym(v,k)
  R := Parent(v[1]);
  s := R!0;
  for T in Subsets({1..#v}, k) do s +:= &*[ v[i] : i in T ]; end for;
  return s;
end function;
E3 := Esym(zz,3); E4 := Esym(zz,4); E5 := Esym(zz,5);
Eq1 := tt*E3 - 2*E4 + 2*tt^4;      // e3 = 2 e4 - 2   (homogeneous of degree 4)
Eq2 := E5 - 2*tt^5;                // e5 = 2          (homogeneous of degree 5)
printf "Eq1 degree %o, Eq2 degree %o\n", TotalDegree(Eq1), TotalDegree(Eq2);
S := Scheme(P4, [Eq1, Eq2]);
printf "dim S = %o ; degree = %o\n", Dimension(S), Degree(S);
printf "S irreducible = %o ; reduced = %o\n", IsIrreducible(S), IsReduced(S);
printf "Numerical invariants of a SMOOTH (4,5) c.i. surface in P^4 (for reference):\n";
printf "   omega = O(4), K^2 = 4*5*(4+5-5)^2 = %o,  p_g = h^0(P^4,O(4)) - 1 = %o,  q = 0\n",
   4*5*(4+5-5)^2, Binomial(4+4,4)-1;
printf "   (S is a complete intersection, hence Gorenstein, so omega_S = O_S(4) regardless.)\n";
tt0 := Cputime();
try
  Sing := SingularSubscheme(S);
  printf "dim(Sing S) = %o  [%o s]\n", Dimension(Sing), Cputime(tt0);
  Hinf := Scheme(P4, tt);
  SingInf := Sing meet Hinf;
  printf "dim(Sing S meet {t=0}) = %o\n", Dimension(SingInf);
  Saff := Difference(Sing, Hinf);
  printf "Sing S is contained in the hyperplane at infinity {t=0}: %o\n", IsEmpty(Saff);
  comps := IrreducibleComponents(Sing);
  printf "Sing S has %o irreducible components, dims %o, degrees %o\n",
     #comps, [Dimension(c) : c in comps], [Degree(ReducedSubscheme(c)) : c in comps];
catch e
  printf "SingularSubscheme analysis failed: %o [%o s]\n", e`Object, Cputime(tt0);
end try;
// the hyperplane section at infinity
Binf := Scheme(P4, [Eq1, Eq2, tt]);
printf "boundary S meet {t=0}: dim %o, degree(reduced) %o, #components %o\n",
   Dimension(Binf), Degree(ReducedSubscheme(Binf)), #IrreducibleComponents(Binf);
// the affine part that actually carries configurations
A4<y1,y2,y3,y4> := AffineSpace(Q,4);
y5 := 6 - y1-y2-y3-y4;
yy := [y1,y2,y3,y4,y5];
F1 := Esym(yy,3) - 2*Esym(yy,4) + 2;
F2 := Esym(yy,5) - 2;
Saf := Scheme(A4,[F1,F2]);
printf "affine S: dim %o, irreducible %o, reduced %o\n", Dimension(Saf), IsIrreducible(Saf), IsReduced(Saf);
try
  SgA := SingularSubscheme(Saf);
  printf "affine S: IsEmpty(Sing) = %o", IsEmpty(SgA);
  if not IsEmpty(SgA) then printf " ; dim = %o, degree = %o", Dimension(SgA), Degree(ReducedSubscheme(SgA)); end if;
  printf "\n";
catch e printf "affine SingularSubscheme failed: %o\n", e`Object; end try;

printf "\n======== PART B: slice genera ========\n";
R2<aa,bb> := PolynomialRing(Q,2);
A2 := AffineSpace(R2);
Pt1<XX> := PolynomialRing(Q);

function SliceCurve(c)
  // returns the defining polynomial of Cab_c and the two discriminants d1,d2
  // c = 0, 1 and 1/2 are the degenerate z-values (T3 = 1 for c = 1/2)
  if c eq 0 or c eq 1 then return R2!0, 0, 0; end if;
  T3 := 3 - 1/(1-c);
  if T3 eq 1 then return R2!0, 0, 0; end if;
  e1 := 6 - c; e4 := 2/c;
  A := (T3-2)/(T3-1);
  B := (T3*(c-5+2/c) - 3*c + 14)/(T3-1);
  F := (e1 - A - aa)*bb^2 + (A*aa^2 - A*e1*aa - B)*bb + e4*(aa - A);
  return F, e1, e4;
end function;

// sanity: c = 2 must reproduce  -3ab^2+10b^2+2a^2b-8ab+3a-2  (up to a constant)
Fc2 := SliceCurve(2);
printf "c=2 slice poly * 3 = %o\n", 3*Fc2;
assert 3*Fc2 eq -3*aa*bb^2 + 10*bb^2 + 2*aa^2*bb - 8*aa*bb + 3*aa - 2;
printf "c=2 cross-check against claude_ov_c7z_xdecide.m: OK\n";

procedure DoC(c)
  F, e1, e4 := SliceCurve(c);
  if F eq 0 then printf "c=%-8o : degenerate (F identically 0)\n", c; return; end if;
  C := Curve(A2, F);
  irr := IsIrreducible(F);
  g0 := -1;
  try g0 := Genus(ProjectiveClosure(C)); catch e ; end try;
  g1 := -1; g2 := -1;
  if irr and g0 ge 0 then
    try
      KC := FunctionField(C);
      FF, mp := AlgorithmicFunctionField(KC);
      aF := mp(KC!(A2.1)); bF := mp(KC!(A2.2));
      d1F := aF^2 - 4*bF;
      d2F := (e1-aF)^2 - 4*e4/bF;
      PF<uu> := PolynomialRing(FF);
      L1 := FunctionField(uu^2 - d1F);
      g1 := Genus(L1);
      PL<vv> := PolynomialRing(L1);
      L2 := FunctionField(vv^2 - (L1!d2F));
      g2 := Genus(L2);
    catch e
      printf "   (cover genus failed for c=%o: %o)\n", c, e`Object;
    end try;
  end if;
  printf "c=%-10o : Cab irred=%-5o genus(Cab)=%-3o genus(3-root cover)=%-3o genus(X_c)=%o\n",
     c, irr, g0, g1, g2;
end procedure;

// the special slice and a spread of others (c = 0, 1 are excluded: e4 or A undefined)
clist := [Q| 2, -1, 3, 4, 5, -2, 3/2, -1/3, 5/2, 7, 2/3, -5, 10, 6, -1/2, -1/9, -7/3, 22/7, 49/36 ];
for c in clist do DoC(c); end for;
printf "SURFACE_DONE\n";
quit;
