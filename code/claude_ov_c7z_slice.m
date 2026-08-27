// Lane 2 (claude_ov_c7z): the contact-7 rational-Weierstrass problem in z-coordinates.
//
//   Q(v) = v^5 + c4 v^4 + (2c4-1) v^3 + (c4+c0-1/2) v^2 + 2 c0 v + c0
// satisfies, identically in (c4,c0):   Q(-1) = -1/2,  Q'(-1) = 3,  Q'(0) = 2 Q(0).
// Hence for the five roots v_i, with  z_i = 1/(1+v_i):
//        sum z_i = 6,     prod z_i = 2,     sum 1/(1-z_i) = 3.
// k rational z_i  <->  quintic factor type with k rational roots  <->  2-rank k-1 (k<5) / 4 (k=5)
//        k=3 -> [2,2,14] (order 56, realized);  k=5 -> [2,2,2,14] (order 112, the target).
//
// SLICE z_1 = c:  the other four w_i satisfy  e1 = 6-c, e4 = 2/c, and P'(1)/P(1) = 3 - 1/(1-c),
// which makes e3 an affine function of e2.  So the slice is a PENCIL of quartics
//        P_c,m(w) = w^4 - e1 w^3 + m w^2 - (A m + B) w + e4 ,
// i.e. m = M_c(w) for a degree-4 rational map M_c, and
//   * "two more rational z" (total 3 -> [2,2,14]) = the correspondence curve M_c(x) = M_c(y);
//   * "four more rational z" (total 5 -> [2,2,2,14]) = a point where P_c,m splits completely.
// This script computes those curves and their genus.
SetColumns(0);
Qx<c> := FunctionField(Rationals());
Pol<w> := PolynomialRing(Qx);

function SliceData(cv)
  // returns e1, e4, A, B for the slice z1 = cv  (cv in a field)
  K := Parent(cv);
  T := 3 - 1/(1-cv);
  e1 := 6 - cv;  e4 := 2/cv;
  // e2*(2-T) + e3*(T-1) = T*(c-5+2/c) - 3c + 14
  A := (T-2)/(T-1);
  B := (T*(cv-5+2/cv) - 3*cv + 14)/(T-1);
  return e1, e4, A, B;
end function;

procedure DoSlice(cv, nm)
  K := Parent(cv);
  e1, e4, A, B := SliceData(cv);
  R<x,y> := PolynomialRing(K,2);
  // m = M(w) = -(w^4 - e1 w^3 - B w + e4)/(w^2 - A w)
  Nm := func<t| -(t^4 - e1*t^3 - B*t + e4)>;
  Dn := func<t| t^2 - A*t>;
  N := Nm(x)*Dn(y) - Nm(y)*Dn(x);
  F := ExactQuotient(N, x-y);
  printf "=== slice %o : e1=%o e4=%o A=%o B=%o\n", nm, e1, e4, A, B;
  printf "  correspondence curve F, bidegree (%o,%o), %o terms\n", Degree(F,x), Degree(F,y), #Terms(F);
  printf "  F = %o\n", F;
  fa := Factorisation(F);
  printf "  #irreducible factors = %o : %o\n", #fa, [<Degree(g[1],x),Degree(g[1],y),g[2]> : g in fa];
  for g in fa do
    C := Curve(AffineSpace(R), g[1]);
    ok, gg := true, -1;
    try gg := Genus(ProjectiveClosure(C)); catch e ok := false; end try;
    if ok then
      printf "    factor bidegree (%o,%o):  GENUS = %o\n", Degree(g[1],x),Degree(g[1],y), gg;
    else
      printf "    factor bidegree (%o,%o):  genus computation failed\n", Degree(g[1],x),Degree(g[1],y);
    end if;
  end for;
end procedure;

// ---- the slice c = 2  (v = -1/2), which carries 5 of the 11 known three-root hits ----
DoSlice(Rationals()!2, "c=2");
// ---- generic slice, over Q(c) ----
DoSlice(c, "generic c");
printf "SLICE_DONE\n";
quit;
