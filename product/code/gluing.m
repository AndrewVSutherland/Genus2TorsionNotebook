/*
  gluing.m -- the elliptic-curve gluing intrinsics used by the split-Jacobian
  workstream in product/code: the subset of A. Sutherland's genus2.m package
  that this repository actually uses, made self-contained.

  Implements the intrinsics

    Genus2Elliptic2  -- (2,2)-gluing of two elliptic curves along their
                        2-torsion, following Bröker-Howe-Lauter-Stevenhagen
                        (MR3349314), based on Howe-Leprevost-Poonen (MR2514865)
    Genus2Elliptic3  -- (3,3)-gluing along 3-torsion, Algorithm 5.4 of [BHLS];
                        the BHLS3 code below is copied from Genus2.magma,
                        https://ewhowe.com/Magma/Genus2.magma
    hyperellred      -- reduced model of a genus 2 curve over Q via PARI/GP

  The small helper functions sprint, CoefficientString, and polredabs replace
  the corresponding utilities of the full package's dependencies (utils.m,
  polredabs.m), so this file has no Magma dependencies; PARI/GP ("gp") must
  be on the PATH for polredabs (used by Genus2Elliptic2) and hyperellred.

  Usage (from product/code):  Attach("gluing.m");
*/

// whitespace-stripped Sprint
function sprint(x)
    return &cat Split(Sprint(x), " \n");
end function;

// "[[coeffs(f)],[coeffs(h)]]" (ascending degree) for y^2 + h(x)y = f(x)
function CoefficientString(C)
    f, h := HyperellipticPolynomials(C);
    cs := func<g | "[" cat Join([sprint(c) : c in Coefficients(g)], ",") cat "]">;
    return "[" cat cs(f) cat "," cat cs(h) cat "]";
end function;

// canonical (polredabs) form of a number field, via PARI/GP
function polredabs(K)
    if Type(K) eq FldRat then return RationalsAsNumberField(); end if;
    f := DefiningPolynomial(K);
    if Degree(f) le 1 then return RationalsAsNumberField(); end if;
    cmd := Sprintf("{print(Vecrev(Vec(polredabs(Pol(%o)))))}",
                   sprint(Reverse(Coefficients(f))));
    s := Pipe("gp -q", cmd);
    R := PolynomialRing(Rationals());
    return NumberField(R!(eval sprint(s)));
end function;

intrinsic Genus2Elliptic2(f::RngUPolElt,g::RngUPolElt) -> SeqEnum[CrvHyp]
{ Given non-isomorphic elliptic curves with isomorphic 2-torsion Galois modules, returns a list of hyperelliptic curves C that is a degree-2 cover of both E1 and E2 using the algorithm given in Reinier-Howe-Lauter-Stevenhagen (MR3349314), which is based on Howe-Leprevost-Poonen (MR2514865). }
    Kf := polredabs(SplittingField(f)); Kg := polredabs(SplittingField(g));
    require DefiningPolynomial(Kf) eq DefiningPolynomial(Kg): Sprintf("Elliptic curves must have the same 2-torsion field, %o ne %o",Kf,Kg);
    Df := Discriminant(f); Dg := Discriminant(g);
    alpha := [r[1]:r in Roots(ChangeRing(f,Kf))]; gamma := [r[1]:r in Roots(ChangeRing(g,Kf))];
    L := [];
    G := Automorphisms(Kf);
    Gaction := func<rs|[[Index(rs,g(r)):r in rs]:g in G]>;
    a:=[[alpha[i]-alpha[j]:j in [1..3]]:i in [1..3]];
    R<t>:=PolynomialRing(Kf);
    F := BaseRing(f);
    for sigma in SymmetricGroup(3) do
        beta := [gamma[i^sigma]:i in [1..3]];
        b:=[[beta[i]-beta[j]:j in [1..3]]:i in [1..3]];
        if alpha[1]*b[3][2] + alpha[2]*b[1][3] + alpha[3]*b[2][1] ne 0 and Gaction(alpha) eq Gaction(beta) then
            a1 := a[3][2]^2/b[3][2]+a[2][1]^2/b[2][1]+a[1][3]^2/b[1][3]; a2:=alpha[1]*b[3][2]+alpha[2]*b[1][3]+alpha[3]*b[2][1];
            b1 := b[3][2]^2/a[3][2]+b[2][1]^2/a[2][1]+b[1][3]^2/a[1][3]; b2:=beta[1]*a[3][2]+beta[2]*a[1][3]+beta[3]*a[2][1];
            A := Dg*a1/a2;  B:=Df*b1/b2;
            Append(~L,ChangeRing(-(A*a[2][1]*a[1][3]*t^2+B*b[2][1]*b[1][3])*(A*a[3][2]*a[2][1]*t^2+B*b[3][2]*b[2][1])*(A*a[1][3]*a[3][2]*t^2+B*b[1][3]*b[3][2]),F));
        end if;
    end for;
    return [ReducedMinimalWeierstrassModel(HyperellipticCurve(f)):f in L];
end intrinsic;

intrinsic Genus2Elliptic2(E1::CrvEll,E2::CrvEll) -> SeqEnum[CrvHyp]
{ Given non-isomorphic elliptic curves with isomorphic 2-torsion Galois modules, returns a list of hyperelliptic curves C that is a degree-2 cover of both E1 and E2 using the algorithm given in Reinier-Howe-Lauter-Stevenhagen (MR3349314), which is based on Howe-Leprevost-Poonen (MR2514865). }
    f := HyperellipticPolynomials(WeierstrassModel(E1));
    g := HyperellipticPolynomials(WeierstrassModel(E2));
    return Genus2Elliptic2(f,g);
end intrinsic;

intrinsic Genus2Elliptic2(E1::MonStgElt,E2::MonStgElt) -> SeqEnum[CrvHyp]
{ Given non-isomorphic elliptic curves with isomorphic 2-torsion Galois modules, returns a list of hyperelliptic curves C that is a degree-2 cover of both E1 and E2 using the algorithm given in Reinier-Howe-Lauter-Stevenhagen (MR3349314), which is based on Howe-Leprevost-Poonen (MR2514865). }
    return Genus2Elliptic2(EllipticCurve(E1),EllipticCurve(E2));
end intrinsic;


// The code below is copied from Genus2.magma downloaded from https://ewhowe.com/Magma/Genus2.magma

function normalize_quintuple(quint)
  // Given a quintuple (a,b,c,d,t) as in the Appendix to [BHLS], return
  // a semi-normal form of the quintuple under the action of k* x k*
  // defined in the Appendix.
  a := quint[1]; b := quint[2]; c := quint[3]; d := quint[4]; t := quint[5];
  if (a eq 0 and c eq 0) or (b eq 0 and d eq 0) then return quint; end if;
  x := a*b ne 0 select a/b else d/c;
  return [x^2*a, x^3*b, c/x^2, d/x^3, x*t];
end function;


function are_equivalent(quint1, quint2);
  // Are two quintuples equivalent under the action of k* x k* defined
  // in the Appendix to [BHLS]?
  q1 := normalize_quintuple(quint1);
  q2 := normalize_quintuple(quint2);
  if [z eq 0 : z in q1] ne [z eq 0 : z in q2] then return false; end if;
  if #[z : z in q1 | z eq 0] lt 2 then
    return (&and[q1[i] eq q2[i] : i in [1..4]]) and IsSquare(q1[5]*q2[5]);
  end if;
  K := Parent(q1[1]);
  y := PolynomialRing(K).1;
  if q1[1] eq 0 and q1[3] eq 0 then
    rts := Roots(y^3 - q2[2]/q1[2]);
    if #rts eq 0 then return false; end if;
    rts := [r[1] : r in rts];
    return &or[q2[2] eq x^3*q1[2] and q1[4] eq x^3*q2[4] and IsSquare(x*q1[5]*q2[5]) : x in rts];
  end if;

  // We must have q1[2] = q1[4] = 0.
  rts := Roots(y^2 - q2[1]/q1[1]);
  if #rts eq 0 then return false; end if;
  rts := [r[1] : r in rts];
  return &or[q2[1] eq x^2*q1[1] and q1[3] eq x^2*q2[3] and IsSquare(x*q1[5]*q2[5]) : x in rts];
end function;

function BHLS3(f,g)
  // Given two monic separable cubics in F_q[x], where q is a prime
  // power not divisible by 2 or 3, produce all of the genus-2 curves
  // that can be obtained by gluing y^2 = f and y^2 = g together along
  // their 3-torsion subgroups.  We use Algorithm 5.4 of [BHLS].

  abcdlist := [];
  abcdtlist := [];
  x := Parent(f).1;
  K := BaseRing(Parent(f));
  R := PolynomialRing(K,2);
  E1 := EllipticCurve(f);  j1 := jInvariant(E1);
  E2 := EllipticCurve(g);  j2 := jInvariant(E2);

  // Want solutions [a:b:c:d] that solve 12*w*y + 16*x*z = 1 (among
  // other equations).  We first assume that a*b != 0 and then scale
  // so that a = b.  Then 12*a*c + 16*a*d = 1, so we may as well set
  // a = 1/(12*c + 16*d).

  c := R.1;  d := R.2;  b := 1/(12*c + 16*d); a := b;
  g1 := 1728*(a^2*c   + 4*a*b*d - 4*b^2*c^2)^3 - j1*(a^3 + b^2)^2 * (c^3 + d^2)  ;
  g2 := 1728*(a  *c^2 + 4*b*c*d - 4*a^2*d^2)^3 - j2*(a^3 + b^2)   * (c^3 + d^2)^2;
  f1 := Numerator(g1);  f2 := Numerator(g2);
  g := UnivariatePolynomial(Resultant(f1,f2,1));
  rg := [r[1] : r in Roots(g) | not 8*r[1] eq 1];
  for dd in rg do
    cpol := GCD(Evaluate(f1,[x,dd]), Evaluate(f2,[x,dd]));
    for cn in Roots(cpol) do
      cc := cn[1];
      if 12*cc + 16*dd ne 0 then
        bb := 1/(12*cc + 16*dd);
        abcdlist cat:= [ [bb,bb,cc,dd] ];
      end if;
    end for;
  end for;

  // Now, the solutions with a*b = 0 but c*d != 0.
  a := R.1;  b := R.2;  c := 1/(12*a + 16*b);  d := c;
  g1 := 1728*(a^2*c   + 4*a*b*d - 4*b^2*c^2)^3 - j1*(a^3 + b^2)^2 * (c^3 + d^2)  ;
  g2 := 1728*(a  *c^2 + 4*b*c*d - 4*a^2*d^2)^3 - j2*(a^3 + b^2)   * (c^3 + d^2)^2;
  f1 := Numerator(g1);  f2 := Numerator(g2);
  g := UnivariatePolynomial(Resultant(f1,f2,1));
  rg := [r[1] : r in Roots(g) | not 8*r[1] eq 1];
  for bb in rg do
    apol := GCD(Evaluate(f1,[x,bb]), Evaluate(f2,[x,bb]));
    for an in Roots(apol) do
      aa := an[1];
      if 12*aa + 16*bb ne 0 and aa*bb eq 0 then
        cc := 1/(12*aa + 16*bb);
        abcdlist cat:= [ [aa,bb,cc,cc] ];
      end if;
    end for;
  end for;

  for abcd in abcdlist do
    a := abcd[1]; b := abcd[2]; c := abcd[3]; d := abcd[4];
    Delta1 := a^3 + b^2;  Delta2 := c^3 + d^2;
    f1 := x^3 + 12*(2*a^2*d - b*c)*x^2 + 12*(16*a*d^2 + 3*c^2)*Delta1*x + 512*Delta1^2*d^3;
    f2 := x^3 + 12*(2*b*c^2 - a*d)*x^2 + 12*(16*b^2*c + 3*a^2)*Delta2*x + 512*Delta2^2*b^3;
    F1 := EllipticCurve(f1); F2 := EllipticCurve(f2);
    assert jInvariant(F1) eq j1 and jInvariant(F2) eq j2;

    // Is there a simultaneous twist of F1 and F2 that will make them
    // isomorphic to E1 and E2?  There is a subtlety to answering this,
    // because if -1 is not a square, there may be *two* classes in
    // k^* mod squares that produce the same twist of the j=1728 curve.
    bool1, t1 := IsQuadraticTwist(E1,F1);
    bool2, t2 := IsQuadraticTwist(E2,F2);
    if bool1 and bool2 then
      if IsSquare(t1*t2) then
        abcdtlist cat:= [[a,b,c,d,t1]];
      end if;
      if (j1 eq 1728 or j2 eq 1728) and not IsSquare(K!-1) then
        if j1 eq 1728 then t := t2; else t := t1; end if;
        if IsSquare(-t1*t2) then
          abcdtlist cat:= [[a,b,c,d,t]];
        end if;
      end if;
    end if;

  end for;

  // Now the exceptional cases.

  if j1 eq j2 then
    eq1 := HyperellipticPolynomials(WeierstrassModel(E1));
    eq2 := HyperellipticPolynomials(WeierstrassModel(E2));
    if j1 eq 0 then
      e1 := Coefficient(eq1,0);
      e2 := Coefficient(eq2,0);
      if #Roots(x^6 - e1*e2/4) gt 0 then
        b := e1;  d := 1/(16*e1);  a := 0;  c := 0;  t := 2;
        abcdtlist cat:= [[a,b,c,d,t]];
      end if;
    end if;

    if j1 eq 1728 then
      e1 := Coefficient(eq1,1);
      e2 := Coefficient(eq2,1);
      if #Roots(x^4 - e1*e2/108) gt 0 then
        a := e1;  c := 1/(12*e1);  b := 0;  d := 0;  t := 2;
        abcdtlist cat:= [[a,b,c,d,t]];
      end if;
    end if;
  end if;

  if #abcdtlist eq 0 then return [Parent(f)|]; end if;

  shortlist := [abcdtlist[1]];

  for q in abcdtlist do
    if not &or[are_equivalent(q,qq) : qq in shortlist] then
      shortlist cat:= [q];
    end if;
  end for;

  hlist := [];
  for abcdt in shortlist do
    a := abcdt[1];  b := abcdt[2];  c := abcdt[3];  d := abcdt[4];  t := abcdt[5];
    hlist cat:= [t * (x^3 + 3*a*x + 2*b) * (2*d*x^3 + 3*c*x^2 + 1)];
  end for;
  return hlist;
end function;

intrinsic Genus2Elliptic3(f::RngUPolElt,g::RngUPolElt) -> SeqEnum[CrvHyp]
{ Given non-isomorphic elliptic curves with isomorphic 2-torsion Galois modules, returns a list of hyperelliptic curves C that is a degree-2 cover of both E1 and E2 using the algorithm given in Reinier-Howe-Lauter-Stevenhagen (MR3349314). }
    hs := BHLS3(f,g);
    return [ReducedMinimalWeierstrassModel(HyperellipticCurve(h)) : h in hs];
end intrinsic;

intrinsic Genus2Elliptic3(E1::CrvEll,E2::CrvEll) -> SeqEnum[CrvHyp]
{ Given non-isomorphic elliptic curves with isomorphic 2-torsion Galois modules, returns a list of hyperelliptic curves C that is a degree-2 cover of both E1 and E2 using the algorithm given in Reinier-Howe-Lauter-Stevenhagen (MR3349314), which is based on Howe-Leprevost-Poonen (MR2514865). }
    f := HyperellipticPolynomials(WeierstrassModel(E1));
    g := HyperellipticPolynomials(WeierstrassModel(E2));
    return Genus2Elliptic3(f,g);
end intrinsic;

intrinsic Genus2Elliptic3(E1::MonStgElt,E2::MonStgElt) -> SeqEnum[CrvHyp]
{ Given non-isomorphic elliptic curves with isomorphic 2-torsion Galois modules, returns a list of hyperelliptic curves C that is a degree-2 cover of both E1 and E2 using the algorithm given in Reinier-Howe-Lauter-Stevenhagen (MR3349314), which is based on Howe-Leprevost-Poonen (MR2514865). }
    return Genus2Elliptic3(EllipticCurve(E1),EllipticCurve(E2));
end intrinsic;

intrinsic hyperellred(C::CrvHyp) -> CrvHyp
{ returns a reduced model of the genus 2 curve y^2+h(x)y=f(x) specified by [coeffs(f),coeffs(h)] using PARI/GP hyeprellminimalmodel. }
    require BaseRing(C) eq Rationals(): "Curve must be defined over Q";
    cmd := Sprintf("{print(hyperellred(apply(Polrev,%o)))}", CoefficientString(C));
    R<x> := PolynomialRing(Rationals());
    fh := eval(sprint(Pipe("gp -q 2>/dev/null", cmd)));
    return HyperellipticCurve(fh[1],fh[2]);
end intrinsic;
