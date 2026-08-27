# R = -8: A-side norm quotients on the lambda-line (Chabauty targets)

Date: 2026-07-02.  Follows agent_m18_416_other_ELS_v4_scan.md (R=-8 the
last open ELS fiber; S_B-side V_4 quotients all rank 1).

## Scripts

```text
code/agent_m18_416_R8_SA_product_quotients.m
code/agent_m18_416_R8_VB_certificate_check.m
data/agent_m18_416_R8_SA_product_quotients.log
data/agent_m18_416_R8_VB_certificate_check.log
```

## Key structural point

On the lambda-line (the conic S^2 = Gamma(X), X = m^2, with the
parametrization X(l), S(l) of the mwsieve note; reduced coordinate
satisfies (126*S)^2 = Gamma):

- h_B is EVEN in m, so V_B descends to the lambda-line (this is why the
  S_B V_4 analysis works);
- h_A is ODD in m (h_A = hAm/m), so V_A does NOT descend; but its
  m-conjugate norm does:
      (y_A * y_A')^2 = (2 alpha_A)^2 - X*bb^2 = d_A(X(l)),
  giving NEW lambda-line necessary conditions involving the A-side.

## Verified quotient list for R = -8

Genus-1 quotients (reproduce the mwsieve note exactly, all rank 1):

```text
X        -> E_m  = [0,0,0,-1263,17138],   rank [1,1], torsion [2,2]
V_B      -> E_mY = [1,1,1,-326,1874],     rank [1,1], torsion [2,2]
X*V_B    -> E_Y  = [0,0,0,-2847,29986],   rank [1,1], torsion [2,2]
```

NEW genus >= 2 quotients (point-finite by Faltings; Chabauty targets):

```text
d_A          : y^2 = deg-8 squarefree  (genus 3 hyperelliptic)
X*d_A        : y^2 = deg-11 squarefree
d_A*V_B      : y^2 = deg-11 squarefree
X*d_A*V_B    : y^2 = deg-12 squarefree
```

Explicit models in the log.  Any of these with provably finitely many
FOUND points forces the lambda values of any halving point into a
finite checked set: a certificate route for R=-8 that does not need a
rank-0 elliptic quotient.

## Validation and a cautionary note

`agent_m18_416_R8_VB_certificate_check.m` validates the V_B descent
against the exact number-field S_B test at 295 sample lambdas (295/295
agree, including the known S_B-pass points l = -5/7, -11/9).

An intermediate version of the quotient script DROPPED THE LEADING UNIT
returned by Magma's Factorization when forming squarefree models, which
produced a spurious rank-0 "killer" quotient for V_B.  This is the
third unit-constant bug in this project (after c4 and C2): treat every
Magma Factorization as monic-factors-times-unit and keep the unit.

## Next step for R = -8

Chabauty on the corrected deg-8 d_A quotient (genus 3): compute rank
bounds of its Jacobian (rank < 3 suffices), then Magma Chabauty to list
its rational points; intersect with the rank-1 elliptic quotients'
constraints.  If the finite lambda set contains only the known
degenerate points (disc_A = 0 at l = -5/7, -11/9 and boundaries), the
R=-8 fiber is certified empty and ALL THREE ELS fibers are closed.

## ERRATUM (same day)

The "d_A character quotients" above are NOT valid necessary conditions.
Because h_A is odd in m, y_A^2 = 2*alpha_A + m*bb generates a NON-ABELIAN
(D_4) quartic extension of Q(lambda); there is no (Z/2)^3 character group,
and a rational halving point does not map to y^2 = d_A(X(l)):
at such a point d_A == V_A' mod squares, which is not forced square.
Do NOT run Chabauty on the deg-8/11/11/12 models above.

The correct A-side necessary curve is the quartic

    C_A :  y^4 - 4*alpha_A(X(l))*y^2 + d_A(X(l)) = 0,

whose rational points automatically have m = (y^2 - 2*alpha_A)/bb
rational; C_A is a double cover of E_m (z = y^2 gives E_m), so
Jac(C_A) ~ E_m x Prym.  The Prym factors are the objects to bound.
See agent_m18_416_R8_CA_curve.m.

## Correct A/B condition curves on the X-line (2026-07-02, later)

Script: `code/agent_m18_416_R8_CA_curve.m` (log: `data/..._R8_CA_curve.log`).

The corrected necessary curves live on the X-line (X = m^2); every
halving point on the R=-8 fiber maps to a rational point of each.

**C_A^X** (the S_A condition), genus 3, NON-hyperelliptic:

```text
6561/1048576*X^8 - 405/8*X^7 - 81/512*X^6*y^2 + 183310047/1024*X^6
+ 639*X^5*y^2 - 361374048*X^5 + X^4*y^4 - 966735*X^4*y^2
+ 3644251097571/8*X^4 + 649264896*X^3*y^2 - 367179160707072*X^3
- 163326699648*X^2*y^2 + 184811265363146688*X^2
- 53104121520366551040*X + 6668902704477000830976 = 0
```

Canonical plane quartic model (coordinates x:y:z):

```text
1048576*x^4 - 165888*x^2*y^2 + 6561*y^4 + 1338900480*x^2*y*z
- 105992064*y^3*z - 2700412452864*x^2*z^2 + 641912739840*y^2*z^2
- 1727319315972096*y*z^3 + 1742534499206430720*z^4 = 0
```

Its z = y^2 quotient is the NEW (sixth) rank-1 elliptic curve

```text
E'_A: v^2 = X^3 - 1824*X^2 + 1016064*X,  min model [0,0,0,-363,41078],
rank bounds [1,1], torsion [2]
```

so Jac(C_A^X) ~ E'_A x Prym (2-dim).  Chabauty on C_A^X needs
rank Jac <= 2, i.e. rank(Prym) <= 1: computing the Prym rank is the
decisive remaining computation.  Known small points of C_A^X: the four
degenerate ones above X = 784, 1296 (d_A = 0, singular genus-2 fiber)
and points at infinity.

**C_B^X** (the S_B condition), genus 1:

```text
-9/8*X^2*y^2 - 81/8*X^2 + X*y^4 + 2178*X*y^2 + 22437*X
- 1143072*y^2 - 10287648 = 0
```

elliptic via (784, 3): minimal model [1,1,1,-326,1874] = E_mY again,
rank [1,1].  No new information from the B side on the X-line.

## R = -8 quotient inventory (all rank 1)

```text
E_m   = [0,0,0,-1263,17138]      (lambda: X character)
E_mY  = [1,1,1,-326,1874]        (lambda: V_B character; also C_B^X)
E_Y   = [0,0,0,-2847,29986]      (lambda: X*V_B character)
E'_A  = [0,0,0,-363,41078]       (X-line: z-quotient of C_A^X)  NEW
```

No rank-zero quotient exists among all elliptic quotients found so far.
The certificate route for R=-8 is genus-3 machinery on C_A^X:
 1. compute rank of the 2-dim Prym of C_A^X -> E'_A (<= 1 needed);
 2. then classical/quadratic Chabauty on the plane quartic;
 3. alternatively TwoCoverDescent-style analysis of the full fiber
    product with C_B.
