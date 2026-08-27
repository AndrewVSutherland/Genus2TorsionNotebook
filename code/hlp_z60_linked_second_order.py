#!/usr/bin/env python3
"""Second-order test of the linked B5=q0 locus at the HLP Z/60 seed.

The linked marked incidence is written with 23 variables and 18 equations.
The order-4 identity eliminates the seven sextic coefficients:

    f = (q0*L)^2 - k4*q0*u^2,

and B5=q0 is imposed in the order-5 norm identity.  Exact two-jets over Q
then give a genuine (non-PGL2) tangent direction and its forced quadratic
term.  Finally we ask whether a trace-zero projective involution near
x -> -x can preserve the resulting binary sextic through order t^2.

Only small dense Fraction matrices (at most 23 by 24) are used.
"""

from fractions import Fraction as Q

from hlp_z60_marked_tangent import rank, rref, solve_square


class Jet:
    """A scalar modulo t^3, with coefficients in Q."""

    __slots__ = ("c",)

    def __init__(self, c0=0, c1=0, c2=0):
        self.c = (Q(c0), Q(c1), Q(c2))

    @staticmethod
    def coerce(a):
        return a if isinstance(a, Jet) else Jet(a)

    def __add__(self, other):
        other = Jet.coerce(other)
        return Jet(*(self.c[i] + other.c[i] for i in range(3)))

    __radd__ = __add__

    def __neg__(self):
        return Jet(*(-a for a in self.c))

    def __sub__(self, other):
        return self + (-Jet.coerce(other))

    def __rsub__(self, other):
        return Jet.coerce(other) - self

    def __mul__(self, other):
        other = Jet.coerce(other)
        a, b = self.c, other.c
        return Jet(a[0] * b[0],
                   a[0] * b[1] + a[1] * b[0],
                   a[0] * b[2] + a[1] * b[1] + a[2] * b[0])

    __rmul__ = __mul__

    def __pow__(self, n):
        assert n >= 0
        ans, base = Jet(1), self
        while n:
            if n & 1:
                ans = ans * base
            base = base * base
            n //= 2
        return ans


def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) +
            (b[i] if i < len(b) else 0) for i in range(n)]


def psub(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) -
            (b[i] if i < len(b) else 0) for i in range(n)]


def pmul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] = out[i + j] + ai * bj
    return out


def ppow(a, n):
    out = [1]
    for _ in range(n):
        out = pmul(out, a)
    return out


def pscale(a, c):
    return [c * ai for ai in a]


def ppad(a, n):
    return a + [0] * (n - len(a))


def unpack(z):
    """Return q0,u,L,k4,q5,A,k5,q3,H,k3 from the 23 variables."""
    q0 = [z[0], z[1], 1]
    u = [z[2], z[3], 1]
    line = [z[4], z[5]]
    k4 = z[6]
    q5 = [z[7], z[8], 1]
    a5 = z[9:15]
    k5 = z[15]
    q3 = [z[16], z[17], 1]
    h3 = z[18:22]
    k3 = z[22]
    return q0, u, line, k4, q5, a5, k5, q3, h3, k3


def curve(z):
    q0, u, line, k4, *_ = unpack(z)
    ell = pmul(q0, line)
    return psub(pmul(ell, ell), pscale(pmul(q0, pmul(u, u)), k4))


def equations(z):
    q0, _u, _line, _k4, q5, a5, k5, q3, h3, k3 = unpack(z)
    f = curve(z)
    e5 = psub(psub(pmul(a5, a5), pmul(f, pmul(q0, q0))),
                pscale(ppow(q5, 5), k5))
    e3 = psub(psub(pmul(h3, h3), f), pscale(ppow(q3, 3), k3))
    return ppad(e5, 11)[:11] + ppad(e3, 7)[:7]


def seed():
    a5 = pscale(pmul([0, 1], [Q(7962624), 0, Q(-1294080), 0,
                               Q(52705)]), Q(37, 320))
    return [
        Q(-1728, 125), 0, Q(-1506, 125), 0, Q(-2775), 0,
        Q(46250000),
        Q(-4608, 395), 0, *a5, Q(341553260289, 4096),
        Q(-316, 25), 0, Q(29600), 0, Q(-2775), 0, Q(46250000),
    ]


def jacobian(z0):
    ans = [[Q(0) for _ in z0] for _ in range(18)]
    for j in range(len(z0)):
        zj = [Jet(a, int(i == j), 0) for i, a in enumerate(z0)]
        ej = equations(zj)
        for i in range(18):
            ans[i][j] = ej[i].c[1]
    return ans


def null_basis(a):
    rr, pivots, red = rref(a)
    assert rr == len(a)
    free = [j for j in range(len(a[0])) if j not in pivots]
    basis = []
    for j in free:
        v = [Q(0)] * len(a[0])
        v[j] = 1
        for i, p in enumerate(pivots):
            v[p] = -red[i][j]
        basis.append(v)
    return pivots, free, basis


def curve_jet(z0, v, w=None):
    if w is None:
        w = [Q(0)] * len(z0)
    zj = [Jet(z0[i], v[i], w[i]) for i in range(len(z0))]
    return curve(zj)


def gauge_vectors(f):
    translation = [Q(0)] * 7
    dilation = [Q(0)] * 7
    projective = [Q(0)] * 7
    for i in range(1, 7):
        translation[i - 1] = i * f[i]
        dilation[i] = i * f[i]
    for i in range(7):
        if i + 1 <= 6:
            projective[i + 1] = (i - 6) * f[i]
    return [f, translation, dilation, projective]


def solve_overdetermined(columns, rhs):
    """Solve columns*x=rhs, returning (consistent, one solution)."""
    m = [[columns[j][i] for j in range(len(columns))] + [rhs[i]]
         for i in range(len(rhs))]
    rr, pivots, red = rref(m)
    n = len(columns)
    if any(all(red[i][j] == 0 for j in range(n)) and red[i][n] != 0
           for i in range(rr)):
        return False, []
    sol = [Q(0)] * n
    for i, p in enumerate(pivots):
        if p < n:
            sol[p] = red[i][n]
    return True, sol


def binomial(n, k):
    if k < 0 or k > n:
        return 0
    out = 1
    for j in range(1, k + 1):
        out = out * (n + 1 - j) // j
    return out


def binary_transform(f, b, c):
    """F(-X+bZ,cX+Z), coefficients ordered by powers of X."""
    out = [Jet(0) for _ in range(7)]
    for i, fi in enumerate(f):
        # (-X+bZ)^i (cX+Z)^(6-i)
        for r in range(i + 1):
            left = binomial(i, r) * ((-1) ** r) * (b ** (i - r))
            for s in range(7 - i):
                power_x = r + s
                right = binomial(6 - i, s) * (c ** s)
                out[power_x] = out[power_x] + fi * left * right
    return out


def involution_residual(f, b, c):
    lhs = binary_transform(f, b, c)
    multiplier = (Jet(1) + b * c) ** 3
    return [lhs[i] - multiplier * f[i] for i in range(7)]


def primitive(v):
    from functools import reduce
    from math import gcd, lcm
    den = reduce(lcm, (a.denominator for a in v), 1)
    nums = [int(a * den) for a in v]
    content = reduce(gcd, (abs(a) for a in nums if a), 0)
    if content:
        nums = [a // content for a in nums]
    first = next((a for a in nums if a), 1)
    return [-a for a in nums] if first < 0 else nums


def main():
    z0 = seed()
    e0 = equations(z0)
    assert all(a == 0 for a in e0)
    j = jacobian(z0)
    pivots, free, basis = null_basis(j)
    f0 = [Q(a) for a in curve(z0)]
    gauge = gauge_vectors(f0)

    print("HLP_Z60_LINKED_SECOND_ORDER")
    print("LINKED_JACOBIAN_RANK", rank(j), "VARIABLES", len(z0),
          "LOCAL_DIMENSION", len(free))
    print("FREE_VARIABLE_INDICES", free)

    chosen = None
    full_df_basis = []
    for free_index, v in zip(free, basis):
        df = [a.c[1] for a in curve_jet(z0, v)]
        full_df_basis.append(df)
        moduli = rank(gauge + [df]) > rank(gauge)
        print("NULL_DIRECTION", free_index, "DF", primitive(df),
              "NON_GAUGE", moduli)
        if chosen is None and moduli:
            chosen = (free_index, v, df)
    assert chosen is not None
    free_index, v, df = chosen

    # The parity-preserving locus has ten even equations in thirteen variables.
    # Its tangent plus the two odd PGL2 directions fills the linked tangent.
    parity_variables = [0, 2, 4, 6, 7, 10, 12, 14, 15,
                        16, 18, 20, 22]
    parity_equations = [0, 2, 4, 6, 8, 10, 11, 13, 15, 17]
    jpar = [[j[i][k] for k in parity_variables]
            for i in parity_equations]
    _ppiv, pfree, pbasis = null_basis(jpar)
    parity_df_basis = []
    for pv in pbasis:
        vv = [Q(0)] * len(z0)
        for k, value in zip(parity_variables, pv):
            vv[k] = value
        parity_df_basis.append([a.c[1] for a in curve_jet(z0, vv)])
    odd_pgl = [gauge[1], gauge[3]]
    full_df_rank = rank(full_df_basis)
    saturation_rank = rank(parity_df_basis + odd_pgl)
    print("PARITY_SYSTEM_RANK", rank(jpar), "VARIABLES",
          len(parity_variables), "LOCAL_DIMENSION", len(pfree))
    print("PARITY_DF_RANK", rank(parity_df_basis))
    print("FULL_LINKED_DF_RANK", full_df_rank)
    print("PARITY_PLUS_ODD_PGL_RANK", saturation_rank)
    print("SATURATION_FILLS_LINKED_TANGENT",
          saturation_rank == full_df_rank == 5)

    # Solve the t^2 equations while keeping all five free variables affine:
    # the selected one equals its seed plus t, and the other four are fixed.
    z_linear = [Jet(z0[i], v[i], 0) for i in range(len(z0))]
    quadratic_rhs = [-a.c[2] for a in equations(z_linear)]
    jpiv = [[j[i][p] for p in pivots] for i in range(18)]
    wpiv = solve_square(jpiv, quadratic_rhs)
    w = [Q(0)] * len(z0)
    for p, value in zip(pivots, wpiv):
        w[p] = value
    assert all(a.c[2] == 0 for a in equations(
        [Jet(z0[i], v[i], w[i]) for i in range(len(z0))]))
    fjet = curve_jet(z0, v, w)

    normal = [Q(0), Q(81125), 0, Q(904800), 0, Q(9916416), 0]
    print("CHOSEN_FREE_VARIABLE", free_index)
    print("FIRST_DF", primitive([a.c[1] for a in fjet]))
    print("FIRST_H4_NORMAL", sum(normal[i] * fjet[i].c[1]
                                  for i in range(7)))

    # First-order preserving involution M=[-1,b;c,1].
    base1 = [a.c[1] for a in involution_residual(
        fjet, Jet(0), Jet(0))]
    col_b = [a.c[1] for a in involution_residual(
        [Jet(a.c[0]) for a in fjet], Jet(0, 1), Jet(0))]
    col_c = [a.c[1] for a in involution_residual(
        [Jet(a.c[0]) for a in fjet], Jet(0), Jet(0, 1))]
    ok1, bc1 = solve_overdetermined([col_b, col_c], [-a for a in base1])
    assert ok1
    b1, c1 = bc1

    # At t^2, b2 and c2 enter through the same two linear columns.
    base2 = [a.c[2] for a in involution_residual(
        fjet, Jet(0, b1), Jet(0, c1))]
    ok2, bc2 = solve_overdetermined([col_b, col_c], [-a for a in base2])
    print("FIRST_INVOLUTION_BC", b1, c1)
    print("SECOND_INVOLUTION_SOLVABLE", ok2)
    if ok2:
        residual2 = [Q(0)] * 7
        print("SECOND_INVOLUTION_BC", bc2[0], bc2[1])
    else:
        # Cancel two independent rows and display the remaining exact
        # obstruction vector.  Its nonzero primitive form is certificate
        # enough that no b2,c2 can extend the involution.
        pivot_rows = None
        for i in range(7):
            for k in range(i + 1, 7):
                m = [[col_b[i], col_c[i]], [col_b[k], col_c[k]]]
                if rank(m) == 2:
                    pivot_rows = (i, k, m)
                    break
            if pivot_rows:
                break
        i, k, m = pivot_rows
        correction = solve_square(m, [-base2[i], -base2[k]])
        residual2 = [base2[r] + col_b[r] * correction[0]
                     + col_c[r] * correction[1] for r in range(7)]
        print("SECOND_BC_FROM_PIVOT_ROWS", i, k,
              correction[0], correction[1])
        print("SECOND_OBSTRUCTION_PRIMITIVE", primitive(residual2))
        print("SECOND_OBSTRUCTION_NONZERO_ROWS",
              [i for i, a in enumerate(residual2) if a])

    fills = saturation_rank == full_df_rank == 5
    if fills and ok2:
        conclusion = "LINKED_GERM_LOCALLY_CONTAINED_IN_H4"
    elif not ok2:
        conclusion = "LEAVES_H4_AT_SECOND_ORDER"
    else:
        conclusion = "NO_SECOND_ORDER_ESCAPE_DETECTED"
    print("CONCLUSION", conclusion)
    print("HLP_Z60_LINKED_SECOND_ORDER_DONE")


if __name__ == "__main__":
    main()
