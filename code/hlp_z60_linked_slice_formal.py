#!/usr/bin/env python3
"""Formal linked-slice test against the split (Humbert-4) branch.

We impose the exact linked conditions

    B5=q0,  q0=x^2-1728/125,  L4(0)=-2775,  k4=46250000.

After eliminating F through the order-4 identity, the remaining marked
order-5 and order-3 system has 19 variables and 18 equations.  It is a
smooth curve at the cyclic-60 HLP point.  This script constructs its unique
formal branch after choosing one local parameter, one coefficient at a
time.  At each order it also tries to lift the involution x |-> -x in the
form

    (X:Z) |-> (-X+b Z : c X+Z).

Failure of the seven binary-sextic invariance equations detects the first
order at which the linked curve leaves the Humbert-4 branch.  Arithmetic is
exact and the matrices are at most 19 by 19.
"""

from fractions import Fraction as Q
from math import comb

from hlp_z60_marked_tangent import rank, rref, solve_square


ORDER = 30


def spad(a, n):
    return list(a) + [Q(0)] * (n - len(a))


def sadd(a, b, n):
    aa, bb = spad(a, n + 1), spad(b, n + 1)
    return [aa[i] + bb[i] for i in range(n + 1)]


def sscale(a, c, n):
    aa = spad(a, n + 1)
    return [Q(c) * aa[i] for i in range(n + 1)]


def smul(a, b, n):
    aa, bb = spad(a, n + 1), spad(b, n + 1)
    out = [Q(0)] * (n + 1)
    for i in range(n + 1):
        for j in range(n + 1 - i):
            out[i + j] += aa[i] * bb[j]
    return out


def spow(a, e, n):
    out = [Q(1)] + [Q(0)] * n
    for _ in range(e):
        out = smul(out, a, n)
    return out


def xpad(a, n):
    zero = [Q(0)] * (n + 1)
    return [list(c) for c in a] + [zero[:] for _ in range(n - len(a))]


def xsub(a, b, sn):
    m = max(len(a), len(b))
    aa, bb = xpad(a, m), xpad(b, m)
    return [sadd(aa[i], sscale(bb[i], -1, sn), sn) for i in range(m)]


def xmul(a, b, sn):
    out = [[Q(0)] * (sn + 1) for _ in range(len(a) + len(b) - 1)]
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] = sadd(out[i + j], smul(ai, bj, sn), sn)
    return out


def xpow(a, e, sn):
    out = [[Q(1)] + [Q(0)] * sn]
    for _ in range(e):
        out = xmul(out, a, sn)
    return out


def xscale(a, s, sn):
    return [smul(c, s, sn) for c in a]


def const_poly(values, sn):
    return [[Q(v)] + [Q(0)] * sn for v in values]


def curve_and_equations(z, sn):
    # z = u0,u1,l1, q50,q51,A0..A5,k5, q30,q31,H0..H3,k3
    u0, u1, l1 = z[0:3]
    q50, q51 = z[3:5]
    aa = z[5:11]
    k5 = z[11]
    q30, q31 = z[12:14]
    hh = z[14:18]
    k3 = z[18]
    q0 = const_poly([Q(-1728, 125), 0, 1], sn)
    uu = [u0, u1, [Q(1)] + [Q(0)] * sn]
    line = [const_poly([Q(-2775)], sn)[0], l1]
    ell = xmul(q0, line, sn)
    ff = xsub(xpow(ell, 2, sn),
              xscale(xmul(q0, xpow(uu, 2, sn), sn),
                     [Q(46250000)] + [Q(0)] * sn, sn), sn)

    q5 = [q50, q51, [Q(1)] + [Q(0)] * sn]
    e5 = xsub(xsub(xpow(aa, 2, sn),
                    xmul(ff, xpow(q0, 2, sn), sn), sn),
              xscale(xpow(q5, 5, sn), k5, sn), sn)
    q3 = [q30, q31, [Q(1)] + [Q(0)] * sn]
    e3 = xsub(xsub(xpow(hh, 2, sn), ff, sn),
              xscale(xpow(q3, 3, sn), k3, sn), sn)
    e5 = xpad(e5, 11)
    e3 = xpad(e3, 7)
    return ff, e5 + e3


def seed():
    a = Q(37, 320)
    return [
        Q(-1506, 125), 0, 0,
        Q(-4608, 395), 0,
        0, a * 7962624, 0, -a * 1294080, 0, a * 52705,
        Q(341553260289, 4096),
        Q(-316, 25), 0,
        Q(29600), 0, Q(-2775), 0, Q(46250000),
    ]


def jacobian(z0):
    out = [[Q(0)] * len(z0) for _ in range(18)]
    for j in range(len(z0)):
        z = [[v, Q(1) if i == j else Q(0)] for i, v in enumerate(z0)]
        _f, ee = curve_and_equations(z, 1)
        for i in range(18):
            out[i][j] = ee[i][1]
    return out


def kernel_vector(jac):
    rr, pivots, red = rref(jac)
    assert rr == 18
    free = next(j for j in range(19) if j not in pivots)
    v = [Q(0)] * 19
    v[free] = 1
    for i, p in enumerate(pivots):
        v[p] = -red[i][free]
    return free, v


def solve_rectangular(mat, rhs):
    aug = [list(mat[i]) + [rhs[i]] for i in range(len(mat))]
    r0 = rank(mat)
    r1 = rank(aug)
    if r1 != r0:
        return None, r0, r1
    _rr, pivots, red = rref(aug)
    assert r0 == len(mat[0])
    sol = [Q(0)] * len(mat[0])
    for i, p in enumerate(pivots):
        if p < len(sol):
            sol[p] = red[i][-1]
    return sol, r0, r1


def formal_marked_branch(max_order):
    z0 = seed()
    jac = jacobian(z0)
    free, z1 = kernel_vector(jac)
    coeffs = [[v] for v in z0]
    for i in range(19):
        coeffs[i].append(z1[i])
    square = jac + [[Q(1) if j == free else Q(0) for j in range(19)]]
    for n in range(2, max_order + 1):
        trial = [spad(c, n + 1) for c in coeffs]
        _f, ee = curve_and_equations(trial, n)
        rhs = [-ee[i][n] for i in range(18)] + [Q(0)]
        zn = solve_square(square, rhs)
        for i in range(19):
            coeffs[i].append(zn[i])
    ff, ee = curve_and_equations(coeffs, max_order)
    assert all(all(v == 0 for v in e) for e in ee)
    return free, coeffs, xpad(ff, 7)


def transformed_binary_sextic(ff, b, c, lam, sn):
    # F(-X+bZ,cX+Z)-lam*F(X,Z), coefficient of X^j Z^(6-j).
    out = [[Q(0)] * (sn + 1) for _ in range(7)]
    minus_one = [Q(-1)] + [Q(0)] * sn
    one = [Q(1)] + [Q(0)] * sn
    for i in range(7):
        for r in range(i + 1):
            left = sscale(spow(b, i - r, sn), comb(i, r) * ((-1) ** r), sn)
            for ss in range(6 - i + 1):
                right = sscale(spow(c, ss, sn), comb(6 - i, ss), sn)
                j = r + ss
                term = smul(ff[i], smul(left, right, sn), sn)
                out[j] = sadd(out[j], term, sn)
    for j in range(7):
        out[j] = sadd(out[j], sscale(smul(lam, ff[j], sn), -1, sn), sn)
    return out


def lift_involution(ff, max_order):
    # First-order matrix for the new b_n,c_n,lambda_n at every n.
    zero = [Q(0), Q(0)]
    one = [Q(1), Q(0)]
    base = transformed_binary_sextic(
        [spad(f, 2) for f in ff], zero, zero, one, 1)
    cols = []
    for which in range(3):
        b = [Q(0), Q(1) if which == 0 else Q(0)]
        c = [Q(0), Q(1) if which == 1 else Q(0)]
        lam = [Q(1), Q(1) if which == 2 else Q(0)]
        val = transformed_binary_sextic(
            [spad(f, 2) for f in ff], b, c, lam, 1)
        cols.append([val[i][1] - base[i][1] for i in range(7)])
    mat = [[cols[j][i] for j in range(3)] for i in range(7)]
    bser, cser, lser = [Q(0)], [Q(0)], [Q(1)]
    for n in range(1, max_order + 1):
        bb, cc, ll = spad(bser, n + 1), spad(cser, n + 1), spad(lser, n + 1)
        val = transformed_binary_sextic(ff, bb, cc, ll, n)
        rhs = [-val[i][n] for i in range(7)]
        sol, r0, r1 = solve_rectangular(mat, rhs)
        if sol is None:
            return n, r0, r1, rhs, bser, cser, lser
        bser.append(sol[0])
        cser.append(sol[1])
        lser.append(sol[2])
    return None, rank(mat), rank(mat), [], bser, cser, lser


def main():
    free, coeffs, ff = formal_marked_branch(ORDER)
    normal = [Q(0), Q(81125), 0, Q(904800), 0, Q(9916416), 0]
    print("HLP_Z60_LINKED_SLICE_FORMAL")
    print("variables 19 equations 18 jacobian_rank 18")
    print("local_parameter_variable_index", free)
    for n in range(1, min(3, ORDER) + 1):
        pair = sum(normal[i] * ff[i][n] for i in range(7))
        print("raw_normal_pairing_order", n, pair)
    failure, r0, r1, rhs, b, c, lam = lift_involution(ff, ORDER)
    print("involution_linear_rank", r0)
    if failure is None:
        print("INVOLUTION_LIFTS_THROUGH_ORDER", ORDER)
        r = Q(-1728, 125)
        relation_b = all(spad(b, ORDER + 1)[i] == r * spad(c, ORDER + 1)[i]
                         for i in range(ORDER + 1))
        mu = sadd([Q(1)], sscale(spow(c, 2, ORDER), r, ORDER), ORDER)
        relation_lambda = spad(lam, ORDER + 1) == spow(mu, 3, ORDER)
        print("B_EQUALS_Q0_CONSTANT_TIMES_C", relation_b)
        print("LAMBDA_EQUALS_Q0_MULTIPLIER_CUBED", relation_lambda)
    else:
        print("INVOLUTION_OBSTRUCTION_ORDER", failure,
              "rank", r0, "augmented_rank", r1)
        print("obstruction_rhs", rhs)
    print("HLP_Z60_LINKED_SLICE_FORMAL_DONE")


if __name__ == "__main__":
    main()
