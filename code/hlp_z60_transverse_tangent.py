#!/usr/bin/env python3
"""Exact tangent calculation at the split HLP cyclic-Z/60 seed.

All arithmetic is over Q using Fraction, and all polynomials are stored from
low to high coefficient.  The marked layers are

  * an order-5 norm identity A5^2-F*B5^2 = k5*q5^5;
  * an order-3 contact identity A3^2-F = k3*q3^3;
  * an order-4 halving identity ell^2-F = k4*u4^2*q2,
    with ell=q2*L and 2D4=[q2,0].

The scale of the order-5 principal function is fixed by B5(0)=1.  The
auxiliary Jacobian blocks are square.  Their invertibility proves that the
marked torsion incidence is etale over the seven sextic coefficients at the
seed.  We then compare its tangent space with the unique Humbert-4 branch.
"""

from fractions import Fraction as Q
from functools import reduce
from math import gcd, lcm


def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def pad(a, n):
    return list(a) + [Q(0)] * (n - len(a))


def add(a, b):
    n = max(len(a), len(b))
    return trim([(a[i] if i < len(a) else 0)
                 + (b[i] if i < len(b) else 0) for i in range(n)])


def sub(a, b):
    n = max(len(a), len(b))
    return trim([(a[i] if i < len(a) else 0)
                 - (b[i] if i < len(b) else 0) for i in range(n)])


def scale(a, c):
    return trim([Q(c) * x for x in a])


def mul(a, b):
    out = [Q(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return trim(out)


def power(a, n):
    out = [Q(1)]
    for _ in range(n):
        out = mul(out, a)
    return out


def divrem(a, b):
    a, b = trim(a), trim(b)
    quo = [Q(0)] * max(1, len(a) - len(b) + 1)
    while a != [0] and len(a) >= len(b):
        j = len(a) - len(b)
        c = a[-1] / b[-1]
        quo[j] = c
        for i, bi in enumerate(b):
            a[i + j] -= c * bi
        a = trim(a)
    return trim(quo), a


def coefficient_matrix(columns, nrows):
    return [[pad(columns[j], nrows)[i] for j in range(len(columns))]
            for i in range(nrows)]


def rref(matrix):
    a = [[Q(x) for x in row] for row in matrix]
    nr = len(a)
    nc = len(a[0]) if nr else 0
    row = 0
    pivots = []
    for col in range(nc):
        pivot = next((i for i in range(row, nr) if a[i][col]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        c = a[row][col]
        a[row] = [x / c for x in a[row]]
        for i in range(nr):
            if i != row and a[i][col]:
                c = a[i][col]
                a[i] = [a[i][j] - c * a[row][j] for j in range(nc)]
        pivots.append(col)
        row += 1
        if row == nr:
            break
    return row, pivots, a


def rank(matrix):
    return rref(matrix)[0]


def solve_square(matrix, rhs):
    n = len(matrix)
    assert all(len(row) == n for row in matrix)
    rr, pivots, red = rref([matrix[i] + [Q(rhs[i])] for i in range(n)])
    assert rr == n and pivots[:n] == list(range(n))
    return [red[i][-1] for i in range(n)]


def primitive_integer_vector(v):
    den = reduce(lcm, (x.denominator for x in v), 1)
    nums = [int(x * den) for x in v]
    content = reduce(gcd, (abs(x) for x in nums if x), 0)
    nums = [x // content for x in nums]
    first = next((x for x in nums if x), 1)
    return [-x for x in nums] if first < 0 else nums


def order5_auxiliary_jacobian(f, q, a, b, k):
    # Variables q0,q1,A0..A5,B1,B2,k; B0=1 fixes function scaling.
    columns = [scale(power(q, 4), -5 * k),
               scale(mul([0, 1], power(q, 4)), -5 * k)]
    columns += [scale(mul(a, [0] * i + [1]), 2) for i in range(6)]
    columns += [scale(mul(mul(f, b), [0, 1]), -2),
                scale(mul(mul(f, b), [0, 0, 1]), -2)]
    columns += [scale(power(q, 5), -1)]
    return coefficient_matrix(columns, 11)


def order3_auxiliary_jacobian(q, a, k):
    # Variables q0,q1,A0..A3,k; the y coefficient is scaled to one.
    columns = [scale(power(q, 2), -3 * k),
               scale(mul([0, 1], power(q, 2)), -3 * k)]
    columns += [scale(mul(a, [0] * i + [1]), 2) for i in range(4)]
    columns += [scale(power(q, 3), -1)]
    return coefficient_matrix(columns, 7)


def order4_auxiliary_jacobian(q2, u, line, ell, k):
    # Variables q20,q21,u0,u1,L0,L1,k.
    columns = []
    for dq in ([1], [0, 1]):
        dell = mul(dq, line)
        columns.append(sub(scale(mul(ell, dell), 2),
                           scale(mul(power(u, 2), dq), k)))
    for du in ([1], [0, 1]):
        columns.append(scale(mul(mul(u, du), q2), -2 * k))
    for dl in ([1], [0, 1]):
        columns.append(scale(mul(ell, mul(q2, dl)), 2))
    columns.append(scale(mul(power(u, 2), q2), -1))
    return coefficient_matrix(columns, 7)


def curve_columns_order5(b):
    b2 = power(b, 2)
    columns = [scale(mul(b2, [0] * j + [1]), -1) for j in range(7)]
    return coefficient_matrix(columns, 11)


def curve_columns_plain():
    return [[Q(-1) if i == j else Q(0) for j in range(7)]
            for i in range(7)]


def main():
    # HLP exact cyclic-[60] seed, low coefficients.
    f = [Q(94277468160), 0, Q(-22332312000), 0,
         Q(1761500625), 0, Q(-46250000)]

    # 12*T60: exact order 5.  Phi5=A5+B5*y vanishes to order five
    # on q5, and Norm(Phi5)=A5^2-F*B5^2=k5*q5^5.
    q5 = [Q(-4608, 395), 0, 1]
    v5 = [0, Q(-164280, 79)]
    b5 = [1, 0, Q(-125, 1728)]
    a5 = [0, Q(66600), 0, Q(-1558625, 144),
          0, Q(48752125, 110592)]
    k5 = Q(1778923230671875, 4076863488)

    # 20*T60: exact order 3.
    q3 = [Q(-316, 25), 0, 1]
    v3 = [Q(-5476)]
    a3 = [Q(-29600), 0, Q(2775)]
    k3 = Q(46250000)

    # 15*T60: exact order 4; its double is [q2,0].
    q2 = [Q(-1728, 125), 0, 1]
    u4 = [Q(-1506, 125), 0, 1]
    v4 = [Q(-24642, 5)]
    line4 = [Q(2775)]
    ell4 = mul(q2, line4)
    k4 = Q(46250000)

    assert sub(sub(power(a5, 2), mul(f, power(b5, 2))),
               scale(power(q5, 5), k5)) == [0]
    assert divrem(add(a5, mul(b5, v5)), q5)[1] == [0]
    assert sub(sub(power(a3, 2), f), scale(power(q3, 3), k3)) == [0]
    assert divrem(add(a3, v3), q3)[1] == [0]
    assert sub(sub(power(ell4, 2), f),
               scale(mul(power(u4, 2), q2), k4)) == [0]
    assert divrem(sub(ell4, v4), u4)[1] == [0]

    j5 = order5_auxiliary_jacobian(f, q5, a5, b5, k5)
    j3 = order3_auxiliary_jacobian(q3, a3, k3)
    j4 = order4_auxiliary_jacobian(q2, u4, line4, ell4, k4)

    # Full system: seven curve variables, followed by the 11+7+7
    # auxiliary variables.  Rows are the 11+7+7 norm coefficients.
    cf5 = curve_columns_order5(b5)
    cfp = curve_columns_plain()
    full = []
    for i in range(11):
        full.append(cf5[i] + j5[i] + [Q(0)] * 14)
    for i in range(7):
        full.append(cfp[i] + [Q(0)] * 11 + j3[i] + [Q(0)] * 7)
    for i in range(7):
        full.append(cfp[i] + [Q(0)] * 18 + j4[i])

    print("HLP Z60 EXACT MARKED TANGENT")
    print("identity5 q", q5, "A", a5, "B", b5, "k", k5)
    print("identity3 q", q3, "A", a3, "k", k3)
    print("identity4 q2", q2, "u", u4, "ell", ell4, "k", k4)
    print("AUXILIARY_BLOCK_RANKS", rank(j5), rank(j3), rank(j4))
    print("FULL_JACOBIAN_RANK", rank(full), "AMBIENT", 32,
          "TANGENT_DIM", 32 - rank(full))
    print("PROJECTION_TO_SEXTIC_RANK", 7)

    # Infinitesimal equation scaling and PGL2 orbit of the binary sextic.
    translation = [Q(0)] * 7       # f'(x)
    dilation = [Q(0)] * 7          # x*f'(x)
    projective = [Q(0)] * 7        # x^2*f'(x)-6*x*f(x)
    for i in range(1, 7):
        translation[i - 1] = i * f[i]
        dilation[i] = i * f[i]
    for i in range(7):
        if i + 1 <= 6:
            projective[i + 1] = (i - 6) * f[i]
    gauge = [f, translation, dilation, projective]
    print("GAUGE_RANK_PGL2_PLUS_Y_SCALING", rank(gauge))
    print("MARKED_MODULI_TANGENT_DIM", 7 - rank(gauge))

    # The unique reduced geometric involution is x -> -x.  Its H_4 branch
    # tangent is the even-sextic subspace plus the PGL2 orbit.
    even = []
    for i in (0, 2, 4, 6):
        e = [Q(0)] * 7
        e[i] = 1
        even.append(e)
    split_tangent = even + [translation, dilation, projective]
    print("SPLIT_TANGENT_RANK", rank(split_tangent))
    print("SPLIT_MODULI_TANGENT_DIM", rank(split_tangent) - rank(gauge))

    # Exact primitive normal to even+PGL2.
    odd_translation = [translation[i] for i in (1, 3, 5)]
    odd_projective = [projective[i] for i in (1, 3, 5)]
    cross = [odd_translation[1] * odd_projective[2]
             - odd_translation[2] * odd_projective[1],
             odd_translation[2] * odd_projective[0]
             - odd_translation[0] * odd_projective[2],
             odd_translation[0] * odd_projective[1]
             - odd_translation[1] * odd_projective[0]]
    odd_normal = primitive_integer_vector(cross)
    normal = [0, odd_normal[0], 0, odd_normal[1], 0, odd_normal[2], 0]
    df = [Q(1), Q(1), 0, 0, 0, 0, 0]
    pairing = sum(Q(normal[i]) * df[i] for i in range(7))
    print("UNIQUE_HUMBERT4_BRANCH_NORMAL", normal)
    print("TRANSVERSE_DF", df, "NORMAL_PAIRING", pairing)
    print("DF_PLUS_GAUGE_RANK", rank(gauge + [df]))

    # Unique rational lift of df through all three invertible blocks.
    rhs5 = pad(mul(df, power(b5, 2)), 11)
    rhs3 = pad(df, 7)
    rhs4 = pad(df, 7)
    sol5 = solve_square(j5, rhs5)
    sol3 = solve_square(j3, rhs3)
    sol4 = solve_square(j4, rhs4)
    names5 = ["dq50", "dq51"] + [f"dA5{i}" for i in range(6)] \
        + ["dB51", "dB52", "dk5"]
    names3 = ["dq30", "dq31"] + [f"dA3{i}" for i in range(4)] \
        + ["dk3"]
    names4 = ["dq20", "dq21", "du40", "du41", "dL40", "dL41", "dk4"]
    print("RATIONAL_LIFT_ORDER5", list(zip(names5, sol5)))
    print("RATIONAL_LIFT_ORDER3", list(zip(names3, sol3)))
    print("RATIONAL_LIFT_ORDER4", list(zip(names4, sol4)))
    print("CONCLUSION rational non-split tangent exists")


if __name__ == "__main__":
    main()
