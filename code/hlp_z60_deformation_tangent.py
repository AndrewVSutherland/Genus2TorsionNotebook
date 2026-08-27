#!/usr/bin/env python3
"""Exact tangent calculation for simultaneous 5-, 3-, and 4-torsion.

The base point is the split HLP genus-2 curve whose rational Jacobian
torsion is cyclic of order 60.  Polynomial coefficient lists are stored
low-to-high and all arithmetic is over Q via Fraction.

The marked identities are

    A^2 - B^2 F = k5 q5^5,
    H^2 - F       = k3 q3^3,
    ell^2 - F     = k4 u4^2 q2,   ell=q2*L.

The first identity is put in the local gauge B(0)=1.  It therefore has
11 auxiliary variables and 11 coefficient equations.  The other two
blocks have seven variables and seven equations each.
"""

from fractions import Fraction as Q
from functools import reduce
from math import gcd, lcm


def trim(a):
    a = [Q(x) for x in a]
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def pad(a, n):
    return trim(a) + [Q(0)] * (n - len(trim(a)))


def add(a, b):
    n = max(len(a), len(b))
    return trim([(a[i] if i < len(a) else 0) +
                 (b[i] if i < len(b) else 0) for i in range(n)])


def sub(a, b):
    n = max(len(a), len(b))
    return trim([(a[i] if i < len(a) else 0) -
                 (b[i] if i < len(b) else 0) for i in range(n)])


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
        shift = len(a) - len(b)
        coeff = a[-1] / b[-1]
        quo[shift] = coeff
        for i, bi in enumerate(b):
            a[i + shift] -= coeff * bi
        a = trim(a)
    return trim(quo), a


def rref(matrix):
    a = [[Q(x) for x in row] for row in matrix]
    nr, nc = len(a), len(a[0])
    row = 0
    pivots = []
    for col in range(nc):
        pivot = next((i for i in range(row, nr) if a[i][col]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        lead = a[row][col]
        a[row] = [x / lead for x in a[row]]
        for i in range(nr):
            if i != row and a[i][col]:
                c = a[i][col]
                a[i] = [a[i][j] - c * a[row][j] for j in range(nc)]
        pivots.append(col)
        row += 1
        if row == nr:
            break
    return row, pivots, a


def matrix_rank(matrix):
    return rref(matrix)[0]


def solve_square(matrix, rhs):
    n = len(matrix)
    assert len(matrix[0]) == n and len(rhs) == n
    rank, pivots, red = rref([matrix[i] + [Q(rhs[i])] for i in range(n)])
    assert rank == n and pivots[:n] == list(range(n))
    return [red[i][-1] for i in range(n)]


def columns_to_matrix(columns, nrows):
    return [[pad(columns[j], nrows)[i] for j in range(len(columns))]
            for i in range(nrows)]


def primitive_normal(columns):
    equations = [[column[i] for i in range(7)] for column in columns]
    rank, pivots, red = rref(equations)
    assert rank == 6
    free = next(j for j in range(7) if j not in pivots)
    v = [Q(0)] * 7
    v[free] = 1
    for i, col in enumerate(pivots):
        v[col] = -red[i][free]
    den = reduce(lcm, (x.denominator for x in v), 1)
    nums = [int(x * den) for x in v]
    content = reduce(gcd, (abs(x) for x in nums if x), 0)
    return [x // content for x in nums]


def norm5_auxiliary_jacobian(f, q, a, b, kappa):
    """Columns in gauge B(0)=1: dq0,dq1,dA0..dA5,dB1,dB2,dk."""
    cols = []
    for dq in ([1], [0, 1]):
        cols.append(scale(mul(power(q, 4), dq), -5 * kappa))
    for i in range(6):
        cols.append(scale(mul(a, [Q(0)] * i + [Q(1)]), 2))
    for db in ([0, 1], [0, 0, 1]):
        cols.append(scale(mul(mul(b, f), db), -2))
    cols.append(scale(power(q, 5), -1))
    return columns_to_matrix(cols, 11)


def contact3_auxiliary_jacobian(q, h, kappa):
    cols = []
    for dq in ([1], [0, 1]):
        cols.append(scale(mul(power(q, 2), dq), -3 * kappa))
    for i in range(4):
        cols.append(scale(mul(h, [Q(0)] * i + [Q(1)]), 2))
    cols.append(scale(power(q, 3), -1))
    return columns_to_matrix(cols, 7)


def halving4_auxiliary_jacobian(q2, u4, line, kappa):
    ell = mul(q2, line)
    cols = []
    for dq in ([1], [0, 1]):
        dell = mul(dq, line)
        cols.append(add(scale(mul(ell, dell), 2),
                        scale(mul(power(u4, 2), dq), -kappa)))
    for du in ([1], [0, 1]):
        cols.append(scale(mul(mul(u4, du), q2), -2 * kappa))
    for dl in ([1], [0, 1]):
        cols.append(scale(mul(ell, mul(q2, dl)), 2))
    cols.append(scale(mul(power(u4, 2), q2), -1))
    return columns_to_matrix(cols, 7)


def main():
    f = [Q(94277468160), 0, Q(-22332312000), 0,
         Q(1761500625), 0, Q(-46250000)]

    # D5=12*P60, with Mumford data
    # q5=x^2-4608/395 and v5=-(164280/79)*x.
    q5 = [Q(-4608, 395), 0, 1]
    v5 = [0, Q(-164280, 79)]
    a5 = [0, Q(66600), 0, Q(-1558625, 144), 0,
          Q(48752125, 110592)]
    b5 = [Q(1), 0, Q(-125, 1728)]
    k5 = Q(1778923230671875, 4076863488)

    # D3=20*P60, with Mumford data q3=x^2-316/25, v3=-5476.
    q3 = [Q(-316, 25), 0, 1]
    v3 = [Q(-5476)]
    h3 = [Q(29600), 0, Q(-2775), 0]
    k3 = Q(46250000)

    # D4=15*P60 halves D2=30*P60.
    q2 = [Q(-1728, 125), 0, 1]
    u4 = [Q(-1506, 125), 0, 1]
    v4 = [Q(-24642, 5)]
    line = [Q(2775), 0]
    ell = mul(q2, line)
    k4 = Q(46250000)

    assert sub(sub(power(a5, 2), mul(power(b5, 2), f)),
               scale(power(q5, 5), k5)) == [0]
    assert divrem(scale(a5, -1), q5)[1] == mul(divrem(b5, q5)[1], v5)
    assert sub(sub(power(h3, 2), f), scale(power(q3, 3), k3)) == [0]
    assert divrem(h3, q3)[1] == v3
    assert sub(sub(power(ell, 2), f),
               scale(mul(power(u4, 2), q2), k4)) == [0]
    assert divrem(ell, u4)[1] == v4

    print("HLP_Z60_SIMULTANEOUS_DEFORMATION_TANGENT")
    print("EXACT_IDENTITIES verified: order5 norm, order3 contact, order4 half")

    j5 = norm5_auxiliary_jacobian(f, q5, a5, b5, k5)
    j3 = contact3_auxiliary_jacobian(q3, h3, k3)
    j4 = halving4_auxiliary_jacobian(q2, u4, line, k4)
    print("AUXILIARY_BLOCK_RANKS", matrix_rank(j5), matrix_rank(j3),
          matrix_rank(j4), "expected 11 7 7")

    # Variables: F(7), norm5(11), contact3(7), halving4(7).
    full = []
    for i in range(11):
        row = [Q(0)] * 32
        # dE/dF_j = -B^2*x^j.
        b2 = power(b5, 2)
        for j in range(7):
            if 0 <= i - j < len(b2):
                row[j] = -b2[i - j]
        row[7:18] = j5[i]
        full.append(row)
    for offset, block in ((18, j3), (25, j4)):
        for i in range(7):
            row = [Q(0)] * 32
            row[i] = -1
            row[offset:offset + 7] = block[i]
            full.append(row)
    rank = matrix_rank(full)
    print("FULL_JACOBIAN rank", rank, "variables 32 equations 25",
          "tangent_dimension", 32 - rank)
    print("PROJECTION_TO_SEXTIC_TANGENT rank 7")

    # The exact algebraic line in sextic space F_t=F+t*(1+x) cuts out a
    # smooth one-dimensional incidence slice through the seed.
    slice_rows = []
    for entries in (((0, 1), (1, -1)), ((2, 1),), ((3, 1),),
                    ((4, 1),), ((5, 1),), ((6, 1),)):
        row = [Q(0)] * 32
        for index, value in entries:
            row[index] = value
        slice_rows.append(row)
    slice_rank = matrix_rank(full + slice_rows)
    print("TRANSVERSE_SLICE F_t=F+t*(1+x) equations 31 rank", slice_rank,
          "local_dimension", 32 - slice_rank)

    # Tangent to the Humbert-4 branch through the visible involution x -> -x.
    translation = [Q(0)] * 7
    scaling = [Q(0)] * 7
    projective = [Q(0)] * 7
    for i in range(1, 7):
        translation[i - 1] = i * f[i]
        scaling[i] = i * f[i]
    for i in range(7):
        if i + 1 <= 6:
            projective[i + 1] += (i - 6) * f[i]
    even = []
    for i in (0, 2, 4, 6):
        e = [Q(0)] * 7
        e[i] = 1
        even.append(e)
    normal = primitive_normal(even + [translation, scaling, projective])
    print("HUMBERT4_NORMAL df0..df6", normal)

    df = [Q(1), Q(1), 0, 0, 0, 0, 0]
    pairing = sum(normal[i] * df[i] for i in range(7))
    print("TRANSVERSE_DIRECTION dF=1+x normal_pairing", pairing)
    assert pairing != 0

    # Unique rational lifts of the chosen curve tangent to all three blocks.
    rhs5 = pad(mul(power(b5, 2), df), 11)
    lift5 = solve_square(j5, rhs5)
    lift3 = solve_square(j3, df)
    lift4 = solve_square(j4, df)
    names5 = ["dq50", "dq51"] + [f"dA{i}" for i in range(6)] + \
             ["dB1", "dB2", "dk5"]
    names3 = ["dq30", "dq31"] + [f"dH{i}" for i in range(4)] + ["dk3"]
    names4 = ["dq20", "dq21", "du40", "du41", "dL0", "dL1", "dk4"]
    print("RATIONAL_TANGENT_LIFT dF=1+x")
    print("ORDER5", ", ".join(f"{n}={v}" for n, v in zip(names5, lift5)))
    print("ORDER3", ", ".join(f"{n}={v}" for n, v in zip(names3, lift3)))
    print("ORDER4", ", ".join(f"{n}={v}" for n, v in zip(names4, lift4)))
    print("HLP_Z60_SIMULTANEOUS_DEFORMATION_TANGENT_DONE")


if __name__ == "__main__":
    main()
