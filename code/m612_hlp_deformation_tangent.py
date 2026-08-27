#!/usr/bin/env python3
"""Exact marked-torsion deformation and tangent calculation at the HLP seed.

This uses only Python's Fraction class.  It reconstructs the two cubic-contact
identities and the 2-halving identity from the Mumford data, forms the general
21-equation incidence system in 28 variables, evaluates its Jacobian, and
tests a rational tangent direction against the three elliptic-involution
branches visible at the even reciprocal seed.
"""

from fractions import Fraction as Q
from functools import reduce
from math import gcd, lcm


def pad(a, n=7):
    return list(a) + [Q(0)] * (n - len(a))


def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def add(a, b):
    n = max(len(a), len(b))
    return trim(
        [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
         for i in range(n)]
    )


def sub(a, b):
    n = max(len(a), len(b))
    return trim(
        [(a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)
         for i in range(n)]
    )


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


def polynomial_text(a):
    return "[" + ", ".join(str(x) for x in a) + "] (low coefficients)"


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


def matrix_rank(matrix):
    return rref(matrix)[0]


def solve_square(matrix, rhs):
    n = len(matrix)
    aug = [matrix[i][:] + [Q(rhs[i])] for i in range(n)]
    rank, pivots, red = rref(aug)
    assert rank == n and pivots[:n] == list(range(n))
    return [red[i][-1] for i in range(n)]


def columns_to_matrix(columns):
    return [[pad(columns[j])[i] for j in range(len(columns))]
            for i in range(7)]


def contact_auxiliary_jacobian(q, h, kappa):
    columns = []
    for dq in ([1], [0, 1]):
        columns.append(scale(mul(power(q, 2), dq), -3 * kappa))
    for dh in ([1], [0, 1], [0, 0, 1], [0, 0, 0, 1]):
        columns.append(scale(mul(h, dh), 2))
    columns.append(scale(power(q, 3), -1))
    return columns_to_matrix(columns)


def halving_auxiliary_jacobian(q0, u, line, kappa):
    ell = mul(q0, line)
    columns = []
    for dq in ([1], [0, 1]):
        dell = mul(dq, line)
        columns.append(add(scale(mul(ell, dell), 2),
                           scale(mul(power(u, 2), dq), -kappa)))
    for du in ([1], [0, 1]):
        columns.append(scale(mul(mul(u, du), q0), -2 * kappa))
    for dl in ([1], [0, 1]):
        columns.append(scale(mul(ell, mul(q0, dl)), 2))
    columns.append(scale(mul(power(u, 2), q0), -1))
    return columns_to_matrix(columns)


def primitive_normal(columns):
    # Find the one-dimensional orthogonal complement of a rank-6 subspace.
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


def main():
    # All polynomials use low-to-high coefficients.
    f = [Q(187392), 0, Q(-118767), 0, Q(-118767), 0, Q(187392)]

    q_a = [Q(1), Q(-61, 8), Q(1)]
    v_a = [Q(-3721), Q(197213, 8)]
    h_a = [Q(-3904, 9), Q(61, 3), Q(-61, 3), Q(3904, 9)]
    k_a = Q(62464, 81)

    q_b = [Q(-13, 48), 0, Q(1)]
    v_b = [Q(18605, 48)]
    h_b = [Q(2623, 6), 0, Q(-183), 0]
    k_b = Q(-187392)

    q0 = [Q(1), 0, Q(1)]
    u = [Q(-32, 29), 0, Q(1)]
    v_g = [0, Q(11163, 29)]
    line = [0, Q(183)]
    ell = mul(q0, line)
    k_g = Q(-153903)

    assert sub(sub(power(h_a, 2), f), scale(power(q_a, 3), k_a)) == [0]
    assert sub(sub(power(h_b, 2), f), scale(power(q_b, 3), k_b)) == [0]
    assert sub(sub(power(ell, 2), f),
               scale(mul(power(u, 2), q0), k_g)) == [0]
    assert divrem(h_a, q_a)[1] == v_a
    assert divrem(h_b, q_b)[1] == v_b
    assert divrem(ell, u)[1] == v_g

    print("EXACT MARKED IDENTITIES")
    print("f   =", polynomial_text(f))
    print("qA  =", polynomial_text(q_a))
    print("HA  =", polynomial_text(h_a), "kA =", k_a)
    print("qB  =", polynomial_text(q_b))
    print("HB  =", polynomial_text(h_b), "kB =", k_b)
    print("q0  =", polynomial_text(q0))
    print("uG  =", polynomial_text(u))
    print("ell =", polynomial_text(ell), "kG =", k_g)

    ja = contact_auxiliary_jacobian(q_a, h_a, k_a)
    jb = contact_auxiliary_jacobian(q_b, h_b, k_b)
    jg = halving_auxiliary_jacobian(q0, u, line, k_g)
    print("AUXILIARY BLOCK RANKS", matrix_rank(ja), matrix_rank(jb),
          matrix_rank(jg))

    # Variables: f0..f6, A(7), B(7), G(7).  There are 21 identities.
    full = []
    for offset, block in ((7, ja), (14, jb), (21, jg)):
        for i in range(7):
            row = [Q(0)] * 28
            row[i] = -1
            row[offset:offset + 7] = block[i]
            full.append(row)
    full_rank = matrix_rank(full)
    print("FULL JACOBIAN rank", full_rank, "variables 28 equations 21",
          "tangent_dimension", 28 - full_rank)
    print("projection_to_f_rank 7 (each auxiliary block is invertible)")

    slice_rows = []
    for entries in (((0, 1), (1, -1)), ((2, 1),), ((3, 1),),
                    ((4, 1),), ((5, 1),), ((6, 1),)):
        row = [Q(0)] * 28
        for index, value in entries:
            row[index] = value
        slice_rows.append(row)
    slice_rank = matrix_rank(full + slice_rows)
    print("TRANSVERSE SLICE F_t=F_seed+t*(1+x): equations 27",
          "Jacobian_rank", slice_rank, "local_dimension", 28 - slice_rank)

    # Visible reduced involutions: -x, 1/x, -1/x.  Their marked Humbert-4
    # branch tangents equal the invariant binary-sextic subspace plus the
    # infinitesimal PGL2 orbit of f.
    translation = [Q(0)] * 7
    scaling = [Q(0)] * 7
    projective = [Q(0)] * 7
    for i in range(1, 7):
        translation[i - 1] = i * f[i]
        scaling[i] = i * f[i]
        if i + 1 <= 6:
            projective[i + 1] += i * f[i]
    for i in range(7):
        if i + 1 <= 6:
            projective[i + 1] -= 6 * f[i]
    orbit = [translation, scaling, projective]

    even = []
    for i in (0, 2, 4, 6):
        e = [Q(0)] * 7
        e[i] = 1
        even.append(e)

    reciprocal = []
    anti_reciprocal = []
    for i in range(4):
        e = [Q(0)] * 7
        e[i] = 1
        if 6 - i != i:
            e[6 - i] = 1
        reciprocal.append(e)
        e = [Q(0)] * 7
        e[i] = 1
        if 6 - i != i:
            e[6 - i] = (-1) ** i
        anti_reciprocal.append(e)

    normals = {
        "x->-x": primitive_normal(even + orbit),
        "x->1/x": primitive_normal(reciprocal + orbit),
        "x->-1/x": primitive_normal(anti_reciprocal + orbit),
    }
    print("HUMBERT-4 BRANCH NORMALS (coefficients on df0..df6)")
    for name, normal in normals.items():
        print(name, normal)

    # A simple rational tangent outside all three branch tangent hyperplanes.
    df = [Q(1), Q(1), 0, 0, 0, 0, 0]
    print("TRANSVERSE CURVE TANGENT df = 1+x")
    for name, normal in normals.items():
        print(name, "normal(df) =", sum(normal[i] * df[i] for i in range(7)))

    names = {
        "A": ["dqa0", "dqa1", "dhA0", "dhA1", "dhA2", "dhA3", "dkA"],
        "B": ["dqb0", "dqb1", "dhB0", "dhB1", "dhB2", "dhB3", "dkB"],
        "G": ["dq00", "dq01", "du0", "du1", "dl0", "dl1", "dkG"],
    }
    print("RATIONAL TANGENT LIFT FOR df=1+x")
    for label, block in (("A", ja), ("B", jb), ("G", jg)):
        solution = solve_square(block, df)
        print(label, ", ".join(f"{n}={v}" for n, v in zip(names[label], solution)))


if __name__ == "__main__":
    main()
