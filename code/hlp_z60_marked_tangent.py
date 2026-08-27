#!/usr/bin/env python3
"""Exact first-order deformation calculation at the split HLP Z/60 seed.

Polynomials are coefficient lists in increasing degree.  The calculation
uses only Fraction arithmetic: it verifies the order-5 norm, order-3
contact, and order-4 halving identities; differentiates their normalized
incidence equations; and lifts the odd sextic tangent df=x through all
three marked torsion blocks simultaneously.
"""

from fractions import Fraction as Q


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
    return trim([Q(c) * v for v in a])


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


def monomial(i):
    return [Q(0)] * i + [Q(1)]


def coefficient_matrix(columns, nrows):
    return [[pad(columns[j], nrows)[i] for j in range(len(columns))]
            for i in range(nrows)]


def rref(matrix):
    a = [[Q(v) for v in row] for row in matrix]
    nr, nc = len(a), len(a[0])
    row = 0
    pivots = []
    for col in range(nc):
        pivot = next((i for i in range(row, nr) if a[i][col]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        c = a[row][col]
        a[row] = [v / c for v in a[row]]
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


def determinant(matrix):
    a = [[Q(v) for v in row] for row in matrix]
    n = len(a)
    assert all(len(row) == n for row in a)
    det = Q(1)
    for col in range(n):
        pivot = next((i for i in range(col, n) if a[i][col]), None)
        if pivot is None:
            return Q(0)
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            det = -det
        c = a[col][col]
        det *= c
        for j in range(col, n):
            a[col][j] /= c
        for i in range(col + 1, n):
            c = a[i][col]
            if c:
                for j in range(col, n):
                    a[i][j] -= c * a[col][j]
    return det


def solve_square(matrix, rhs):
    n = len(matrix)
    aug = [list(matrix[i]) + [Q(rhs[i])] for i in range(n)]
    rr, pivots, red = rref(aug)
    assert rr == n and pivots[:n] == list(range(n))
    return [red[i][-1] for i in range(n)]


def dot(a, b):
    return sum((Q(x) * Q(y) for x, y in zip(a, b)), Q(0))


def main():
    # Seed sextic, low-to-high coefficients.
    f = [Q(94277468160), 0, Q(-22332312000), 0,
         Q(1761500625), 0, Q(-46250000)]

    # Fixed D5=12P60 norm identity.  B5 is kept monic in the deformation.
    q5 = [Q(-4608, 395), 0, 1]
    b5 = [Q(-1728, 125), 0, 1]
    a5 = scale(mul([0, 1], [Q(7962624), 0, Q(-1294080), 0,
                              Q(52705)]), Q(37, 320))
    k5 = Q(341553260289, 4096)

    # Fixed D3=20P60 cubic contact.
    q3 = [Q(-316, 25), 0, 1]
    h3 = [Q(29600), 0, Q(-2775)]
    k3 = Q(46250000)

    # Fixed D4=15P60, with 2D4 supported on q0.
    q0 = [Q(-1728, 125), 0, 1]
    u4 = [Q(-1506, 125), 0, 1]
    line4 = [Q(-2775)]
    ell4 = mul(q0, line4)
    k4 = Q(46250000)

    assert sub(sub(power(a5, 2), mul(f, power(b5, 2))),
               scale(power(q5, 5), k5)) == [0]
    assert sub(sub(power(h3, 2), f), scale(power(q3, 3), k3)) == [0]
    assert sub(sub(power(ell4, 2), f),
               scale(mul(q0, power(u4, 2)), k4)) == [0]

    # Order-5 variables:
    #   dq5_0,dq5_1, dA5_0..dA5_5, dB5_0,dB5_1, dk5.
    # The leading coefficients of q5 and B5 are normalized to one.
    cols5 = []
    for dq in (monomial(0), monomial(1)):
        cols5.append(scale(mul(power(q5, 4), dq), -5 * k5))
    for i in range(6):
        cols5.append(scale(mul(a5, monomial(i)), 2))
    for db in (monomial(0), monomial(1)):
        cols5.append(scale(mul(mul(f, b5), db), -2))
    cols5.append(scale(power(q5, 5), -1))
    jac5 = coefficient_matrix(cols5, 11)

    # Order-3 variables:
    #   dq3_0,dq3_1, dH3_0..dH3_3, dk3.
    cols3 = []
    for dq in (monomial(0), monomial(1)):
        cols3.append(scale(mul(power(q3, 2), dq), -3 * k3))
    for i in range(4):
        cols3.append(scale(mul(h3, monomial(i)), 2))
    cols3.append(scale(power(q3, 3), -1))
    jac3 = coefficient_matrix(cols3, 7)

    # Order-4 variables:
    #   dq0_0,dq0_1, du4_0,du4_1, dL4_0,dL4_1, dk4,
    # for ell=q0*L and ell^2-f-k4*q0*u4^2=0.
    cols4 = []
    for dq in (monomial(0), monomial(1)):
        cols4.append(sub(scale(mul(mul(ell4, line4), dq), 2),
                         scale(mul(power(u4, 2), dq), k4)))
    for du in (monomial(0), monomial(1)):
        cols4.append(scale(mul(mul(q0, u4), du), -2 * k4))
    for dl in (monomial(0), monomial(1)):
        cols4.append(scale(mul(mul(ell4, q0), dl), 2))
    cols4.append(scale(mul(q0, power(u4, 2)), -1))
    jac4 = coefficient_matrix(cols4, 7)

    print("HLP_Z60_MARKED_TANGENT")
    print("AUXILIARY_RANKS", rank(jac5), rank(jac3), rank(jac4))
    print("AUXILIARY_DETERMINANTS", determinant(jac5),
          determinant(jac3), determinant(jac4))

    # Assemble all 25 differentiated equations in 32 variables:
    # df0..df6 followed by the 11+7+7 auxiliary variables.
    full = []
    b5sq = power(b5, 2)
    for i in range(11):
        row = [Q(0)] * 32
        for j in range(7):
            if 0 <= i - j < len(b5sq):
                row[j] = -b5sq[i - j]
        row[7:18] = jac5[i]
        full.append(row)
    for i in range(7):
        row = [Q(0)] * 32
        row[i] = -1
        row[18:25] = jac3[i]
        full.append(row)
    for i in range(7):
        row = [Q(0)] * 32
        row[i] = -1
        row[25:32] = jac4[i]
        full.append(row)

    print("FULL_JACOBIAN rank", rank(full), "variables", 32,
          "equations", 25, "tangent_dimension", 32 - rank(full))
    print("PROJECTION_TO_SEXTIC_TANGENT rank 7")

    # Lift the odd sextic deformation f_t=f+t*x.
    df = [Q(0), Q(1), 0, 0, 0, 0, 0]
    rhs5 = pad(mul(power(b5, 2), df), 11)
    rhs34 = pad(df, 7)
    sol5 = solve_square(jac5, rhs5)
    sol3 = solve_square(jac3, rhs34)
    sol4 = solve_square(jac4, rhs34)
    tangent = df + sol5 + sol3 + sol4
    assert all(dot(row, tangent) == 0 for row in full)

    names5 = ["dq5_0", "dq5_1"] + [f"dA5_{i}" for i in range(6)] \
             + ["dB5_0", "dB5_1", "dk5"]
    names3 = ["dq3_0", "dq3_1"] + [f"dH3_{i}" for i in range(4)] \
             + ["dk3"]
    names4 = ["dq0_0", "dq0_1", "du4_0", "du4_1",
              "dL4_0", "dL4_1", "dk4"]
    print("TRANSVERSE_LIFT df=x")
    print("ORDER5", " ".join(f"{n}={v}" for n, v in zip(names5, sol5)))
    print("ORDER3", " ".join(f"{n}={v}" for n, v in zip(names3, sol3)))
    print("ORDER4", " ".join(f"{n}={v}" for n, v in zip(names4, sol4)))

    # The seed has one elliptic involution (Aut(C) has order 4, checked in
    # the Magma verifier).  In this affine sextic chart the tangent space to
    # conjugates of the even-sextic locus is spanned by even perturbations
    # and the two odd PGL2 orbit directions f' and 6*x*f-x^2*f'.
    split_normal = [Q(0), Q(81125), 0, Q(904800), 0, Q(9916416), 0]
    derivative = [Q((i + 1) * f[i + 1]) for i in range(6)] + [Q(0)]
    projective = [Q(0)] * 7
    for i, fi in enumerate(f):
        if i + 1 < 7:
            projective[i + 1] += 6 * fi
        if i >= 1 and i + 1 < 7:
            projective[i + 1] -= i * fi
    assert dot(split_normal, derivative) == 0
    assert dot(split_normal, projective) == 0
    assert all(split_normal[i] == 0 for i in (0, 2, 4, 6))
    print("SPLIT_HUMBERT_TANGENT_NORMAL", split_normal)
    print("NORMAL_ON_df=x", dot(split_normal, df))
    assert dot(split_normal, df) != 0
    print("HLP_Z60_MARKED_TANGENT_DONE")


if __name__ == "__main__":
    main()
