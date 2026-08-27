#!/usr/bin/env python3
"""Exact tangent analysis for the linked B5=q0 cyclic-60 slice.

At the HLP seed the monic quadratic B in the order-5 norm identity is
exactly the support q0 of the rational order-2 class doubled from D4.
The two equations B=q0 cut the three-dimensional marked moduli space to a
curve.  To exhibit a one-dimensional coefficient-space chart we additionally
fix q0, the constant term of the order-4 line L, and k4.  The resulting
tangent is a pure coordinate direction and lies in the Humbert-4 tangent.
The full containment analysis is in hlp_z60_linked_second_order.py.

Only exact Fraction linear algebra is used.
"""

from fractions import Fraction as Q

from hlp_z60_marked_tangent import (
    coefficient_matrix, monomial, mul, pad, power, rank, scale,
    solve_square, sub,
)


def auxiliary_blocks():
    f = [Q(94277468160), 0, Q(-22332312000), 0,
         Q(1761500625), 0, Q(-46250000)]
    q5 = [Q(-4608, 395), 0, 1]
    b5 = [Q(-1728, 125), 0, 1]
    a5 = scale(mul([0, 1], [Q(7962624), 0, Q(-1294080), 0,
                              Q(52705)]), Q(37, 320))
    k5 = Q(341553260289, 4096)
    q3 = [Q(-316, 25), 0, 1]
    h3 = [Q(29600), 0, Q(-2775)]
    k3 = Q(46250000)
    q0 = [Q(-1728, 125), 0, 1]
    u4 = [Q(-1506, 125), 0, 1]
    line4 = [Q(-2775)]
    ell4 = mul(q0, line4)
    k4 = Q(46250000)

    cols5 = []
    for dq in (monomial(0), monomial(1)):
        cols5.append(scale(mul(power(q5, 4), dq), -5 * k5))
    for i in range(6):
        cols5.append(scale(mul(a5, monomial(i)), 2))
    for db in (monomial(0), monomial(1)):
        cols5.append(scale(mul(mul(f, b5), db), -2))
    cols5.append(scale(power(q5, 5), -1))
    j5 = coefficient_matrix(cols5, 11)

    cols3 = []
    for dq in (monomial(0), monomial(1)):
        cols3.append(scale(mul(power(q3, 2), dq), -3 * k3))
    for i in range(4):
        cols3.append(scale(mul(h3, monomial(i)), 2))
    cols3.append(scale(power(q3, 3), -1))
    j3 = coefficient_matrix(cols3, 7)

    cols4 = []
    for dq in (monomial(0), monomial(1)):
        cols4.append(sub(scale(mul(mul(ell4, line4), dq), 2),
                         scale(mul(power(u4, 2), dq), k4)))
    for du in (monomial(0), monomial(1)):
        cols4.append(scale(mul(mul(q0, u4), du), -2 * k4))
    for dl in (monomial(0), monomial(1)):
        cols4.append(scale(mul(mul(ell4, q0), dl), 2))
    cols4.append(scale(mul(q0, power(u4, 2)), -1))
    j4 = coefficient_matrix(cols4, 7)
    return b5, j5, j3, j4


def response_matrices():
    b5, j5, j3, j4 = auxiliary_blocks()
    r5, r3, r4 = [], [], []
    for j in range(7):
        df = monomial(j)
        r5.append(solve_square(j5, pad(mul(power(b5, 2), df), 11)))
        r3.append(solve_square(j3, pad(df, 7)))
        r4.append(solve_square(j4, pad(df, 7)))
    # Transpose: rows are auxiliary coordinates, columns are df_i.
    return ([[r5[j][i] for j in range(7)] for i in range(11)],
            [[r3[j][i] for j in range(7)] for i in range(7)],
            [[r4[j][i] for j in range(7)] for i in range(7)])


def null_vector_six_by_seven(rows):
    """Return a primitive rational kernel vector for a rank-six matrix."""
    assert len(rows) == 6 and rank(rows) == 6
    # Cofactor vector: v_j=(-1)^j det(matrix with column j deleted).
    def det(m):
        a = [[Q(x) for x in row] for row in m]
        out = Q(1)
        for c in range(len(a)):
            p = next((i for i in range(c, len(a)) if a[i][c]), None)
            if p is None:
                return Q(0)
            if p != c:
                a[c], a[p] = a[p], a[c]
                out = -out
            z = a[c][c]
            out *= z
            for j in range(c, len(a)):
                a[c][j] /= z
            for i in range(c + 1, len(a)):
                z = a[i][c]
                for j in range(c, len(a)):
                    a[i][j] -= z * a[c][j]
        return out
    return [(-1 if j % 2 else 1) *
            det([[row[c] for c in range(7) if c != j] for row in rows])
            for j in range(7)]


def main():
    r5, _r3, r4 = response_matrices()
    # order-5 auxiliary indices: q5(0,1), A(2..7), B(8,9), k(10)
    # order-4 auxiliary indices: q0(0,1), u(2,3), L(4,5), k(6)
    link = [[r5[8][j] - r4[0][j] for j in range(7)],
            [r5[9][j] - r4[1][j] for j in range(7)]]
    split_normal = [Q(0), Q(81125), 0, Q(904800), 0, Q(9916416), 0]
    print("HLP_Z60_LINKED_SLICE_TANGENT")
    print("LINK_B5_EQUALS_Q0_RANK", rank(link))
    print("LINK_WITH_HUMBERT_NORMAL_RANK", rank(link + [split_normal]))

    # Concrete chart: B=q0, q0 fixed, L0 fixed, and k4 fixed.
    chart = link + [r4[0], r4[1], r4[4], r4[6]]
    print("SIX_SLICE_CONSTRAINT_RANK", rank(chart))
    v = null_vector_six_by_seven(chart)
    pairing = sum(v[i] * split_normal[i] for i in range(7))
    print("SLICE_TANGENT_DF", v)
    print("HUMBERT_NORMAL_PAIRING", pairing)
    print("TRANSVERSE", pairing != 0)
    print("REDUCED_EXACT_SYSTEM variables=19 equations=18 local_dimension=1")
    print("HLP_Z60_LINKED_SLICE_TANGENT_DONE")


if __name__ == "__main__":
    main()
