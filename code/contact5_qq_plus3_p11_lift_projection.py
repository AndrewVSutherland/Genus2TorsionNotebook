#!/usr/bin/env python3
"""Bounded 11-adic lift projection for the live qq_tminus3 branch.

This follows the two boundary cubic-contact points

    (m,U,V) = (2,10,0), (9,10,0) mod 11

which the tangent diagnostic identifies as the points capable of moving in
every r direction.  It enumerates their lifts modulo 11^2, then uses the
linearized lift equation to enumerate all children modulo 11^3.  The output
records only projection counts on the r-line, not a claim about all higher
11-adic levels.
"""

from collections import defaultdict

P = 11
LIVE_R = (3, 5, 7, 8)
OPEN_CONTACT_POINTS = ((2, 10, 0), (9, 10, 0))


def inv(a, modulus):
    return pow(a % modulus, -1, modulus)


def t_qq(r, modulus):
    den = (r * r - 2) ** 2 * (r * r - 2 * r + 2) % modulus
    num = (
        r**6
        - 2 * r**5
        + 2 * r**4
        - 4 * r**3
        + 4 * r * r
        - 8 * r
        + 8
    ) % modulus
    return -num * inv(den, modulus) % modulus


def family_coefficients(t, modulus):
    inv2 = inv(2, modulus)
    inv4 = inv(4, modulus)
    b = (t * t - 1) * inv2 % modulus
    return (
        1 % modulus,
        2 * t % modulus,
        (t * t + 2 * b) % modulus,
        2 * t * b % modulus,
        b * b % modulus,
        -(t + 1) ** 4 * inv4 % modulus,
    )


def contact_equations(point, modulus):
    r, m, u, v = point
    t = t_qq(r, modulus)
    coeffs = family_coefficients(t, modulus)
    inv2m = inv(2 * m, modulus)
    ncoef = (3 * m * m * u + coeffs[5]) * inv2m % modulus
    rcoef = (
        3 * m * m * (u * u + v) + coeffs[4] - ncoef * ncoef
    ) * inv2m % modulus
    scoef = (
        m * m * (u**3 + 6 * u * v)
        + coeffs[3]
        - 2 * ncoef * rcoef
    ) * inv2m % modulus
    return (
        (
            rcoef * rcoef
            + 2 * ncoef * scoef
            - coeffs[2]
            - 3 * m * m * (u * u * v + v * v)
        )
        % modulus,
        (2 * rcoef * scoef - coeffs[1] - 3 * m * m * u * v * v)
        % modulus,
        (scoef * scoef - coeffs[0] - m * m * v**3) % modulus,
    )


def jacobian_columns(base):
    """Jacobian mod 11 from first differences evaluated mod 11^2."""
    modulus = P * P
    base_value = contact_equations(base, modulus)
    columns = []
    for j in range(4):
        moved = list(base)
        moved[j] += P
        moved_value = contact_equations(tuple(moved), modulus)
        columns.append(
            tuple(
                ((moved_value[i] - base_value[i]) // P) % P
                for i in range(3)
            )
        )
    return columns


def direction_map(columns):
    result = defaultdict(list)
    for dr in range(P):
        for dm in range(P):
            for du in range(P):
                for dv in range(P):
                    direction = (dr, dm, du, dv)
                    image = tuple(
                        sum(columns[j][i] * direction[j] for j in range(4)) % P
                        for i in range(3)
                    )
                    result[image].append(direction)
    return result


def lifts_mod_p2(r0, contact_point):
    modulus = P * P
    m0, u0, v0 = contact_point
    lifts = []
    for dr in range(P):
        for dm in range(P):
            for du in range(P):
                for dv in range(P):
                    point = (
                        r0 + P * dr,
                        m0 + P * dm,
                        u0 + P * du,
                        v0 + P * dv,
                    )
                    if contact_equations(point, modulus) == (0, 0, 0):
                        lifts.append(point)
    return lifts


def main():
    print("CONTACT5_QQ_P11_LIFT_PROJECTION_START")
    for r0 in LIVE_R:
        p2_points = []
        p3_projection = defaultdict(int)
        for contact_point in OPEN_CONTACT_POINTS:
            base = (r0,) + contact_point
            assert contact_equations(base, P) == (0, 0, 0)
            columns = jacobian_columns(base)
            images = direction_map(columns)
            points = lifts_mod_p2(r0, contact_point)
            p2_points.extend(points)
            for point in points:
                values = contact_equations(point, P**3)
                rhs = tuple(-(value // (P * P)) % P for value in values)
                for direction in images.get(rhs, ()):
                    r3 = point[0] + P * P * direction[0]
                    p3_projection[r3] += 1

        p2_projection = {point[0] for point in p2_points}
        expected_p2 = {r0 + P * a for a in range(P)}
        expected_p3 = {
            r0 + P * a + P * P * b for a in range(P) for b in range(P)
        }
        assert p2_projection == expected_p2
        assert set(p3_projection) == expected_p3
        multiplicities = sorted(set(p3_projection.values()))
        print(
            "RESIDUE",
            r0,
            "p2_points",
            len(p2_points),
            "p2_r_classes",
            len(p2_projection),
            "p3_r_classes",
            len(p3_projection),
            "of",
            P * P,
            "p3_projection_multiplicities",
            multiplicities,
        )
    print("VERDICT no_parameter_congruence_gain_through_mod_11^3")
    print("CONTACT5_QQ_P11_LIFT_PROJECTION_DONE")


if __name__ == "__main__":
    main()
