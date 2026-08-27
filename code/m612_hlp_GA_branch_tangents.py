#!/usr/bin/env python3
"""Exact first-order lifts of the three marked branches on the G_A slice.

The auxiliary Jacobians at the HLP seed are invertible, so these rational
solutions certify that t is an etale local parameter on each selected
branch.  They also distinguish the A and B points on the common degree-40
contact-support cover and the G point on the degree-120 halving cover.
"""

from fractions import Fraction as Q

from m612_hlp_deformation_tangent import (
    contact_auxiliary_jacobian,
    halving_auxiliary_jacobian,
    solve_square,
)


def show(label, names, values):
    print(label, ", ".join(f"{name}={value}" for name, value in zip(names, values)))


def main():
    q_a = [Q(1), Q(-61, 8), Q(1)]
    h_a = [Q(-3904, 9), Q(61, 3), Q(-61, 3), Q(3904, 9)]
    k_a = Q(62464, 81)

    q_b = [Q(-13, 48), 0, Q(1)]
    h_b = [Q(2623, 6), 0, Q(-183), 0]
    k_b = Q(-187392)

    q0 = [Q(1), 0, Q(1)]
    u = [Q(-32, 29), 0, Q(1)]
    line = [0, Q(183)]
    k_g = Q(-153903)

    ja = contact_auxiliary_jacobian(q_a, h_a, k_a)
    jb = contact_auxiliary_jacobian(q_b, h_b, k_b)
    jg = halving_auxiliary_jacobian(q0, u, line, k_g)

    # G_A=2+x-x^2+x^3+x^4+x^5+x^6, low coefficients first.
    df = [Q(2), Q(1), Q(-1), Q(1), Q(1), Q(1), Q(1)]
    sa = solve_square(ja, df)
    sb = solve_square(jb, df)
    sg = solve_square(jg, df)

    print("M612_HLP_GA_MARKED_BRANCH_TANGENTS")
    print("dF", df)
    names_a = ["dV", "dU", "dh0", "dh1", "dh2", "dh3", "dk"]
    names_g = ["dc0", "dc1", "du0", "du1", "dl0", "dl1", "dk"]
    show("contact_A", names_a, sa)
    show("contact_B", names_a, sb)
    show("halving_G", names_g, sg)

    # Coordinates used by the eliminant charts H^2=w*(1+r1*x+...)^2 and
    # L^2=wG*(x+z)^2.
    for label, h, sol in (("A", h_a, sa), ("B", h_b, sb)):
        h0 = h[0]
        dh0, dh1, dh2, dh3 = sol[2:6]
        dr = [
            (dhi * h0 - hi * dh0) / (h0 * h0)
            for hi, dhi in zip(h[1:], (dh1, dh2, dh3))
        ]
        print(
            f"normalized_contact_{label}",
            f"dU={sol[1]}",
            f"dV={sol[0]}",
            f"dw={2*h0*dh0}",
            f"dr1={dr[0]}",
            f"dr2={dr[1]}",
            f"dr3={dr[2]}",
        )

    dz = sg[4] / line[1]
    dw_g = 2 * line[1] * sg[5]
    print(
        "normalized_halving_G",
        f"dc0={sg[0]}",
        f"dc1={sg[1]}",
        f"du0={sg[2]}",
        f"du1={sg[3]}",
        f"dz={dz}",
        f"dw={dw_g}",
    )


if __name__ == "__main__":
    main()
