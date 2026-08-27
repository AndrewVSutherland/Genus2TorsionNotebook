#!/usr/bin/env sage
"""Bounded boundary audit for the contact-9 [3,18] and contact-6 [6,12] routes.

This does two deliberately local computations:

* resolve the first visible p=7 boundary charts of the one-parameter
  contact-9 family and apply a Neron special-fibre 3-primary capacity test;
* classify every rank-five point of the seven-variable contact-6
  intersection over F_5 and enumerate its complete affine space of lifts
  modulo 25, recording which boundary initial forms can be left on the
  first blowup.

It is not a global resolution of either parameter space.
"""

from collections import Counter, defaultdict
from itertools import product

from sage.all import (
    GF,
    HyperellipticCurve,
    Matrix,
    PolynomialRing,
    QQ,
    ZZ,
    matrix,
    vector,
)


def vp(value, p):
    value = QQ(value)
    if value == 0:
        return +Infinity
    return ZZ(value.numerator()).valuation(p) - ZZ(value.denominator()).valuation(p)


def contact9_family():
    parameter_ring = PolynomialRing(QQ, "t")
    t0 = parameter_ring.gen()
    function_field = parameter_ring.fraction_field()
    t = function_field(t0)
    polynomial_ring = PolynomialRing(function_field, "x")
    x = polynomial_ring.gen()
    root = 1 - t**2
    h0 = 1 - QQ(9)/2*root + QQ(63)/8*root**2 - QQ(105)/16*root**3
    contact_a = (t**9 - h0)/root**4
    h = (
        1 - QQ(9)/2*x + QQ(63)/8*x**2 - QQ(105)/16*x**3
        + contact_a*x**4
    )
    f, remainder = (h**2 + (x - 1)**9).quo_rem(x**4)
    assert remainder == 0 and f.degree() == 5
    return parameter_ring, function_field, polynomial_ring, t, x, f


def specialize_rational_polynomial(poly, value, target_ring):
    return target_ring([coefficient(value) for coefficient in poly])


def three_part(n):
    return ZZ(n).valuation(3)


def order_at_zero(value, polynomial_ring):
    value = polynomial_ring.fraction_field()(value)
    if value == 0:
        return +Infinity
    numerator = polynomial_ring(value.numerator())
    denominator = polynomial_ring(value.denominator())
    return numerator.valuation() - denominator.valuation()


def contact9_audit():
    print("=== contact-9 [3,18] boundary audit at p=7 ===")
    print("effective parameter t=eps*s; the two signs give the same t-chart")
    parameter_ring, function_field, polynomial_ring, t, x, f = contact9_family()
    discriminant = f.discriminant()
    print("discriminant_factor", discriminant.factor())

    field = GF(7)
    finite_ring = PolynomialRing(field, "x")
    print("finite_special_fibres")
    special = {}
    for residue in range(6):
        fp = specialize_rational_polynomial(f, QQ(residue), finite_ring)
        special[residue] = fp
        print(f" t={residue} factor={fp.factor()}")

    print("first_exceptional_discriminant_valuations")
    for residue in (0, 1, 3, 4, 5):
        values = []
        for direction in range(7):
            point = QQ(residue + 7*direction)
            values.append(vp(discriminant(point), 7))
        print(f" t0={residue} directions_u=0..6 valuations={values}")

    print("ordinary_node_charts")
    for residue in (3, 4, 5):
        fp = special[residue]
        repeated = fp.gcd(fp.derivative()).monic()
        assert repeated.degree() == 1
        alpha = -repeated[0]
        normalization_cubic = (fp // repeated**2).monic()
        assert normalization_cubic.degree() == 3
        elliptic = EllipticCurve(
            field,
            [
                0,
                normalization_cubic[2],
                0,
                normalization_cubic[1],
                normalization_cubic[0],
            ],
        )
        node_value = normalization_cubic(alpha)
        split = node_value.is_square()
        torus_order = 6 if split else 8
        identity_order = elliptic.cardinality()*torus_order
        thicknesses = [
            vp(discriminant(QQ(residue + 7*direction)), 7)
            for direction in range(7)
        ]
        component_three_bound = max(
            0 if thickness is Infinity else three_part(thickness)
            for thickness in thicknesses
        )
        capacity = three_part(identity_order) + component_three_bound
        node_point = elliptic(alpha, node_value.sqrt())
        gluing_class = 2*node_point
        gluing_divisible_by_three = any(
            3*point == gluing_class for point in elliptic.points()
        )
        print(
            f" t0={residue} repeated_root={alpha} normalization={normalization_cubic} "
            f"E_invariants={elliptic.abelian_group().invariants()} "
            f"E_order={elliptic.cardinality()} split_node={split} "
            f"J0_order={identity_order} J0_3valuation={three_part(identity_order)} "
            f"thicknesses={thicknesses} component_3_bound={component_three_bound} "
            f"gluing_class_order={gluing_class.order()} "
            f"gluing_class_3divisible={gluing_divisible_by_three} "
            f"target_3valuation=3 capacity_pass={capacity >= 3}"
        )

    # At the pole t=-1, d=t+1 and z=x/d.  The integral equation is
    # Y^2=d^4*f(d*z).  Its special quartic is the normalization of a
    # one-node genus-two fibre.  The fifth branch point meets infinity at
    # distance d^9, hence the split node has thickness 18.
    d_ring = PolynomialRing(QQ, "d")
    d0 = d_ring.gen()
    d_field = d_ring.fraction_field()
    d = d_field(d0)
    z_ring = PolynomialRing(d_field, "z")
    z = z_ring.gen()
    pole_f = z_ring([d_field(coefficient(-1 + d)) for coefficient in f])
    pole_g = z_ring(d**4 * pole_f(d*z))
    quartic_ring = PolynomialRing(QQ, "z")
    quartic = quartic_ring([coefficient(0) for coefficient in pole_g])
    quartic_mod7 = PolynomialRing(field, "z")(quartic)
    elliptic_count = HyperellipticCurve(quartic_mod7).count_points(1)[0]
    assert pole_g[5] == d**9
    pole_split = quartic_mod7[4].is_square()
    pole_identity_order = elliptic_count*(6 if pole_split else 8)
    pole_thickness = 18
    pole_capacity = three_part(pole_identity_order) + three_part(pole_thickness)
    print("pole_chart_tminus1")
    print(
        f" integral_model=d^4*f(d*z) special_quartic={quartic} "
        f"factor_mod7={quartic_mod7.factor()} discriminant_mod7={quartic_mod7.discriminant()}"
    )
    print(
        f" normalization_order={elliptic_count} split_node={pole_split} "
        f"branch_collision_order=9 node_thickness={pole_thickness} "
        f"J0_3valuation={three_part(pole_identity_order)} "
        f"component_3valuation={three_part(pole_thickness)} "
        f"target_3valuation=3 capacity_pass={pole_capacity >= 3}"
    )

    # At t=infinity put d=1/t and multiply by d^2.  The four-root cluster
    # has depth 1/4, so it is not resolved by an unramified first blowup.
    infinity_f = z_ring([d_field(coefficient(1/d)) for coefficient in f])
    infinity_g = z_ring(d**2 * infinity_f)
    coefficient_orders = [order_at_zero(coefficient, d_ring) for coefficient in infinity_g]
    infinity_special = PolynomialRing(QQ, "x")([
        coefficient(0) if order_at_zero(coefficient, d_ring) >= 0 else 0
        for coefficient in infinity_g
    ])
    print("infinity_chart")
    print(f" coefficient_d_orders_after_d2_scaling={coefficient_orders}")
    print(f" special_polynomial={infinity_special}")
    print(" root_newton_slopes=-1/4 (four roots), 2 (one root); ramification_degree=4")

    print("contact9_decision")
    print(" no_go ordinary-node charts t=3 and t=5: 3-primary order at most 3^2")
    print(" continue ordinary-node chart t=4: special Neron identity has 3-primary order 3^3")
    print(" continue pole chart t=-1: split thickness-18 node supplies component 3^2")
    print(" deeper_only t=0,t=1 and t=infinity: first special fibres are not ordinary-node charts")


def contact6_symbolic_data():
    ring = PolynomialRing(ZZ, ("a", "b", "L", "U", "v", "q", "r"))
    a, b, contact_scale, contact_u, contact_v, half_q, half_r = ring.gens()
    xring = PolynomialRing(ring, "x")
    x = xring.gen()
    h6 = 1 + a*x + b*x**2 + x**3
    f = h6**2 - (x - 1)**6
    c1, c2, c3, c4, c5 = [f[i] for i in range(1, 6)]
    scale_square = contact_scale**2
    contact_b = c5*scale_square + 3*contact_u
    delta = (
        4*c4*scale_square + 12*(contact_u**2 + contact_v**2)
        - contact_b**2
    )
    f3 = (
        contact_b*delta + 16*contact_v**3 - 8*c3*scale_square
        - 8*contact_u**3 - 48*contact_u*contact_v**2
    )
    f2 = (
        delta**2 + 64*contact_b*contact_v**3 - 64*c2*scale_square
        - 192*(contact_u**2*contact_v**2 + contact_v**4)
    )
    f1 = delta*contact_v**3 - 4*c1*scale_square - 12*contact_u*contact_v**4
    center_a = b + 3
    center_c = a + b + 2
    k2 = (
        center_a**2*center_c**2*half_q
        - 2*center_a**2*center_c*half_q**2
        + center_a**2*half_q**4
        - 2*center_a*center_c**3
        - 2*center_a*center_c**2*half_q**2
        + 2*center_a*center_c**2*half_q
        + 4*center_a*center_c*half_q**2
        + 2*center_c**3*half_q
        - 2*center_c**3*half_r
        - center_c**2*half_q*half_r**2
        - 4*center_c**2*half_q
    )
    k3x8 = (
        center_a**4*center_c - 12*center_a**3*center_c
        - 12*center_a**2*center_c**2 + 32*center_a**2*center_c*half_q
        - 2*center_a**2*center_c*half_r**2
        - 16*center_a**2*center_c*half_r + 28*center_a**2*center_c
        + 16*center_a**2*half_q**2*half_r
        - 16*center_a*center_c**2*half_r - 8*center_a*center_c**2
        - 4*center_a*center_c*half_r**2 + 32*center_a*center_c*half_r
        - 16*center_a*center_c + 4*center_c**3
        - 4*center_c**2*half_r**2 - 16*center_c**2
        + center_c*half_r**4 + 8*center_c*half_r**2 + 16*center_c
    )
    equations = [ring(poly) for poly in (f1, f2, f3, k2, k3x8)]
    jacobian = Matrix(ring, [
        [equation.derivative(variable) for variable in ring.gens()]
        for equation in equations
    ])
    return ring, equations, jacobian


def contact6_data(field, values):
    a, b, contact_scale, contact_u, contact_v, half_q, half_r = values
    ring = PolynomialRing(field, "x")
    x = ring.gen()
    h6 = 1 + a*x + b*x**2 + x**3
    f = h6**2 - (x - 1)**6
    c1, c2, c3, c4, c5 = [f[i] for i in range(1, 6)]
    contact_b = c5*contact_scale**2 + 3*contact_u
    delta = (
        4*c4*contact_scale**2 + 12*(contact_u**2 + contact_v**2)
        - contact_b**2
    )
    f3 = (
        contact_b*delta + 16*contact_v**3 - 8*c3*contact_scale**2
        - 8*contact_u**3 - 48*contact_u*contact_v**2
    )
    f2 = (
        delta**2 + 64*contact_b*contact_v**3 - 64*c2*contact_scale**2
        - 192*(contact_u**2*contact_v**2 + contact_v**4)
    )
    f1 = delta*contact_v**3 - 4*c1*contact_scale**2 - 12*contact_u*contact_v**4
    A = b + 3
    C = a + b + 2
    k2 = (
        A**2*C**2*half_q - 2*A**2*C*half_q**2 + A**2*half_q**4
        - 2*A*C**3 - 2*A*C**2*half_q**2 + 2*A*C**2*half_q
        + 4*A*C*half_q**2 + 2*C**3*half_q - 2*C**3*half_r
        - C**2*half_q*half_r**2 - 4*C**2*half_q
    )
    k3x8 = (
        A**4*C - 12*A**3*C - 12*A**2*C**2 + 32*A**2*C*half_q
        - 2*A**2*C*half_r**2 - 16*A**2*C*half_r + 28*A**2*C
        + 16*A**2*half_q**2*half_r - 16*A*C**2*half_r - 8*A*C**2
        - 4*A*C*half_r**2 + 32*A*C*half_r - 16*A*C + 4*C**3
        - 4*C**2*half_r**2 - 16*C**2 + C*half_r**4
        + 8*C*half_r**2 + 16*C
    )
    return (f1, f2, f3, k2, k3x8), (ring, x, f, A, C, x**2 + contact_u*x + contact_v**2)


def boundary_labels(field_data, values):
    ring, x, f, A, C, contact_q = field_data
    _a, _b, L, U, v, q, _r = values
    labels = []
    if A == 0:
        labels.append("A0")
    if C == 0:
        labels.append("C0")
    if L == 0:
        labels.append("L0")
    if v == 0:
        labels.append("v0")
    if v**3 == 1:
        labels.append("v3eq1")
    if U**2 == 4*v**2:
        labels.append("contact_q_double")
    if q == 0:
        labels.append("half_q0")
    if f.degree() != 5:
        labels.append("degree_drop")
    elif f.discriminant() == 0:
        labels.append("singular_curve")
    if contact_q.gcd(f).degree() > 0:
        labels.append("contact_q_meets_f")
    return tuple(sorted(labels))


def integer_boundary_values(values):
    a, b, L, U, v, q, _r = map(ZZ, values)
    ring = PolynomialRing(ZZ, "x")
    x = ring.gen()
    h6 = 1 + a*x + b*x**2 + x**3
    f = h6**2 - (x - 1)**6
    contact_q = x**2 + U*x + v**2
    return {
        "A0": b + 3,
        "C0": a + b + 2,
        "L0": L,
        "v0": v,
        "v3eq1": v**3 - 1,
        "contact_q_double": U**2 - 4*v**2,
        "half_q0": q,
        "degree_drop": b + 3,
        "singular_curve": f.discriminant(),
        "contact_q_meets_f": f.resultant(contact_q),
    }


def contact6_audit():
    print("\n=== contact-6 [6,12] rank-five boundary audit at p=5 ===")
    integer_ring, integer_equations, integer_jacobian = contact6_symbolic_data()
    field = GF(5)
    solutions = []
    for first_five in product(field, repeat=5):
        seed = tuple(first_five) + (field(0), field(0))
        equations, _ = contact6_data(field, seed)
        if any(equations[index] for index in range(3)):
            continue
        for q, r in product(field, repeat=2):
            values = tuple(first_five) + (q, r)
            equations, data = contact6_data(field, values)
            if equations[3] or equations[4]:
                continue
            point = [ZZ(value) for value in values]
            jacobian = matrix(field, [
                [entry(*point) for entry in row]
                for row in integer_jacobian.rows()
            ])
            if jacobian.rank() == 5:
                solutions.append((values, data, jacobian))

    print(f"rank_five_points={len(solutions)}")
    signature_points = Counter()
    signature_lifts = Counter()
    signature_full_exit = Counter()
    signature_smooth_exit = Counter()
    signature_exit_sets = defaultdict(Counter)
    signature_second_exit_sets = defaultdict(Counter)
    signature_second_full = Counter()
    point_records = []

    for index, (values, data, jacobian) in enumerate(solutions, start=1):
        labels = boundary_labels(data, values)
        signature = "+".join(labels)
        signature_points[signature] += 1
        point = [ZZ(value) for value in values]
        constants = []
        for equation in integer_equations:
            value = ZZ(equation(*point))
            assert value % 5 == 0
            constants.append(field(-(value // 5)))
        rhs = vector(field, constants)
        assert jacobian.augment(rhs.column()).rank() == 5
        particular = jacobian.solve_right(rhs)
        kernel_basis = jacobian.right_kernel().basis()
        assert len(kernel_basis) == 2
        lift_records = []
        smooth_sample = None
        second_order_full_count = 0
        second_order_sample = None
        for c0, c1 in product(field, repeat=2):
            correction = particular + c0*kernel_basis[0] + c1*kernel_basis[1]
            lift = tuple(point[i] + 5*ZZ(correction[i]) for i in range(7))
            assert all(ZZ(equation(*lift)) % 25 == 0 for equation in integer_equations)
            boundary_values = integer_boundary_values(lift)
            exited = tuple(
                label for label in labels if ZZ(boundary_values[label]) % 25 != 0
            )
            full_exit = len(exited) == len(labels)
            smooth_exit = (
                ZZ(boundary_values["A0"]) % 25 != 0
                and ZZ(boundary_values["singular_curve"]) % 25 != 0
            )
            signature_lifts[signature] += 1
            signature_full_exit[signature] += full_exit
            signature_smooth_exit[signature] += smooth_exit
            signature_exit_sets[signature][exited] += 1
            lift_records.append((exited, full_exit, smooth_exit))
            if smooth_exit and smooth_sample is None:
                smooth_sample = tuple(value % 25 for value in lift)
            if smooth_exit:
                second_constants = []
                for equation in integer_equations:
                    equation_value = ZZ(equation(*lift))
                    assert equation_value % 25 == 0
                    second_constants.append(field(-(equation_value // 25)))
                second_rhs = vector(field, second_constants)
                assert jacobian.augment(second_rhs.column()).rank() == 5
                second_particular = jacobian.solve_right(second_rhs)
                for d0, d1 in product(field, repeat=2):
                    second_correction = (
                        second_particular + d0*kernel_basis[0] + d1*kernel_basis[1]
                    )
                    lift125 = tuple(
                        lift[i] + 25*ZZ(second_correction[i]) for i in range(7)
                    )
                    assert all(
                        ZZ(equation(*lift125)) % 125 == 0
                        for equation in integer_equations
                    )
                    second_values = integer_boundary_values(lift125)
                    second_exited = tuple(
                        label for label in labels
                        if ZZ(second_values[label]) % 125 != 0
                    )
                    second_full = len(second_exited) == len(labels)
                    second_order_full_count += second_full
                    signature_second_full[signature] += second_full
                    signature_second_exit_sets[signature][second_exited] += 1
                    if second_full and second_order_sample is None:
                        second_order_sample = tuple(value % 125 for value in lift125)
        full_count = sum(record[1] for record in lift_records)
        smooth_count = sum(record[2] for record in lift_records)
        point_records.append((values, signature, full_count, smooth_count))
        print(
            f" rank5_point_{index} values={tuple(ZZ(v) for v in values)} "
            f"signature={signature} lifts=25 full_transverse={full_count} "
            f"smooth_curve_first_order={smooth_count} smooth_sample_mod25={smooth_sample} "
            f"second_order_full={second_order_full_count} "
            f"second_order_sample_mod125={second_order_sample}"
        )

    print("rank5_signature_summary")
    for signature in sorted(signature_points):
        print(
            f" signature={signature} points={signature_points[signature]} "
            f"lifts={signature_lifts[signature]} "
            f"full_transverse={signature_full_exit[signature]} "
            f"smooth_curve_first_order={signature_smooth_exit[signature]}"
        )
        print(f"  exit_sets={dict(signature_exit_sets[signature])}")
        if signature_second_exit_sets[signature]:
            print(
                f"  second_order_full={signature_second_full[signature]} "
                f"second_exit_sets={dict(signature_second_exit_sets[signature])}"
            )

    minimal = sorted(
        signature_points,
        key=lambda signature: (
            -signature_full_exit[signature],
            len(signature.split("+")),
            signature,
        ),
    )
    print("contact6_decision")
    print(f" signatures_with_full_transverse_lift={[s for s in minimal if signature_full_exit[s]]}")
    print(
        " first_order_trapped_signatures="
        f"{[s for s in minimal if not signature_full_exit[s]]}"
    )
    print(" selected_minimal_signatures", minimal[:5])


def main():
    contact9_audit()
    contact6_audit()


if __name__ == "__main__":
    main()
