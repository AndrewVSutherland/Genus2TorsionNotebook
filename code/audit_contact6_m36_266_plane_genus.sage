#!/usr/bin/env sage
"""Independent Sage/Singular genus audit for the extracted plane factors."""

import argparse

from sage.all import AffineSpace, Curve, GF, PolynomialRing, QQ
from sage.misc.sage_eval import sage_eval


def load_factors(path, field):
    ambient = AffineSpace(field, 2, names=("U", "v"))
    ring = ambient.coordinate_ring()
    source_ring = PolynomialRing(QQ, 2, names=("U", "v"))
    source_U, source_v = source_ring.gens()
    expressions = [
        line.strip().split("factor_polynomial ", 1)[1]
        for line in open(path, encoding="ascii")
        if line.startswith("factor_polynomial ")
    ]
    factors = []
    for expression in expressions:
        source = source_ring(
            sage_eval(expression, locals={"UU": source_U, "vv": source_v})
        )
        factors.append(ring(source) if field is QQ else ring(source.change_ring(field)))
    return ambient, factors


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        default="data/audit_contact6_m36_266_plane_extract_Q.txt",
    )
    parser.add_argument("--prime", type=int, default=0)
    parser.add_argument("--degrees", default="11,21")
    args = parser.parse_args()

    field = QQ if args.prime == 0 else GF(args.prime)
    ambient, factors = load_factors(args.input, field)
    requested = {int(value) for value in args.degrees.split(",")}

    print(f"GENUS_AUDIT field={field}", flush=True)
    for polynomial in factors:
        degree = polynomial.total_degree()
        if degree not in requested:
            continue
        factorization = polynomial.factor()
        print(
            f"factor degree={degree} deg_U={polynomial.degree(ambient.coordinate_ring().gen(0))} "
            f"terms={len(polynomial.monomials())} factorization={factorization}",
            flush=True,
        )
        curve = Curve(polynomial, ambient).projective_closure()
        arithmetic_genus = curve.arithmetic_genus()
        print(
            f"projective degree={degree} arithmetic_genus={arithmetic_genus} "
            f"singular={curve.is_singular()}",
            flush=True,
        )
        geometric_genus = curve.geometric_genus()
        print(
            f"geometric_genus={geometric_genus} delta_total={arithmetic_genus-geometric_genus}",
            flush=True,
        )


if __name__ == "__main__":
    main()
