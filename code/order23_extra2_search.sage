# Condition-driven search for rational 2-torsion on the Kuru-Sadek [23] family.
# Usage: sage code/order23_extra2_search.sage 2000

from sage.all import *
import sys

R.<x> = PolynomialRing(QQ)

GOOD_Q = {
    7: {4},
    11: {3, 5},
    13: {4, 9, 10},
    17: {2, 8, 13, 15},
    19: {5, 6, 7, 11, 17},
    23: {4, 8, 12},
    29: {4, 5, 6, 7, 9, 13, 16, 20, 23, 24, 25},
    31: {2, 5, 7, 8, 9, 10, 14, 16, 18, 20, 25},
    37: {3, 4, 7, 10, 12, 16, 21, 26, 27, 28, 30, 33, 34},
    41: {4, 8, 9, 16, 20, 21, 25, 33, 36, 37, 39},
    43: {6, 9, 10, 11, 15, 17, 21, 25, 35, 36, 40, 41},
    47: {2, 3, 4, 7, 8, 12, 16, 21, 25, 27, 32, 36, 37, 42},
    53: {6, 7, 11, 13, 15, 17, 24, 25, 36, 37, 40, 42, 43, 49},
    59: {3, 4, 7, 9, 12, 15, 16, 17, 22, 25, 27, 28, 29, 45, 49, 51, 53, 57},
    61: {3, 9, 13, 15, 16, 19, 22, 34, 42, 45, 48, 49, 52, 57},
    67: {9, 10, 14, 16, 17, 19, 21, 22, 23, 25, 26, 29, 33, 35, 36, 39, 40, 54, 56, 59, 60, 65},
    71: {2, 4, 5, 6, 8, 9, 10, 12, 15, 16, 18, 24, 29, 30, 32, 43, 45, 48, 50, 54, 57, 60},
    73: {3, 6, 9, 16, 18, 19, 25, 27, 32, 35, 36, 38, 41, 49, 50, 61, 64, 65, 67, 70},
    79: {2, 5, 8, 9, 11, 16, 18, 19, 20, 21, 22, 25, 26, 36, 38, 42, 44, 45, 49, 50, 51, 52, 55, 64, 67, 72, 76},
    83: {3, 4, 7, 9, 10, 11, 12, 16, 17, 21, 23, 25, 28, 29, 30, 33, 36, 37, 38, 40, 48, 51, 61, 63, 64, 68, 69, 70, 75, 77, 81},
    89: {2, 4, 5, 8, 9, 10, 16, 17, 21, 25, 34, 36, 40, 42, 45, 49, 50, 53, 55, 57, 64, 67, 71, 72, 73, 79, 80, 84, 85},
    97: {3, 4, 8, 9, 12, 22, 24, 25, 27, 31, 32, 33, 35, 36, 44, 47, 48, 50, 53, 62, 64, 65, 66, 70, 73, 75, 79, 81, 85, 89, 95},
}


def boundary_q_int(q, p):
    q %= p
    if q in (0, 1, p - 1):
        return True
    if p != 5 and (5*pow(q, 4, p) + 10*pow(q, 2, p) + 1) % p == 0:
        return True
    return False


ALLOWED_T = {}
for p, good in GOOD_Q.items():
    allowed = set()
    for a in range(p):
        q = (a*a) % p
        if boundary_q_int(q, p) or q in good:
            allowed.add(a)
    ALLOWED_T[p] = allowed


def passes_t_filter_numden(num, den):
    for p, allowed in ALLOWED_T.items():
        dm = den % p
        if dm == 0:
            continue
        a = (num % p) * inverse_mod(dm, p) % p
        if a not in allowed:
            return False
    return True


def kuru_sadek23(t):
    t = QQ(t)
    beta = (t**2 + 1)**2/(4*t**2)
    sbeta = (t**2 + 1)/(2*t)
    s = (t**2 - 1)/(2*t)
    alpha = beta - s**5/(beta*sbeta)
    lam = (alpha - 1)**4/((alpha - beta)**2*alpha)
    expr = (x**3*(x - alpha)**2
            - (x - 1)*((x - 1)**4 - lam*(x - beta)**2*x))/(2*(x - alpha)*(x - beta))
    num = R(expr.numerator())
    den = R(expr.denominator())
    A, rem = num.quo_rem(den)
    if rem != 0:
        raise ArithmeticError('non-polynomial A')
    return R(A**2 - lam*x**4*(x - 1))


B = Integer(sys.argv[1]) if len(sys.argv) > 1 else Integer(500)
print('Kuru-Sadek [23] extra-2 search')
print('height', B)
print('allowed_t_sizes', {p: len(v) for p, v in ALLOWED_T.items()})

checked = survivors = factored = singular = hits = 0
for den in range(1, B + 1):
    for num in range(-B, B + 1):
        if gcd(num, den) != 1:
            continue
        if num == 0 or num == den or num == -den:
            continue
        checked += 1
        if not passes_t_filter_numden(num, den):
            continue
        survivors += 1
        t = QQ(num)/QQ(den)
        try:
            f = kuru_sadek23(t)
        except Exception as exc:
            print('SKIP_ERROR', t, repr(exc))
            continue
        if f.discriminant() == 0:
            singular += 1
            continue
        factored += 1
        roots = f.roots(QQ, multiplicities=True)
        if roots:
            hits += 1
            print('HIT t=%s roots=%s' % (t, roots))
            print('f=%s' % f)

print('checked', checked)
print('survivors', survivors)
print('factored', factored)
print('singular', singular)
print('hits', hits)
print('DONE')
