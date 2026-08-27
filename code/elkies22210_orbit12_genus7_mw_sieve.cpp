// Fast finite-field coefficient-box sieve for the genus-7 cover
//
//   v^2 = G(q),  z^2 = D(q),     s = 59/49.
//
// The input lattice is the certified full Mordell--Weil basis used by
// elkies22210_orbit12_genus7_mw_sieve.m.  This program performs only safe
// necessary local tests.  Feed its short survivor list to the Magma script
// for exact reconstruction and cubic/radicand tests.
//
// Build and run, for example:
//
//   c++ -O3 -std=c++17 code/elkies22210_orbit12_genus7_mw_sieve.cpp \
//       -o /tmp/elkies22210_g7_sieve
//   /tmp/elkies22210_g7_sieve 2000

#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

struct Point {
    int x = 0;
    int y = 0;
    bool infinity = true;
};

static int mod(long long a, int p) {
    a %= p;
    if (a < 0) a += p;
    return static_cast<int>(a);
}

static int power_mod(int a, int n, int p) {
    long long r = 1;
    long long b = mod(a, p);
    while (n > 0) {
        if (n & 1) r = r * b % p;
        b = b * b % p;
        n >>= 1;
    }
    return static_cast<int>(r);
}

static int inverse_mod(int a, int p) {
    if (mod(a, p) == 0) throw std::runtime_error("division by zero");
    return power_mod(a, p - 2, p);
}

static int rational_mod(long long numerator, long long denominator, int p) {
    return mod(numerator * inverse_mod(mod(denominator, p), p), p);
}

struct EllipticData {
    int p;
    int a2;
    int a4;
    int a6;
    int s;
    bool d_at_zero_square;
    std::vector<unsigned char> d_square;
    std::vector<Point> p1_multiples;
    std::vector<Point> p2_multiples;
    std::array<Point, 4> torsion_multiples;
};

static Point negate(const Point &a, int p) {
    if (a.infinity) return a;
    return Point{a.x, mod(-a.y, p), false};
}

static Point add(const Point &a, const Point &b, const EllipticData &e) {
    if (a.infinity) return b;
    if (b.infinity) return a;
    const int p = e.p;
    if (a.x == b.x && mod(a.y + b.y, p) == 0) return Point{};

    int lambda;
    if (a.x == b.x && a.y == b.y) {
        if (a.y == 0) return Point{};
        const int numerator = mod(3LL * a.x * a.x + 2LL * e.a2 * a.x + e.a4, p);
        lambda = mod(1LL * numerator * inverse_mod(2 * a.y, p), p);
    } else {
        lambda = mod(1LL * (b.y - a.y) * inverse_mod(b.x - a.x, p), p);
    }
    const int x3 = mod(1LL * lambda * lambda - e.a2 - a.x - b.x, p);
    const int y3 = mod(-a.y + 1LL * lambda * (a.x - x3), p);
    return Point{x3, y3, false};
}

static int discriminant_value(int q, int s, int p) {
    const int oneq = mod(1 + q, p);
    const int qm = mod(q - s, p);
    const int oneq2 = mod(1LL * oneq * oneq, p);
    const int oneq4 = mod(1LL * oneq2 * oneq2, p);
    const int s2 = mod(1LL * s * s, p);
    const int s3 = mod(1LL * s2 * s, p);
    long long d = 0;
    d += 1LL * oneq2 * s2;
    d -= 4LL * s3;
    d += 4LL * oneq4 * qm;
    d -= 27LL * oneq2 * qm * qm;
    d -= 18LL * oneq2 * s * qm;
    return mod(d, p);
}

static bool is_square(int a, int p) {
    a = mod(a, p);
    return a == 0 || power_mod(a, (p - 1) / 2, p) == 1;
}

static std::vector<Point> multiples(const Point &p0, int bound,
                                    const EllipticData &e) {
    std::vector<Point> out(2 * bound + 1);
    out[bound] = Point{};
    for (int n = 1; n <= bound; ++n) {
        out[bound + n] = add(out[bound + n - 1], p0, e);
        out[bound - n] = negate(out[bound + n], e.p);
    }
    return out;
}

static EllipticData make_data(int p, int bound) {
    EllipticData e;
    e.p = p;
    e.s = rational_mod(59, 49, p);
    e.a2 = mod(1LL * e.s * e.s - 2LL * e.s - 2, p);
    e.a4 = mod(-4LL * e.s * e.s * (e.s + 1), p);
    e.a6 = mod(4LL * e.s * e.s * (e.s + 1) * (e.s + 1), p);

    e.d_square.resize(p);
    for (int q = 0; q < p; ++q)
        e.d_square[q] = static_cast<unsigned char>(
            is_square(discriminant_value(q, e.s, p), p));
    e.d_at_zero_square = e.d_square[0] != 0;

    const Point p1{rational_mod(864, 1225, p),
                   rational_mod(182088, 42875, p), false};
    const Point p2{rational_mod(-944, 441, p),
                   rational_mod(-367688, 64827, p), false};
    const Point t4{0, rational_mod(12744, 2401, p), false};
    e.p1_multiples = multiples(p1, bound, e);
    e.p2_multiples = multiples(p2, bound, e);
    e.torsion_multiples[0] = Point{};
    for (int k = 1; k < 4; ++k)
        e.torsion_multiples[k] = add(e.torsion_multiples[k - 1], t4, e);
    return e;
}

static bool passes(const Point &point, const EllipticData &e) {
    if (point.infinity) return e.d_at_zero_square;
    if (point.x == 0) return true; // q=infinity; leading coefficient D is 4.
    const int numerator = mod(2LL * e.s * (e.s + 1), e.p);
    const int q = mod(1LL * numerator * inverse_mod(point.x, e.p), e.p);
    return e.d_square[q] != 0;
}

int main(int argc, char **argv) {
    const int bound = argc > 1 ? std::atoi(argv[1]) : 2000;
    if (bound < 0) {
        std::cerr << "bound must be nonnegative\n";
        return 2;
    }
    const std::vector<int> primes = {
        11,13,17,19,23,29,31,37,41,43,47,53,61,67,71,73,79,83,
        89,97,101,103,107,109,113,127,131,137,139,149,151,157,
        163,167,173,179,181,191,193,197,199,211,223,227,229,233,
        239,241,251,257,263,269,271,277,281,283,293,307,311,313,
        317,331,337,347,349
    };

    std::vector<EllipticData> local;
    local.reserve(primes.size());
    for (int p : primes) local.push_back(make_data(p, bound));

    const std::uint64_t width = static_cast<std::uint64_t>(2 * bound + 1);
    const std::uint64_t total = 4 * width * width;
    std::uint64_t survivors = 0;
    std::vector<std::array<int, 3>> survivor_tuples;

    for (int n1 = -bound; n1 <= bound; ++n1) {
        for (int n2 = -bound; n2 <= bound; ++n2) {
            for (int k = 0; k < 4; ++k) {
                bool keep = true;
                for (const auto &e : local) {
                    Point r = add(e.p1_multiples[n1 + bound],
                                  e.p2_multiples[n2 + bound], e);
                    r = add(r, e.torsion_multiples[k], e);
                    if (!passes(r, e)) {
                        keep = false;
                        break;
                    }
                }
                if (keep) {
                    ++survivors;
                    survivor_tuples.push_back({n1, n2, k});
                }
            }
        }
    }

    std::cout << "ELKIES22210_ORBIT12_GENUS7_FAST_MW_SIEVE\n";
    std::cout << "coefficient_bound " << bound << "\n";
    std::cout << "mw_elements_covered " << total << "\n";
    std::cout << "good_primes " << primes.size() << "\n";
    std::cout << "modular_survivors " << survivors << "\n";
    for (const auto &c : survivor_tuples)
        std::cout << "SURVIVOR " << c[0] << ' ' << c[1] << ' ' << c[2] << '\n';
    std::cout << "DONE\n";
    return 0;
}
