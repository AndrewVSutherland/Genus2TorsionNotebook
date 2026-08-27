#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

using i128 = __int128_t;

static i128 abs128(i128 x) { return x < 0 ? -x : x; }

static i128 gcd128(i128 a, i128 b) {
    a = abs128(a);
    b = abs128(b);
    while (b != 0) {
        i128 r = a % b;
        a = b;
        b = r;
    }
    return a;
}

static std::string to_string128(i128 x) {
    if (x == 0) return "0";
    bool neg = x < 0;
    if (neg) x = -x;
    std::string s;
    while (x > 0) {
        int digit = static_cast<int>(x % 10);
        s.push_back(static_cast<char>('0' + digit));
        x /= 10;
    }
    if (neg) s.push_back('-');
    std::reverse(s.begin(), s.end());
    return s;
}

static bool is_square128(i128 n, i128& root) {
    if (n < 0) return false;
    i128 lo = 0;
    i128 hi = 1;
    while (hi <= n / hi) hi <<= 1;
    while (lo + 1 < hi) {
        i128 mid = (lo + hi) >> 1;
        if (mid <= n / mid) lo = mid;
        else hi = mid;
    }
    root = lo;
    return lo * lo == n;
}

static long long mod128(i128 x, long long m) {
    i128 r = x % static_cast<i128>(m);
    if (r < 0) r += m;
    return static_cast<long long>(r);
}

static std::vector<bool> square_residue_table(long long modulus) {
    std::vector<bool> ok(static_cast<size_t>(modulus), false);
    for (long long x = 0; x < modulus; ++x) {
        ok[static_cast<size_t>((static_cast<__int128_t>(x) * x) % modulus)] = true;
    }
    return ok;
}

static bool make_curve_tuple(i128 a, i128 b, i128 c, i128 d_num, i128 d_den,
                             long long max_abs, std::vector<i128>& curve) {
    if (d_den == 0) return false;
    if (d_den < 0) {
        d_den = -d_den;
        d_num = -d_num;
    }

    i128 gden = gcd128(d_num, d_den);
    d_num /= gden;
    d_den /= gden;

    std::vector<i128> signed_tuple = {a * d_den, b * d_den, c * d_den, d_num};
    i128 g = 0;
    for (i128 x : signed_tuple) g = gcd128(g, x);
    if (g == 0) return false;

    curve.clear();
    for (i128 x : signed_tuple) {
        i128 y = abs128(x / g);
        if (y == 0 || y > max_abs) return false;
        curve.push_back(y);
    }
    std::sort(curve.begin(), curve.end());
    for (int i = 1; i < 4; ++i) {
        if (curve[i] == curve[i - 1]) return false;
    }
    return true;
}

static std::string tuple_key(const std::vector<i128>& v) {
    std::ostringstream os;
    os << "[";
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) os << ",";
        os << to_string128(v[i]);
    }
    os << "]";
    return os.str();
}

static void maybe_record(const std::vector<i128>& curve,
                         std::unordered_set<std::string>& seen,
                         std::vector<std::string>& rows) {
    std::string key = tuple_key(curve);
    if (seen.insert(key).second) {
        rows.push_back(key + "\n");
    }
}

static bool is_square_mod_filter(
        i128 n,
        const std::vector<std::pair<long long, std::vector<bool>>>& square_residues) {
    if (n < 0) return false;
    for (const auto& item : square_residues) {
        long long p = item.first;
        const std::vector<bool>& table = item.second;
        if (!table[static_cast<size_t>(mod128(n, p))]) return false;
    }
    return true;
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "usage: " << argv[0]
                  << " H OUTPUT_FILE MAX_ABS [A_START A_END]\n";
        return 2;
    }

    const long long p = 13;
    long long H = std::atoll(argv[1]);
    std::string output = argv[2];
    long long max_abs = std::atoll(argv[3]);
    long long A_start = -H;
    long long A_end = H;
    if (argc >= 6) {
        A_start = std::max(-H, std::atoll(argv[4]));
        A_end = std::min(H, std::atoll(argv[5]));
    }

    std::ofstream out(output);
    if (!out) {
        std::cerr << "could not open output file: " << output << "\n";
        return 1;
    }

    std::vector<std::pair<long long, std::vector<bool>>> square_residues;
    for (long long ell : {3LL, 5LL, 7LL, 11LL, 17LL, 19LL, 23LL, 29LL, 31LL, 37LL}) {
        square_residues.push_back({ell, square_residue_table(ell)});
    }

    std::unordered_set<std::string> seen;
    std::vector<std::string> rows;
    std::vector<i128> curve;

    i128 checked = 0;
    i128 mod_square = 0;
    i128 exact_square = 0;
    i128 tuple_rows = 0;

    for (long long A = A_start; A <= A_end; ++A) {
        if (A == 0) continue;
        for (long long B = -H; B <= H; ++B) {
            if (B == 0) continue;
            if (A + B == 0) continue;

            long long kmin = (-H + A + B + p - 1) / p;
            long long kmax = (H + A + B) / p;
            if (kmin < -H) kmin = -H;
            if (kmax > H) kmax = H;

            for (long long K = kmin; K <= kmax; ++K) {
                long long C = -A - B + p * K;
                if (C == 0 || C < -H || C > H) continue;
                checked += 1;

                i128 disc16 = -static_cast<i128>(A) * B
                    * (p * static_cast<i128>(K) - A)
                    * (static_cast<i128>(A) + B)
                    * (p * static_cast<i128>(K) - B)
                    * (p * static_cast<i128>(K) - A - B);

                if (!is_square_mod_filter(disc16, square_residues)) continue;
                mod_square += 1;

                i128 r = 0;
                if (!is_square128(disc16, r)) continue;
                exact_square += 1;

                i128 a = p * static_cast<i128>(A);
                i128 b = p * static_cast<i128>(B);
                i128 c = p * static_cast<i128>(C);
                i128 S = static_cast<i128>(A) + B + C; // = p*K
                i128 T = static_cast<i128>(A) * B + static_cast<i128>(A) * C
                    + static_cast<i128>(B) * C;
                i128 ABC = static_cast<i128>(A) * B * C;

                if (S == 0) {
                    // Equation is p*T^2 = 4*ABC*d.
                    i128 d_num = p * T * T;
                    i128 d_den = 4 * ABC;
                    if (make_curve_tuple(a, b, c, d_num, d_den, max_abs, curve)) {
                        maybe_record(curve, seen, rows);
                        tuple_rows += 1;
                    }
                    continue;
                }

                i128 den = p * S * S;
                for (int eps : {-1, 1}) {
                    i128 num = 2 * p * ABC - p * p * S * T + eps * 2 * p * r;
                    if (make_curve_tuple(a, b, c, num, den, max_abs, curve)) {
                        maybe_record(curve, seen, rows);
                        tuple_rows += 1;
                    }
                }
            }
        }

        if ((A - A_start) % 25 == 0) {
            std::cerr << "A=" << A
                      << " checked=" << to_string128(checked)
                      << " mod_square=" << to_string128(mod_square)
                      << " exact_square=" << to_string128(exact_square)
                      << " unique=" << seen.size() << "\n";
        }
    }

    std::sort(rows.begin(), rows.end());
    for (const std::string& row : rows) out << row;

    std::cerr << "SUMMARY H=" << H
              << " max_abs=" << max_abs
              << " A_range=" << A_start << ".." << A_end
              << " checked=" << to_string128(checked)
              << " mod_square=" << to_string128(mod_square)
              << " exact_square=" << to_string128(exact_square)
              << " tuple_rows=" << to_string128(tuple_rows)
              << " unique=" << seen.size()
              << " output=" << output << "\n";
    return 0;
}
