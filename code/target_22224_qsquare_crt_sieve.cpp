// Boundary-aware rational sieve for the q-square [2,2,2,24] cover.
//
// Reads the exhaustive finite masks written by
// target_22224_qsquare_crt_sieve.m and scans positive reduced A,B of
// bounded numerator/denominator height.  The S_3 and square-root sign
// symmetries are applied before consulting a target mask.
//
// A prime is used rigorously as follows.  If A,B,C,rho,sigma and the
// branch discriminant are affine/open, the finite marked class must occur
// in the target mask.  Otherwise the candidate is retained as a boundary
// reduction.  Thus the empty p=31 affine mask becomes a necessary
// p=31-boundary condition rather than a false global obstruction.
//
// Build/run from the repository root:
//   c++ -O3 -std=c++17 -o /tmp/target_22224_qsquare_crt_sieve \
//     code/target_22224_qsquare_crt_sieve.cpp
//   /tmp/target_22224_qsquare_crt_sieve 187 \
//     results/target_22224_qsquare_crt_sieve \
//     results/target_22224_qsquare_crt_sieve_h187_candidates.tsv

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using i128 = __int128_t;
using u128 = __uint128_t;

struct Fraction {
    int n;
    int d;
    std::array<int,4> num_mod{};
    std::array<int,4> den_mod{};
};

struct PrimeMask {
    int p=0;
    std::vector<unsigned char> target;
};

static int mod_norm(long long x,int p) {
    int r=static_cast<int>(x%p);
    return r<0 ? r+p : r;
}

static int mod_pow(int a,int e,int p) {
    long long r=1,b=mod_norm(a,p);
    while(e) {
        if(e&1) r=r*b%p;
        b=b*b%p;
        e>>=1;
    }
    return static_cast<int>(r);
}

static int mod_inv(int a,int p) {
    return mod_pow(a,p-2,p);
}

static bool mod_square_or_zero(int a,int p) {
    a=mod_norm(a,p);
    return a==0 || mod_pow(a,(p-1)/2,p)==1;
}

static std::string i128_string(i128 x) {
    if(x==0) return "0";
    bool neg=x<0;
    u128 y=neg ? static_cast<u128>(-x) : static_cast<u128>(x);
    std::string s;
    while(y) {
        s.push_back(char('0'+y%10));
        y/=10;
    }
    if(neg) s.push_back('-');
    std::reverse(s.begin(),s.end());
    return s;
}

static unsigned long long sqrt_floor_i128(i128 x) {
    if(x<=0) return 0;
    long double xd=static_cast<long double>(x);
    unsigned long long r=static_cast<unsigned long long>(std::sqrt(xd));
    while(static_cast<u128>(r)*r>static_cast<u128>(x)) --r;
    while(static_cast<u128>(r+1)*(r+1)<=static_cast<u128>(x)) ++r;
    return r;
}

static bool square_i128(i128 x,unsigned long long &r) {
    if(x<0) return false;
    r=sqrt_floor_i128(x);
    return static_cast<u128>(r)*r==static_cast<u128>(x);
}

static std::vector<std::string> split_tab(const std::string &s) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string v;
    while(std::getline(ss,v,'\t')) out.push_back(v);
    return out;
}

static PrimeMask load_mask(const std::string &stem,int p) {
    PrimeMask pm;
    pm.p=p;
    pm.target.assign(p*p,0);
    std::string name=stem+"_p"+std::to_string(p)+"_masks.tsv";
    std::ifstream in(name);
    if(!in) {
        std::cerr << "cannot open mask " << name << "\n";
        std::exit(2);
    }
    std::string line;
    std::getline(in,line);
    while(std::getline(in,line)) {
        auto f=split_tab(line);
        if(f.size()<13) continue;
        int A=std::stoi(f[1]);
        int B=std::stoi(f[2]);
        int target_presentations=std::stoi(f[10]);
        if(target_presentations>0) pm.target[A*p+B]=1;
    }
    return pm;
}

static bool raw_target(const PrimeMask &pm,int A,int B) {
    int p=pm.p;
    A=mod_norm(A,p); B=mod_norm(B,p);
    return A!=0 && B!=0 && pm.target[A*p+B];
}

// Union over all S_3 permutations of (A,B,C) and all sign patterns with
// product +1.  The missing global sign is the rho-sheet involution already
// present in the finite enumeration.
static bool symmetry_target(const PrimeMask &pm,int A,int B,int C) {
    int p=pm.p;
    std::array<int,3> x{A,B,C};
    const int perms[6][3]={{0,1,2},{0,2,1},{1,0,2},
                           {1,2,0},{2,0,1},{2,1,0}};
    const int signs[4][3]={{1,1,1},{1,-1,-1},
                           {-1,1,-1},{-1,-1,1}};
    for(const auto &sg:signs) {
        std::array<int,3> y{
            mod_norm(sg[0]*x[0],p),mod_norm(sg[1]*x[1],p),
            mod_norm(sg[2]*x[2],p)};
        for(const auto &pe:perms) {
            if(raw_target(pm,y[pe[0]],y[pe[1]])) return true;
        }
    }
    return false;
}

enum class PrimeStatus { Boundary, Target, NonSquare, SmoothNoTarget };

static const char *status_name(PrimeStatus s) {
    switch(s) {
        case PrimeStatus::Boundary: return "boundary";
        case PrimeStatus::Target: return "target";
        case PrimeStatus::NonSquare: return "nonsquare";
        default: return "smooth_no_target";
    }
}

static PrimeStatus prime_status(const PrimeMask &pm,const Fraction &A,
                                const Fraction &B,int prime_index) {
    int p=pm.p;
    int an=A.num_mod[prime_index], ad=A.den_mod[prime_index];
    int bn=B.num_mod[prime_index], bd=B.den_mod[prime_index];
    // A zero/pole (and hence some zero/pole among A,B,C) is a toric boundary.
    if(an==0 || ad==0 || bn==0 || bd==0) return PrimeStatus::Boundary;
    int ar=static_cast<int>((long long)an*mod_inv(ad,p)%p);
    int br=static_cast<int>((long long)bn*mod_inv(bd,p)%p);
    int cr=mod_inv(static_cast<int>((long long)ar*br%p),p);
    int a2=static_cast<int>((long long)ar*ar%p);
    int b2=static_cast<int>((long long)br*br%p);
    int c2=static_cast<int>((long long)cr*cr%p);
    int R=mod_norm(a2+b2+c2-3,p);
    int S=mod_norm((long long)mod_inv(a2,p)+mod_inv(b2,p)+mod_inv(c2,p)-3,p);
    if(!mod_square_or_zero(R,p) || !mod_square_or_zero(S,p))
        return PrimeStatus::NonSquare;
    // rho=0 or sigma=0 is a sheet boundary.  It is retained, not treated
    // as an affine target.
    if(R==0 || S==0) return PrimeStatus::Boundary;
    int d2=static_cast<int>((long long)S*mod_inv(R,p)%p); // (d/s)^2
    std::array<int,4> roots{a2,b2,c2,d2};
    std::sort(roots.begin(),roots.end());
    if(std::adjacent_find(roots.begin(),roots.end())!=roots.end())
        return PrimeStatus::Boundary;
    if(symmetry_target(pm,ar,br,cr)) return PrimeStatus::Target;
    return PrimeStatus::SmoothNoTarget;
}

static bool fraction_leq(const Fraction &a,const Fraction &b) {
    return static_cast<long long>(a.n)*b.d<=static_cast<long long>(b.n)*a.d;
}

// B <= C=1/(AB), avoiding construction of C.
static bool B_leq_C(const Fraction &A,const Fraction &B) {
    i128 lhs=static_cast<i128>(A.n)*B.n*B.n;
    i128 rhs=static_cast<i128>(A.d)*B.d*B.d;
    return lhs<=rhs;
}

struct ExactSquareData {
    bool ok=false;
    unsigned long long rho_num=0,rho_den=1;
    unsigned long long sigma_num=0,sigma_den=1;
};

static ExactSquareData exact_square_data(const Fraction &A,const Fraction &B) {
    i128 an=A.n,ad=A.d,bn=B.n,bd=B.d;
    i128 an2=an*an,ad2=ad*ad,bn2=bn*bn,bd2=bd*bd;
    i128 D=an2*ad2*bn2*bd2;
    i128 Rnum=bd2*an2*an2*bn2
             +ad2*an2*bn2*bn2
             +ad2*ad2*bd2*bd2-3*D;
    i128 Snum=ad2*ad2*bd2*bn2
             +ad2*bd2*bd2*an2
             +an2*an2*bn2*bn2-3*D;
    unsigned long long rr=0,ss=0;
    ExactSquareData z;
    if(!square_i128(Rnum,rr) || !square_i128(Snum,ss) || rr==0 || ss==0)
        return z;
    unsigned long long denroot=static_cast<unsigned long long>(A.n)*A.d*B.n*B.d;
    unsigned long long gr=std::gcd(rr,denroot);
    unsigned long long gs=std::gcd(ss,denroot);
    z.ok=true;
    z.rho_num=rr/gr; z.rho_den=denroot/gr;
    z.sigma_num=ss/gs; z.sigma_den=denroot/gs;
    return z;
}

int main(int argc,char **argv) {
    int H=argc>1 ? std::stoi(argv[1]) : 187;
    std::string stem=argc>2 ? argv[2] :
      "results/target_22224_qsquare_crt_sieve";
    std::string output=argc>3 ? argv[3] :
      "results/target_22224_qsquare_crt_sieve_h187_candidates.tsv";
    std::string summary_name=argc>4 ? argv[4] :
      "results/target_22224_qsquare_crt_sieve_h187.log";
    const std::array<int,4> primes{29,31,37,41};
    std::array<PrimeMask,4> masks{
        load_mask(stem,29),load_mask(stem,31),
        load_mask(stem,37),load_mask(stem,41)};

    std::ofstream summary(summary_name);
    if(!summary) {
        std::cerr << "cannot open summary " << summary_name << "\n";
        return 2;
    }

    std::vector<Fraction> vals;
    for(int d=1;d<=H;++d) for(int n=1;n<=H;++n) {
        if(std::gcd(n,d)!=1) continue;
        Fraction f{n,d};
        for(int k=0;k<4;++k) {
            f.num_mod[k]=n%primes[k];
            f.den_mod[k]=d%primes[k];
        }
        vals.push_back(f);
    }
    std::sort(vals.begin(),vals.end(),[](const Fraction &a,const Fraction &b){
        i128 lhs=static_cast<i128>(a.n)*b.d;
        i128 rhs=static_cast<i128>(b.n)*a.d;
        if(lhs!=rhs) return lhs<rhs;
        return a.d<b.d;
    });

    std::ofstream out(output);
    if(!out) {
        std::cerr << "cannot open output " << output << "\n";
        return 2;
    }
    out << "A_num\tA_den\tB_num\tB_den\tC_num\tC_den"
           "\trho_num\trho_den\tsigma_num\tsigma_den"
           "\tp29\tp31\tp37\tp41\n";

    unsigned long long ordered_pairs=0,trivial_or_exact_collision=0;
    unsigned long long modular_survivors=0,exact_double_square=0;
    std::array<unsigned long long,4> nonsquare_reject{};
    std::array<unsigned long long,4> smooth_reject{};
    std::array<unsigned long long,4> boundary_pass{};
    std::array<unsigned long long,4> target_pass{};
    unsigned long long written=0;

    for(std::size_t i=0;i<vals.size();++i) {
      const Fraction &A=vals[i];
      for(std::size_t j=i;j<vals.size();++j) {
        const Fraction &B=vals[j];
        if(!B_leq_C(A,B)) break;
        ++ordered_pairs;
        // Exact equality among A,B,C gives a repeated branch root.  B=1
        // is precisely the elementary {r,1,1/r} double-square component.
        bool AeqB=static_cast<long long>(A.n)*B.d==
                  static_cast<long long>(B.n)*A.d;
        bool BeqC=static_cast<i128>(A.n)*B.n*B.n==
                  static_cast<i128>(A.d)*B.d*B.d;
        bool Bis1=B.n==B.d;
        if(AeqB || BeqC || Bis1) {
            ++trivial_or_exact_collision;
            continue;
        }

        std::array<PrimeStatus,4> statuses;
        bool pass=true;
        for(int k=0;k<4;++k) {
            statuses[k]=prime_status(masks[k],A,B,k);
            if(statuses[k]==PrimeStatus::NonSquare) {
                ++nonsquare_reject[k]; pass=false; break;
            }
            if(statuses[k]==PrimeStatus::SmoothNoTarget) {
                ++smooth_reject[k]; pass=false; break;
            }
            if(statuses[k]==PrimeStatus::Boundary) ++boundary_pass[k];
            else ++target_pass[k];
        }
        if(!pass) continue;
        ++modular_survivors;

        // The common denominator for R and S is already a square, so this
        // exact test reduces to two integer-square tests.
        ExactSquareData ex=exact_square_data(A,B);
        if(!ex.ok) continue;
        ++exact_double_square;
        long long Cnum=static_cast<long long>(A.d)*B.d;
        long long Cden=static_cast<long long>(A.n)*B.n;
        long long gc=std::gcd(Cnum,Cden);
        Cnum/=gc; Cden/=gc;
        auto status_char=[](PrimeStatus s)->char {
            return s==PrimeStatus::Boundary ? 'B' :
                   (s==PrimeStatus::Target ? 'T' : 'X');
        };
        out << A.n << '\t' << A.d << '\t' << B.n << '\t' << B.d << '\t'
            << Cnum << '\t' << Cden << '\t'
            << ex.rho_num << '\t' << ex.rho_den << '\t'
            << ex.sigma_num << '\t' << ex.sigma_den;
        for(int k=0;k<4;++k) out << '\t' << status_char(statuses[k]);
        out << '\n';
        ++written;
    }}
    out.close();

    // Positive control: the record is in the height-187 box in canonical
    // order (13/187,17/7,77/13), but must be rejected by the p=31 open mask.
    Fraction RA{13,187},RB{17,7};
    for(int k=0;k<4;++k) {
        RA.num_mod[k]=RA.n%primes[k]; RA.den_mod[k]=RA.d%primes[k];
        RB.num_mod[k]=RB.n%primes[k]; RB.den_mod[k]=RB.d%primes[k];
    }
    ExactSquareData rex=exact_square_data(RA,RB);
    std::ostringstream control;
    control << "RECORD_HEIGHT187_CONTROL exact_double_square " << (rex.ok?1:0);
    for(int k=0;k<4;++k)
        control << " p" << primes[k] << "_" << status_name(prime_status(masks[k],RA,RB,k));
    std::cout << control.str() << "\n";
    summary << control.str() << "\n";

    std::ostringstream headline;
    headline << "QSQUARE_BOUNDARY_SIEVE"
             << " height " << H
             << " positive_values " << vals.size()
             << " ordered_pairs " << ordered_pairs
             << " exact_collision_or_trivial " << trivial_or_exact_collision
             << " modular_survivors " << modular_survivors
             << " exact_double_square_survivors " << exact_double_square
             << " written " << written;
    std::cout << headline.str() << "\n";
    summary << headline.str() << "\n";
    for(int k=0;k<4;++k) {
        std::ostringstream line;
        line << "PRIME_MASK " << primes[k]
             << " nonsquare_reject " << nonsquare_reject[k]
             << " smooth_no_target_reject " << smooth_reject[k]
             << " boundary_pass " << boundary_pass[k]
             << " target_pass " << target_pass[k];
        std::cout << line.str() << "\n";
        summary << line.str() << "\n";
    }
    std::cout << "OUTPUT " << output << "\n";
    summary << "OUTPUT " << output << "\n";
    summary.close();
    return 0;
}
