// Direct two-parameter rational search on the explicit Orbit-B halving cover
// of the Clebsch--Klein family.  This is not a bounded enumeration of integral
// CK points.
//
// Normalize the first marked coordinate to r1=1.  Choose arbitrary rational
// x=r2 and y=r3.  The two CK equations determine u=r4,v=r5 by
//
//   u+v = S = -1-x-y,
//   uv  = P = (x+y)(1+x)(1+y)/(1+x+y).
//
// Hence the CK completion is rational exactly when S^2-4P is a square.  The
// recovered tuple is then tested literally against
//
//   G0=-(1-y^2)(1-u^2)(1-v^2),
//   G3=(y^2-x^2)(1-u^2)(1-v^2),
//   G4=(u^2-x^2)(1-y^2)(1-v^2),
//   G5=(v^2-x^2)(1-y^2)(1-u^2).
//
// The run at height B tests every reduced nonzero x,y with numerator and
// denominator at most B (excluding squares already equal to r1^2), while u,v
// are solved exactly with no height bound.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

using i64=std::int64_t;
using i128=__int128_t;
using u128=__uint128_t;

static u128 ua(i128 x){return x<0?(u128)(-x):(u128)x;}
static u128 gcd128(u128 a,u128 b){while(b){u128 r=a%b;a=b;b=r;}return a;}
static std::string s128(i128 x){if(!x)return"0";bool n=x<0;u128 u=ua(x);std::string s;while(u){s.push_back(char('0'+u%10));u/=10;}if(n)s.push_back('-');std::reverse(s.begin(),s.end());return s;}

struct Q{
    i128 n,d;
    Q(i128 a=0,i128 b=1):n(a),d(b){
        if(!d)std::abort(); if(d<0){n=-n;d=-d;}
        u128 g=gcd128(ua(n),(u128)d);n/=i128(g);d/=i128(g);
    }
};
static Q add(const Q&a,const Q&b){u128 g=gcd128(a.d,b.d);i128 ad=a.d/i128(g),bd=b.d/i128(g);return Q(a.n*bd+b.n*ad,ad*b.d);}
static Q neg(const Q&a){return Q(-a.n,a.d);} static Q sub(const Q&a,const Q&b){return add(a,neg(b));}
static Q mul(const Q&a,const Q&b){u128 g1=gcd128(ua(a.n),b.d),g2=gcd128(ua(b.n),a.d);return Q((a.n/i128(g1))*(b.n/i128(g2)),(a.d/i128(g2))*(b.d/i128(g1)));}
static Q divide(const Q&a,const Q&b){if(!b.n)std::abort();return mul(a,Q(b.d,b.n));}
static Q sqr(const Q&a){return mul(a,a);} static bool eq(const Q&a,const Q&b){return a.n==b.n&&a.d==b.d;}
static std::string qs(const Q&a){return s128(a.n)+"/"+s128(a.d);}

static bool sq64[64],sq63[63],sq65[65];
static void initres(){for(unsigned z=0;z<320;++z){sq64[z*z%64]=1;sq63[z*z%63]=1;sq65[z*z%65]=1;}}
static bool isqrt128(u128 n,u128&r){long double d=std::sqrt((long double)n);std::uint64_t z=(std::uint64_t)d;while((u128)z*z>n)--z;while((u128)(z+1)*(z+1)<=n)++z;r=z;return(u128)z*z==n;}
static bool squareq(const Q&a){
    if(a.n<=0)return false;u128 n=a.n,d=a.d;
    if(!sq64[n%64]||!sq64[d%64]||!sq63[n%63]||!sq63[d%63]||!sq65[n%65]||!sq65[d%65])return false;
    u128 rn,rd;return isqrt128(n,rn)&&isqrt128(d,rd);
}
static bool sqrtq(const Q&a,Q&r){
    if(a.n<=0)return false;u128 n=a.n,d=a.d;
    if(!sq64[n%64]||!sq64[d%64]||!sq63[n%63]||!sq63[d%63]||!sq65[n%65]||!sq65[d%65])return false;
    u128 rn,rd;if(!isqrt128(n,rn)||!isqrt128(d,rd))return false;r=Q(i128(rn),i128(rd));return true;
}

struct RV{i64 n,d;Q q;};

int main(int argc,char**argv){
    int B=argc>1?std::atoi(argv[1]):40;
    if(B<2||B>200){std::cerr<<"require 2 <= B <= 200\n";return 2;}
    initres();
    std::vector<RV> vv;
    for(i64 d=1;d<=B;++d)for(i64 n=-B;n<=B;++n){
        if(!n||std::abs(n)==d||std::gcd(std::abs(n),d)!=1)continue;
        vv.push_back({n,d,Q(n,d)});
    }
    std::uint64_t pairs=0,distinct=0,denok=0,dpos=0,dsq=0,
                  smoothck=0,gpositive=0,g0sq=0,g3sq=0,g4sq=0,g5sq=0,hits=0;
    std::uint64_t positive_rank12=0,positive_rank34=0,positive_rank_other=0;
    std::vector<std::string> hitlist,near;
    const Q one(1),two(2),four(4);
    for(const RV&xx:vv){const Q&x=xx.q;
      for(const RV&yy:vv){const Q&y=yy.q;++pairs;
        if(eq(sqr(x),sqr(y)))continue;++distinct;
        Q z=add(x,y),den=add(one,z);if(!den.n)continue;++denok;
        Q S=neg(den),P=divide(mul(mul(z,add(one,x)),add(one,y)),den);
        Q D=sub(sqr(S),mul(four,P)),sd;if(D.n<=0)continue;++dpos;
        if(!sqrtq(D,sd))continue;++dsq;
        Q u=divide(add(S,sd),two),v=divide(sub(S,sd),two);
        if(!u.n||!v.n)continue;
        std::array<Q,5> r={one,x,y,u,v};bool sm=true;
        for(int i=0;i<5;++i)for(int j=i+1;j<5;++j)if(eq(sqr(r[i]),sqr(r[j])))sm=false;
        if(!sm)continue;++smoothck;
        Q p1(0),p3(0);for(const Q&w:r){p1=add(p1,w);p3=add(p3,mul(sqr(w),w));}
        if(p1.n||p3.n){std::cerr<<"CK check failed\n";return 3;}

        Q a2=sqr(x),a3=sqr(y),a4=sqr(u),a5=sqr(v);
        Q G0=neg(mul(mul(sub(one,a3),sub(one,a4)),sub(one,a5)));
        Q G3=mul(mul(sub(a3,a2),sub(one,a4)),sub(one,a5));
        Q G4=mul(mul(sub(a4,a2),sub(one,a3)),sub(one,a5));
        Q G5=mul(mul(sub(a5,a2),sub(one,a3)),sub(one,a4));
        if(G0.n<=0||G3.n<=0||G4.n<=0||G5.n<=0)continue;++gpositive;
        // Independent verification of the real ordering sieve: the marked
        // pair {a1,a2} must occupy sorted positions {1,2} or {3,4}.
        std::array<Q,5> av={one,a2,a3,a4,a5};
        auto lessq=[](const Q&A,const Q&B){return A.n*B.d<B.n*A.d;};
        std::sort(av.begin(),av.end(),lessq);
        int pos1=-1,pos2=-1;
        for(int z=0;z<5;++z){if(eq(av[z],one))pos1=z+1;if(eq(av[z],a2))pos2=z+1;}
        if(pos1>pos2)std::swap(pos1,pos2);
        if(pos1==1&&pos2==2)++positive_rank12;
        else if(pos1==3&&pos2==4)++positive_rank34;
        else ++positive_rank_other;
        bool b0=squareq(G0),b3=squareq(G3),b4=squareq(G4),b5=squareq(G5);
        if(b0)++g0sq;if(b3)++g3sq;if(b4)++g4sq;if(b5)++g5sq;
        if(near.size()<20&&(b0||b3||b4||b5)){
            std::ostringstream o;o<<"PARTIAL x="<<qs(x)<<" y="<<qs(y)<<" u="<<qs(u)<<" v="<<qs(v)
              <<" squares="<<b0<<","<<b3<<","<<b4<<","<<b5;near.push_back(o.str());
        }
        if(!(b0&&b3&&b4&&b5))continue;++hits;
        std::ostringstream o;o<<"HIT r=[1,"<<qs(x)<<","<<qs(y)<<","<<qs(u)<<","<<qs(v)<<"]";hitlist.push_back(o.str());
      }
    }
    std::cout<<"ELKIES22210_CK_ORBITB_DIRECT\nparameter_height "<<B
      <<"\nrational_values "<<vv.size()<<"\nordered_xy "<<pairs
      <<"\ndistinct_input_squares "<<distinct<<"\nnonzero_elimination_denominator "<<denok
      <<"\ndisc_positive "<<dpos<<"\ndisc_square "<<dsq<<"\nsmooth_ck "<<smoothck
      <<"\nall_G_positive "<<gpositive
      <<"\npositive_marked_pair_ranks_12_34_other "<<positive_rank12<<","<<positive_rank34<<","<<positive_rank_other
      <<"\nindividual_G_squares "<<g0sq<<","<<g3sq<<","<<g4sq<<","<<g5sq<<"\n";
    for(const auto&s:near)std::cout<<s<<"\n";
    std::cout<<"hits "<<hits<<"\n";for(const auto&s:hitlist)std::cout<<s<<"\n";std::cout<<"DONE\n";
}
