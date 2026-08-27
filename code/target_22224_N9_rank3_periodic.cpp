// Fast periodic enumeration of the N9 rank-3 Mordell--Weil lattice.
//
// For fixed (m,n,torsion coset), selected local profiles give short lists
// of admissible k residues.  Generalized CRT combines those lists into
// classes modulo lcm(ord(G3 mod p)); only those k are enumerated.  Every
// resulting triple is then checked against every exported Magma profile.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using i64 = long long;

static std::vector<std::string> split(const std::string &s, char c) {
    std::vector<std::string> z; std::string t; std::stringstream ss(s);
    while (std::getline(ss,t,c)) z.push_back(t); return z;
}
static i64 mod(i64 a,i64 m){a%=m; if(a<0)a+=m; return a;}
static i64 egcd(i64 a,i64 b,i64 &x,i64 &y){
    if(!b){x=1;y=0;return a;} i64 u,v,g=egcd(b,a%b,u,v); x=v; y=u-(a/b)*v; return g;
}
static bool crt(i64 a,i64 m,i64 b,i64 n,i64 &r,i64 &M){
    i64 x,y,g=egcd(m,n,x,y),d=b-a; if(d%g) return false;
    i64 ng=n/g,k=ng==1?0:mod((d/g)*x,ng); M=(m/g)*n; r=mod(a+m*k,M); return true;
}
static i64 floordiv(i64 a,i64 b){
    i64 q=a/b,r=a%b; if(r && ((r>0)!=(b>0))) --q; return q;
}
static i64 ceildiv(i64 a,i64 b){return -floordiv(-a,b);}

struct Profile {
    int p=0,o1=0,o2=0,o3=0,nt=0;
    std::vector<unsigned char> ok;
    size_t ix(int ti,i64 m,i64 n,i64 k) const {
        return (((size_t)(ti-1)*o1+(size_t)mod(m,o1))*o2+(size_t)mod(n,o2))*o3+(size_t)mod(k,o3);
    }
    bool get(int ti,i64 m,i64 n,i64 k) const {return ok[ix(ti,m,n,k)]!=0;}
};

struct Opt {
    int N=100;
    std::string input="results/target_22224_N9_rank3_residue_classes.tsv";
    std::string output="results/target_22224_N9_rank3_periodic_N100.tsv";
    std::vector<int> crtps={67,89,43,61,13};
};
static Opt parse(int ac,char **av){
    Opt o;
    for(int i=1;i<ac;i++){
        std::string a=av[i]; auto need=[&](){if(++i>=ac)throw std::runtime_error("missing option");return std::string(av[i]);};
        if(a=="--N")o.N=std::stoi(need());
        else if(a=="--input")o.input=need();
        else if(a=="--output")o.output=need();
        else if(a=="--crt-primes"){
            o.crtps.clear(); for(auto&s:split(need(),','))o.crtps.push_back(std::stoi(s));
        } else throw std::runtime_error("unknown option "+a);
    }
    return o;
}

int main(int ac,char **av){
 try {
    Opt o=parse(ac,av); std::ifstream in(o.input); if(!in)throw std::runtime_error("cannot open input");
    std::string line; std::getline(in,line); std::map<int,Profile> byp; int maxTi=0; uint64_t rows=0;
    while(std::getline(in,line)){
        if(line.empty())continue; auto z=split(line,'\t'); if(z.size()<9)continue;
        int p=std::stoi(z[0]),ti=std::stoi(z[1]),m=std::stoi(z[2]),n=std::stoi(z[3]),k=std::stoi(z[4]);
        int o1=std::stoi(z[5]),o2=std::stoi(z[6]),o3=std::stoi(z[7]); maxTi=std::max(maxTi,ti);
        auto it=byp.find(p);
        if(it==byp.end()){
            Profile P;P.p=p;P.o1=o1;P.o2=o2;P.o3=o3;P.nt=2;
            P.ok.assign((size_t)2*o1*o2*o3,0); it=byp.emplace(p,std::move(P)).first;
        }
        if(ti>it->second.nt)throw std::runtime_error("unexpected torsion coset");
        it->second.ok[it->second.ix(ti,m,n,k)]=1; ++rows;
    }
    if(maxTi!=2)throw std::runtime_error("expected two torsion cosets");
    std::vector<const Profile*> crtprof,allprof;
    for(auto &[p,P]:byp)allprof.push_back(&P);
    for(int p:o.crtps){auto it=byp.find(p);if(it==byp.end())throw std::runtime_error("missing CRT prime "+std::to_string(p));crtprof.push_back(&it->second);}
    std::ofstream out(o.output); if(!out)throw std::runtime_error("cannot open output"); out<<"m\tn\tk\ttorsion_coset\n";
    uint64_t pairs=0,empty=0,crtclasses=0,kenum=0,surv=0;
    std::cout<<"N9_RANK3_PERIODIC_START N="<<o.N<<" profiles="<<byp.size()<<" rows="<<rows<<" crt_primes=";
    for(int p:o.crtps)std::cout<<p<<","; std::cout<<"\n";
    for(i64 m=-o.N;m<=o.N;m++) for(i64 n=-o.N;n<=o.N;n++) for(int ti=1;ti<=2;ti++){
        ++pairs; std::vector<i64> cls={0}; i64 M=1;
        for(const Profile *P:crtprof){
            std::vector<int> ks; for(int k=0;k<P->o3;k++)if(P->get(ti,m,n,k))ks.push_back(k);
            if(ks.empty()){cls.clear();break;}
            std::vector<i64> nxt; i64 nextM=std::lcm<i64>(M,P->o3);
            for(i64 a:cls)for(int b:ks){i64 r,L;if(crt(a,M,b,P->o3,r,L))nxt.push_back(r);}
            std::sort(nxt.begin(),nxt.end());nxt.erase(std::unique(nxt.begin(),nxt.end()),nxt.end());cls.swap(nxt);M=nextM;
            if(cls.empty())break;
        }
        if(cls.empty()){++empty;continue;} crtclasses+=cls.size();
        for(i64 r:cls){
            i64 lo=ceildiv(-o.N-r,M),hi=floordiv(o.N-r,M);
            for(i64 q=lo;q<=hi;q++){
                i64 k=r+M*q;++kenum;bool good=true;
                for(const Profile *P:allprof)if(!P->get(ti,m,n,k)){good=false;break;}
                if(!good)continue;++surv;out<<m<<'\t'<<n<<'\t'<<k<<'\t'<<ti<<'\n';
            }
        }
    }
    std::cout<<"N9_RANK3_PERIODIC_DONE mn_cosets="<<pairs<<" empty="<<empty
             <<" crt_classes="<<crtclasses<<" k_enumerated="<<kenum
             <<" survivors="<<surv<<" output="<<o.output<<"\n";
    return 0;
 } catch(const std::exception&e){std::cerr<<"ERROR "<<e.what()<<"\n";return 2;}
}
