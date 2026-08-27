#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

struct Profile {
    int p=0,o1=0,o2=0;
    std::vector<unsigned char> ok;
    long long count=0;
    inline bool allows(int ti,long long m,long long n) const {
        int a=(int)(m%o1); if(a<0)a+=o1;
        int b=(int)(n%o2); if(b<0)b+=o2;
        return ok[((ti*o1)+a)*o2+b];
    }
};
struct Lattice {int ti; long long rm,mm,rn,mn;};

static std::vector<std::string> split(const std::string& s){
    std::vector<std::string> a; std::stringstream ss(s);std::string x;
    while(std::getline(ss,x,'\t'))a.push_back(x);return a;
}
static long long floordiv(long long a,long long b){
    long long q=a/b,r=a%b;if(r && ((r>0)!=(b>0)))--q;return q;
}
static long long ceildiv(long long a,long long b){return -floordiv(-a,b);}

int main(int argc,char**argv){
    long long N=argc>1?std::stoll(argv[1]):1000000;
    std::string lattice=argc>2?argv[2]:"results/target_22224_n11_rank2_periodic_lattice_filtered2.tsv";
    std::string profiles=argc>3?argv[3]:"results/target_22224_n11_rank2_profiles_p199_all.tsv";
    std::string output=argc>4?argv[4]:"results/target_22224_n11_rank2_modular_N1000000_p199.tsv";
    const std::unordered_set<int> used={29,37,41,101,107,137,191};
    std::unordered_map<int,Profile> pm;
    std::ifstream pf(profiles);std::string line;std::getline(pf,line);
    while(std::getline(pf,line)){
        auto a=split(line);if(a.size()<8)continue;
        int p=std::stoi(a[0]);if(used.count(p))continue;
        int ti=std::stoi(a[1]),m=std::stoi(a[2]),n=std::stoi(a[3]);
        int o1=std::stoi(a[4]),o2=std::stoi(a[5]);
        auto &q=pm[p];if(!q.p){q.p=p;q.o1=o1;q.o2=o2;q.ok.assign(5*o1*o2,0);}
        auto idx=((ti*o1)+m)*o2+n;if(!q.ok[idx]){q.ok[idx]=1;q.count++;}
    }
    std::vector<Profile> ps;for(auto &kv:pm)ps.push_back(std::move(kv.second));
    std::sort(ps.begin(),ps.end(),[](const Profile&a,const Profile&b){
        return (__int128)a.count*b.o1*b.o2 < (__int128)b.count*a.o1*a.o2;
    });
    std::vector<Lattice> ls;std::ifstream lf(lattice);std::getline(lf,line);
    while(std::getline(lf,line)){auto a=split(line);if(a.size()<5)continue;
        ls.push_back({std::stoi(a[0]),std::stoll(a[1]),std::stoll(a[2]),std::stoll(a[3]),std::stoll(a[4])});}
    std::ofstream out(output);out<<"m\tn\ttorsion_coset\n";
    unsigned long long tested=0,surv=0;
    for(const auto&c:ls){
        long long im0=ceildiv(-N-c.rm,c.mm),im1=floordiv(N-c.rm,c.mm);
        long long in0=ceildiv(-N-c.rn,c.mn),in1=floordiv(N-c.rn,c.mn);
        for(long long im=im0;im<=im1;im++){
            long long m=c.rm+c.mm*im;
            for(long long in=in0;in<=in1;in++){
                long long n=c.rn+c.mn*in;tested++;bool good=true;
                for(const auto&q:ps)if(!q.allows(c.ti,m,n)){good=false;break;}
                if(good){surv++;out<<m<<'\t'<<n<<'\t'<<c.ti<<'\n';}
            }
        }
    }
    std::cout<<"bound="<<N<<" lattices="<<ls.size()<<" profiles="<<ps.size()
             <<" tested="<<tested<<" survivors="<<surv<<" output="<<output<<"\n";
}
