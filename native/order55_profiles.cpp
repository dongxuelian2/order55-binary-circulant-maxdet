// Folded profile enumeration only. This does NOT enumerate binary lifts.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>
using u64=std::uint64_t;
using u128=unsigned __int128;
constexpr u64 prime=2305843009213696591ULL;
u64 mul(u64 a,u64 b){return (u128)a*b%prime;}
u64 power(u64 a,int b){u64 r=1;while(b){if(b&1)r=mul(r,a);a=mul(a,a);b>>=1;}return r;}
u64 pow64(u64 a,u64 b){u64 r=1;while(b){if(b&1)r=mul(r,a);a=mul(a,a);b>>=1;}return r;}
int m,k,cap;
std::map<int,u64> thresholds;
std::array<int,11>a{};
u64 roots[11][11];
u64 nodes=0,leaves=0,canonical=0,surviving=0,covered=0;
void visit(int pos,int sum,int square){
  ++nodes;
  int left=m-pos;
  if(sum>k || sum+left*cap<k)return;
  if(!left){
    ++leaves;
    auto need=thresholds.find(square);
    if(need==thresholds.end())return;
    int stabilizer=0;
    for(int u=1;u<m;++u)for(int t=0;t<m;++t){
      int comparison=0;
      for(int j=0;j<m;++j){int b=a[(u*j+t)%m];if(b!=a[j]){comparison=b<a[j]?-1:1;break;}}
      if(comparison<0)return;
      if(!comparison)++stabilizer;
    }
    ++canonical;
    u64 norm=1;
    for(int j=1;j<m;++j){
      u64 z=0;for(int t=0;t<m;++t){z=(z+mul(a[t],roots[j][t]))%prime;}
      norm=mul(norm,z);
    }
    // All profile norms < (m*cap*cap)^((m-1)/2) < prime for m=5,11.
    if(norm<need->second)return;
    ++surviving;
    int orbit=m*(m-1)/stabilizer;covered+=orbit;
    std::cout<<"{\"m\":"<<m<<",\"k\":"<<k<<",\"squares\":"<<square
       <<",\"norm\":"<<norm<<",\"orbit_size\":"<<orbit<<",\"profile\":[";
    for(int j=0;j<m;++j){if(j)std::cout<<',';std::cout<<a[j];}
    std::cout<<"]}\n";
    return;
  }
  int remaining=k-sum;
  int q=remaining/left,r=remaining%left;
  int min_square=square+(left-r)*q*q+r*(q+1)*(q+1);
  if(min_square>thresholds.rbegin()->first || square+remaining*cap<thresholds.begin()->first)return;
  int lo=std::max(0,remaining-(left-1)*cap), hi=std::min(cap,remaining);
  if(pos)lo=std::max(lo,a[0]); // every orbit has a representative with min at zero
  for(int v=lo;v<=hi;++v){a[pos]=v;visit(pos+1,sum+v,square+v*v);}
}
int main(int argc,char**argv){
 if(argc!=4)return 2;
 m=std::stoi(argv[2]);k=std::stoi(argv[3]);if((m!=5&&m!=11)||k<1||k>54)return 2;cap=55/m;
 std::ifstream in(argv[1]);int mm,kk,ss;u64 need;
 while(in>>mm>>kk>>ss>>need)if(mm==m&&kk==k)thresholds[ss]=need;
 if(thresholds.empty())return 0;
 u64 w=1;for(u64 b=2;w==1;++b)w=pow64(b,(prime-1)/m);
 for(int j=0;j<m;++j)for(int t=0;t<m;++t)roots[j][t]=power(w,j*t);
 auto start=std::chrono::steady_clock::now();visit(0,0,0);
 std::cerr<<"m="<<m<<" k="<<k<<" nodes="<<nodes<<" min_at_zero_leaves="<<leaves
 <<" threshold_canonical="<<canonical<<" survivors="<<surviving<<" covered_profiles="<<covered
 <<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
}
