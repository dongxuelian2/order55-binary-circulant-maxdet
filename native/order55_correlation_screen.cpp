// Exact finite-field determinant of every configured integer correlation profile.
// A profile is an overapproximation: binary realizability is NOT asserted.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
using u64=std::uint64_t;using u128=unsigned __int128;
constexpr u64 P=2305843009213696591ULL,Q=2305843009213697141ULL;
u64 mul(u64 a,u64 b,u64 p){return (u128)a*b%p;}
u64 power(u64 a,u64 b,u64 p){u64 r=1;while(b){if(b&1)r=mul(r,a,p);a=mul(a,a,p);b>>=1;}return r;}
std::string decimal(u128 x){std::string s;do{s+=char('0'+x%10);x/=10;}while(x);std::reverse(s.begin(),s.end());return s;}
u128 parse(std::string s){u128 x=0;for(char c:s)x=10*x+c-'0';return x;}
int k,base,values[27],freqs[27],counts[28];
u64 table[2][27][28],p[2]={P,Q};
std::set<std::array<int,3>> allowed5;
std::set<std::array<int,6>> allowed11;
u64 leaves=0,folded_pass=0,det_tests=0,above=0; std::int64_t coslo[28],coshi[28]; bool use_intervals=false; u128 product_threshold;
u128 best=0,incumbent;u64 inv=power(P,Q-2,Q);
std::array<int,27> winner{};
void evaluate(){
 ++leaves;
 std::array<int,3> h5{};h5[0]=k;
 std::array<int,6> h11{};h11[0]=k;
 for(int s=1;s<=27;++s){int a=s%5,b=s%11;h5[std::min(a,5-a)]+=values[s-1]*(a?1:2);h11[std::min(b,11-b)]+=values[s-1]*(b?1:2);}
 // For nonzero residues, pairing aggregates one term into each symmetric slot.
 if(!allowed5.count(h5)||!allowed11.count(h11))return;
 ++folded_pass;
 // Unit action on the 27 independent +/- shifts has effective group order 20.
 for(int u=2;u<28;++u){if(u%5==0||u%11==0)continue;
   for(int s=1;s<=27;++s){int t=s*u%55;t=std::min(t,55-t);if(values[t-1]<values[s-1])return;if(values[t-1]>values[s-1])break;}
 }
 ++det_tests;
 if(use_intervals){
   u128 upper=(u128)1<<64;
   for(int j=1;j<=27;++j){
     std::int64_t q=(std::int64_t)(k-base)<<24;
     for(int s=1;s<=27;++s){int d=values[s-1]-base,t=j*s%55;t=std::min(t,55-t);q+=d*(d>=0?coshi[t]:coslo[t]);}
     if(q<=0)return;
     upper=(upper*(u64)q+(((u128)1<<28)-1))>>28;
   }
   if(upper<product_threshold)return;
 }
 u64 residue[2];
 for(int field=0;field<2;++field){u64 z=55-k;
   for(int j=0;j<27;++j){u64 e=(k-base+p[field])%p[field];
     for(int s=1;s<=27;++s){int d=values[s-1]-base;
       if(d>0)e=(e+mul(d,table[field][j][s],p[field]))%p[field];
       else if(d<0)e=(e+p[field]-mul(-d,table[field][j][s],p[field]))%p[field];
     }z=mul(z,e,p[field]);
   }residue[field]=z;
 }
 u64 delta=(residue[1]+Q-residue[0])%Q;
 u128 value=residue[0]+(u128)P*mul(delta,inv,Q);
 if(value>(u128)P*Q/2)value=(u128)P*Q-value;
 if(value>best){best=value;std::copy(values,values+27,winner.begin());}
 if(value>=incumbent){++above;std::cout<<"{\"k\":"<<k<<",\"absolute_profile_product\":\""<<decimal(value)<<"\",\"correlations\":[";
 for(int i=0;i<27;++i){if(i)std::cout<<',';std::cout<<values[i];}std::cout<<"]}\n";}
}
void visit(int pos){if(pos==27){evaluate();return;}for(int v=0;v<=k;++v)if(counts[v]){--counts[v];values[pos]=v;visit(pos+1);++counts[v];}}
int main(int argc,char**argv){
 if(argc!=6&&argc!=7)return 2;k=std::stoi(argv[1]);base=k*(k-1)/54;incumbent=parse(argv[2]);
 if(argc==7){std::ifstream ci(argv[6]);for(int i=0;i<28;++i)ci>>coslo[i]>>coshi[i];if(!ci)return 5;use_intervals=true;
  u128 denominator=(u128)(55-k)<<44;product_threshold=(incumbent+denominator-1)/denominator;}
 for(int m:{5,11}){std::ifstream f(std::string(argv[4])+std::to_string(m)+"_k"+std::to_string(k)+".txt");std::string line;
 while(std::getline(f,line)){std::istringstream in(line);if(m==5){std::array<int,3>x;for(auto&v:x)in>>v;allowed5.insert(x);}else{std::array<int,6>x;for(auto&v:x)in>>v;allowed11.insert(x);}}}
 if(allowed5.empty()||allowed11.empty())return 3;
 for(int f=0;f<2;++f){u64 w=1;for(u64 b=2;;++b){w=power(b,(p[f]-1)/55,p[f]);if(power(w,5,p[f])!=1&&power(w,11,p[f])!=1)break;}
 for(int j=1;j<=27;++j)for(int s=1;s<=27;++s)table[f][j-1][s]=(power(w,j*s%55,p[f])+power(w,(55-j*s%55)%55,p[f]))%p[f];}
 auto start=std::chrono::steady_clock::now();std::ifstream parts(argv[3]);std::string line;while(std::getline(parts,line)){std::fill(counts,counts+28,0);std::istringstream in(line);int v,n=0;while(in>>v){++counts[v];++n;}if(n!=27)return 4;visit(0);}
 std::ofstream audit(argv[5]);audit<<"{\"k\":"<<k<<",\"ordered_profiles\":"<<leaves<<",\"folded_matching_profiles\":"<<folded_pass
 <<",\"canonical_determinants\":"<<det_tests<<",\"above_screen\":"<<above<<",\"max_absolute_profile_product\":\""<<decimal(best)<<"\",\"seconds\":"
 <<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
}
