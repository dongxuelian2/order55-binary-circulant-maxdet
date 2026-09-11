// Independent fold enumeration by bounded-alphabet FKM necklaces.
// Unlike production, there is no affine-minimum test; every rotation class
// is generated, so all decimation orientations are already present.
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
using U=std::uint64_t;using W=unsigned __int128;
constexpr U P=2305843009213697141ULL;
U mul(U a,U b){return W(a)*b%P;}U powmod(U a,U e){U r=1;while(e){if(e&1)r=mul(r,a);a=mul(a,a);e>>=1;}return r;}
int m,k,cap,a[12];U count=0,survived=0;
std::map<int,U>threshold;U table[12][12];
void necklace(int t,int period,int total,int square){
 int left=m-t+1,remaining=k-total;if(remaining<0||remaining>left*cap)return;
 if(left){int b=remaining/left,r=remaining%left;if(square+(left-r)*b*b+r*(b+1)*(b+1)>threshold.rbegin()->first)return;}
 if(t>m){if(m%period||total!=k)return;++count;auto it=threshold.find(square);if(it==threshold.end())return;
  U norm=1;for(int f=1;f<m;++f){U v=0;for(int i=1;i<=m;++i)v=(v+mul(a[i],table[f][i-1]))%P;norm=mul(norm,v);}
  if(norm<it->second)return;++survived;std::cout<<"{\"m\":"<<m<<",\"k\":"<<k<<",\"norm\":"<<norm<<",\"profile\":[";
  for(int i=1;i<=m;++i){if(i>1)std::cout<<',';std::cout<<a[i];}std::cout<<"]}\n";return;
 }
 a[t]=a[t-period];necklace(t+1,period,total+a[t],square+a[t]*a[t]);
 for(int v=a[t-period]+1;v<=cap;++v){a[t]=v;necklace(t+1,t,total+v,square+v*v);}
}
int main(int argc,char**argv){
 if(argc!=5)return 2;m=std::stoi(argv[2]);k=std::stoi(argv[3]);cap=std::stoi(argv[4]);if(m<2||m>11||(P-1)%m)return 3;
 std::ifstream in(argv[1]);int mm,kk,square;U need;while(in>>mm>>kk>>square>>need)if(mm==m&&kk==k)threshold[square]=need;
 if(threshold.empty())return 0;U root=1;for(U b=2;root==1;++b)root=powmod(b,(P-1)/m);
 for(int f=0;f<m;++f)for(int i=0;i<m;++i)table[f][i]=powmod(root,f*i);
 necklace(1,1,0,0);std::cerr<<"necklaces="<<count<<" survived="<<survived<<'\n';
}

