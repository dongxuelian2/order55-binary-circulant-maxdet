// Exact small-order SIGN regression: direct Gray-code sign spectrum,
// independent multiset profile traversal, and separate binary split lifts.
#include <algorithm>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
using U=std::uint64_t;using W=unsigned __int128;
constexpr U P=(U(1)<<61)-1;
U mul(U a,U b){return W(a)*b%P;}
U power(U a,U b){U r=1;while(b){if(b&1)r=mul(r,a);a=mul(a,a);b>>=1;}return r;}
U root(int n){for(U a=2;;++a){U z=power(a,(P-1)/n);bool ok=power(z,n)==1;for(int d=1;d<n;++d)if(n%d==0&&power(z,d)==1)ok=false;if(ok)return z;}}
std::string word(U x,int n){std::string s;for(int i=0;i<n;++i)s+=char('0'+((x>>i)&1));return s;}
U rotate(U x,int j,int n){return ((x<<j)|(x>>(n-j)))&((U(1)<<n)-1);}
std::string profile(U x,int n){std::string s;s+=char(std::popcount(x));for(int j=1;j<=n/2;++j)s+=char(std::popcount(x&rotate(x,j,n)));return s;}
int main(int argc,char**argv){
 if(argc<4)return 2;std::string mode=argv[1];int n=std::stoi(argv[2]);int h=n/2;
 if(n!=15&&n!=21)return 3;
 U z=root(n);std::vector<std::vector<U>>t(n,std::vector<U>(n));
 for(int i=0;i<n;++i)for(int j=0;j<n;++j)t[i][j]=power(z,i*j);
 if(mode=="brute"){
  std::vector<U>f(n),maxk(n+1);f[0]=P-n;U best=0,mask=0;std::vector<U>winners;
  for(U step=0;step<(U(1)<<n);++step){
   if(step){int i=std::countr_zero(step);U next=step^(step>>1);bool on=next>>i&1;for(int j=0;j<n;++j){U delta=2*t[i][j]%P;f[j]=on?(f[j]+delta)%P:(f[j]+P-delta)%P;}mask=next;}
   U det=1;for(U v:f)det=mul(det,v);det=std::min(det,P-det);
   int k=std::popcount(mask);maxk[k]=std::max(maxk[k],det);
   if(det>best){best=det;winners.clear();}if(det==best)winners.push_back(mask);
  }
  std::ofstream out(argv[3]);out<<"{\"n\":"<<n<<",\"all_words\":"<<(U(1)<<n)<<",\"raw\":"<<best<<",\"normalized\":"<<(best>>(n-1))<<",\"max_by_weight\":[";
  for(int k=0;k<=n;++k){if(k)out<<',';out<<maxk[k];}out<<"],\"winners\":[";
  for(unsigned i=0;i<winners.size();++i){if(i)out<<',';out<<'\"'<<word(winners[i],n)<<'\"';}out<<"]}\n";
 }else if(mode=="profiles"){
  if(argc!=6)return 4;std::ifstream in(argv[3]);std::ofstream out(argv[4]);U threshold=std::stoull(argv[5]),visited=0,survived=0;int k;
  while(in>>k){std::vector<int>c(h);for(auto&v:c)in>>v;if(!in)return 5;
   do {++visited;U v=n-2*k;for(int j=1;j<=h;++j){U q=k;for(int s=1;s<=h;++s)q=(q+mul(c[s-1],(t[j][s]+t[j][n-s])%P))%P;v=mul(v,q);}v=std::min(v,P-v);
    if(v>=threshold){++survived;out<<k<<' '<<v;for(int c0:c)out<<' '<<c0;out<<'\n';}
   }while(std::next_permutation(c.begin(),c.end()));
  }
  std::cerr<<"profiles "<<n<<' '<<visited<<' '<<survived<<'\n';
 }else if(mode=="lifts"){
  if(argc!=5)return 6;std::ifstream in(argv[3]);std::ofstream out(argv[4]);
  std::unordered_map<std::string,std::pair<U,U>>targets;int k;U value;
  while(in>>k>>value){std::string key(1,char(k));for(int j=1;j<=h;++j){int c;in>>c;key+=char(c);}targets[key]={value,0};}
  // Disjoint halves of CRT columns; every low-weight bit word occurs once.
  int p=3,q=n/3,split=q/2;std::vector<int>left,right;
  for(int col=0;col<q;++col)for(int row=0;row<p;++row) {
   int index=0;while(index%p!=row||index%q!=col)++index;
   (col<split?left:right).push_back(index);
  }
  auto expand=[](U x,const std::vector<int>&v){U b=0;for(unsigned i=0;i<v.size();++i)if(x>>i&1)b|=U(1)<<v[i];return b;};
  std::vector<U>L,R;for(U a=0;a<(U(1)<<left.size());++a)L.push_back(expand(a,left));
  for(U a=0;a<(U(1)<<right.size());++a)R.push_back(expand(a,right));
  U considered=0,solutions=0;
  for(U a:L)for(U b:R){U w=a|b;int weight=std::popcount(w);if(!weight||weight>h)continue;++considered;auto key=profile(w,n);auto it=targets.find(key);if(it==targets.end())continue;++it->second.second;++solutions;out<<word(w,n)<<' '<<it->second.first<<'\n';}
  U realized=0;for(const auto&entry:targets)realized+=entry.second.second>0;
  std::cerr<<"lifts "<<n<<" all reduced "<<considered<<" target profiles "<<targets.size()<<" realized "<<realized<<" solutions "<<solutions<<'\n';
 }else return 7;
}

