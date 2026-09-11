// Heuristic discovery only: floating point never prunes a proof branch.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <numeric>
#include <random>
#include <string>
#include <thread>
#include <vector>
using A=std::array<double,27>;
A re[55],im[55];
std::vector<std::string> seeds;
std::mutex guard;
std::string output;
double seconds;
int firstk,lastk;
void search(int k,int replica) {
 std::mt19937_64 gen(551109ULL+1009*k+replica*1000003ULL);
 std::uniform_real_distribution<double> uniform(0,1);
 std::array<int,55>a{};
 A r{},s{};
 auto rebuild=[&](){r.fill(0);s.fill(0);for(int i=0;i<55;++i)if(a[i])for(int j=0;j<27;++j){r[j]+=re[i][j];s[j]+=im[i][j];}};
 auto score=[&](int u,int v,int w=-1,int z=-1){double q=55-2*k;for(int j=0;j<27;++j){double rr=r[j]+re[v][j]-re[u][j],ss=s[j]+im[v][j]-im[u][j];if(w>=0){rr+=re[z][j]-re[w][j];ss+=im[z][j]-im[w][j];}q*=rr*rr+ss*ss;}return q;};
 auto current=[&](){double q=55-2*k;for(int j=0;j<27;++j)q*=r[j]*r[j]+s[j]*s[j];return q;};
 auto swap=[&](int u,int v){a[u]=0;a[v]=1;for(int j=0;j<27;++j){r[j]+=re[v][j]-re[u][j];s[j]+=im[v][j]-im[u][j];}};
 double best=-1;
 std::array<int,55> champion{};
 long starts=0,moves=0;
 auto start=std::chrono::steady_clock::now();
 auto record=[&](double value){if(value>best*(1+1e-12)){best=value;champion=a;std::string word;for(int b:a)word+=char('0'+b);
   std::lock_guard<std::mutex> lock(guard);std::ofstream out(output,std::ios::app);
   out<<std::setprecision(17)<<"{\"k\":"<<k<<",\"replica\":"<<replica<<",\"seed\":"<<551109ULL+1009*k+replica*1000003ULL<<",\"word\":\""<<word<<"\",\"score\":"<<value<<",\"starts\":"<<starts<<",\"moves\":"<<moves<<"}\n";
 }};
 while(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<seconds) {
  ++starts;
  if(starts%5==0&&best>0)a=champion;
  else if(starts<=long(seeds.size())) {auto word=seeds[starts-1];for(int i=0;i<55;++i)a[i]=word[i]-'0';int w=std::accumulate(a.begin(),a.end(),0);while(w!=k){int i=gen()%55;if(w>k&&a[i]){a[i]=0;--w;}else if(w<k&&!a[i]){a[i]=1;++w;}}}
  else {a.fill(0);std::array<int,55>perm;std::iota(perm.begin(),perm.end(),0);std::shuffle(perm.begin(),perm.end(),gen);for(int i=0;i<k;++i)a[perm[i]]=1;}
  rebuild();double value=current();
  // Annealing with fixed-weight moves; temperatures act on log objective.
  for(int t=0;t<4500;++t) {
   int u,v;do{u=gen()%55;}while(!a[u]);do{v=gen()%55;}while(a[v]);
   double next=score(u,v),temp=0.35*std::pow(0.003/0.35,double(t)/4500);
   if(next>value|| (next>0&&uniform(gen)<std::pow(next/value,1/temp))) {swap(u,v);value=next;++moves;}
  }
  // Complete best one-swap descent.
  for(;;){int u0=-1,v0=-1;double next=value*(1+1e-12);for(int u=0;u<55;++u)if(a[u])for(int v=0;v<55;++v)if(!a[v]){double q=score(u,v);if(q>next){next=q;u0=u;v0=v;}}if(u0<0)break;swap(u0,v0);value=next;++moves;}
  record(value);
  // Complete two-swap on new or competitive local optima.
  if(value>best*0.995) {
   int us[2]={-1,-1},vs[2]={-1,-1};double next=value*(1+1e-12);
   for(int u=0;u<55;++u)if(a[u])for(int w=u+1;w<55;++w)if(a[w])
    for(int v=0;v<55;++v)if(!a[v])for(int z=v+1;z<55;++z)if(!a[z]){double q=score(u,v,w,z);if(q>next){next=q;us[0]=u;us[1]=w;vs[0]=v;vs[1]=z;}}
   if(us[0]>=0){swap(us[0],vs[0]);swap(us[1],vs[1]);rebuild();record(current());}
  }
 }
 std::lock_guard<std::mutex> lock(guard);std::cerr<<"k="<<k<<" replica="<<replica<<" starts="<<starts<<" best="<<std::setprecision(17)<<best<<"\n";
}
int main(int argc,char**argv) {
 if(argc!=7)return 2;
 firstk=std::stoi(argv[1]);lastk=std::stoi(argv[2]);seconds=std::stod(argv[3]);int replicas=std::stoi(argv[4]);output=argv[5];
 if(firstk<1||lastk>27||firstk>lastk||replicas<1)return 3;
 std::ifstream in(argv[6]);std::string word;while(in>>word)if(word.size()==55&&word.find_first_not_of("01")==std::string::npos)seeds.push_back(word);
 for(int i=0;i<55;++i)for(int j=0;j<27;++j){double angle=2*std::acos(-1.0)*i*(j+1)/55;re[i][j]=std::cos(angle);im[i][j]=std::sin(angle);}
 std::vector<std::thread>threads;for(int k=firstk;k<=lastk;++k)for(int replica=0;replica<replicas;++replica)threads.emplace_back(search,k,replica);
 for(auto&t:threads)t.join();
}

