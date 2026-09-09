// Exact 11 columns x 5 bits MITM for fixed margins and full autocorrelation.
#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using u64=std::uint64_t;
constexpr u64 mask55=(1ULL<<55)-1;
struct Entry{std::uint32_t key;u64 word;bool operator<(const Entry&b)const{return key<b.key;}};
std::vector<Entry> table;
int r[5],s[11],target[55],partial[7]{};
u64 nodes=0,left_count=0,right_count=0,joined=0,solutions=0;
int split; std::string task_id;
std::uint32_t key(bool complement){std::uint32_t v=0;for(int i=0;i<7;++i){int n=partial[i];if(complement)n=(i<5?r[i]:target[(i-4)*11])-n;if(n<0||n>15)return UINT32_MAX;v|=n<<(4*i);}return v;}
void visit(int col,int end,u64 word,bool left){
 ++nodes;
 if(col==end){
  auto sig=key(!left);if(sig==UINT32_MAX)return;
  if(left){table.push_back({sig,word});++left_count;return;}
  ++right_count;
  auto it=std::lower_bound(table.begin(),table.end(),Entry{sig,0});
  for(;it!=table.end()&&it->key==sig;++it){
   ++joined;u64 a=word|it->word;bool ok=true;
   for(int shift=1;shift<=27;++shift){if(shift==11||shift==22)continue;
    auto rotated=((a<<shift)|(a>>(55-shift)))&mask55;
    if(std::popcount(a&rotated)!=target[shift]){ok=false;break;}
   }
   if(ok){++solutions;std::string bits(55,'0');for(int i=0;i<55;++i)if(a>>i&1)bits[i]='1';std::cout<<task_id<<' '<<bits<<'\n';}
  }return;
 }
 for(unsigned mask=0;mask<32;++mask){if(std::popcount(mask)!=s[col])continue;
  bool ok=true;for(int i=0;i<5;++i)if(partial[i]+int(mask>>i&1)>r[i])ok=false;if(!ok)continue;
  int a=std::popcount(mask&(((mask<<1)|(mask>>4))&31));
  int b=std::popcount(mask&(((mask<<2)|(mask>>3))&31));
  if(partial[5]+a>target[11]||partial[6]+b>target[22])continue;
  u64 updated=word;
  for(int i=0;i<5;++i){partial[i]+=mask>>i&1;if(mask>>i&1)updated|=1ULL<<((11*i+45*col)%55);}
  partial[5]+=a;partial[6]+=b;visit(col+1,end,updated,left);
  for(int i=0;i<5;++i)partial[i]-=mask>>i&1;partial[5]-=a;partial[6]-=b;
 }
}
int main(int argc,char**argv){
 if(argc!=4&&argc!=5)return 2;
 bool batch=argc==5;std::ifstream in(argv[1]);split=std::stoi(argv[2]);if(split<1||split>10)return 3;
 std::ofstream audit(argv[3]);
 do{
  if(batch){if(!(in>>task_id))break;}else task_id="single";
  for(auto&v:r)in>>v;for(auto&v:s)in>>v;for(auto&v:target)in>>v;if(!in)return 4;
  table.clear();nodes=left_count=right_count=joined=solutions=0;std::fill(partial,partial+7,0);
  auto start=std::chrono::steady_clock::now();visit(0,split,0,true);std::sort(table.begin(),table.end());visit(split,11,0,false);
  audit<<"{\"task\":\""<<task_id<<"\",\"split\":"<<split<<",\"nodes\":"<<nodes<<",\"left_entries\":"<<left_count<<",\"right_entries\":"<<right_count
  <<",\"joined_words\":"<<joined<<",\"solutions\":"<<solutions<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
  audit.flush();std::cout.flush();
 }while(batch);
}