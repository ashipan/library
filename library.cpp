//template
#include <bits/stdc++.h>
#include<boost/multiprecision/cpp_int.hpp>
#include<boost/multiprecision/cpp_dec_float.hpp>
namespace mp=boost::multiprecision;
#define mulint mp::cpp_int
#define mulfloat mp::cpp_dec_float_100
#define rep(i, a, n) for(int i = a; i < n; i++)
using namespace std;
using ll = long long;
using P = pair<int, int>;

const int INF = 1001001001;
const ll LINF = 1001002003004005006ll;
const int mod = 1000000007;
const int mod = 998244353;

//Aho-Corasick
struct Aho {
  vector<unordered_map<char, int>> to;
  vector<int> cnt, fail;
  Aho(): to(1), cnt(1) {}
  int add(const string& s) {
    int v = 0;
    for (char c : s) {
      if (!to[v].count(c)) {
        to[v][c] = (int)to.size();
        to.push_back(unordered_map<char, int>());
        cnt.push_back(0);
      }
      v = to[v][c];
    }
    cnt[v]++;
    return v;
  }
  void init() {
    fail = vector<int>(to.size(), -1);
    queue<int> q;
    q.push(0);
    while (!q.empty()) {
      int v = q.front(); q.pop();
      for (auto [c, u] : to[v]) {
        fail[u] = (*this)(fail[v], c);
        cnt[u] += cnt[fail[u]];
        q.push(u);
      }
    }
  }
  int operator()(int v, char c) const {
    while (v != -1) {
      auto it = to[v].find(c);
      if (it != to[v].end()) return it->second;
      v = fail[v];
    }
    return 0;
  }
  int operator[](int v) const { return cnt[v];}
};


//Array
template<typename T>
struct Array{
  int n;
  int maximum_count = 0, minimul_count = 0;
  vector<T> a;
  vector<int> updown;
  Array(vector<T>& _a): n((int)_a.size()), updown(n) { init(_a);}
  void init(vector<T>& _a){
    a.resize((int)_a.size());
    a = _a;
  }
  void differential(){
    vector<T> b;
    vector<int> res;
    int id = 1, state = 0;
    for(int i = 0; i < (int)a.size(); i++) if(i == 0 || a[i] != a[i-1]) b.emplace_back(a[i]);
    res.resize((int)b.size());
    for(int i = 1; i < n-1; i++){
      if(b[i-1] < b[i]){
        if(b[i] < b[i+1]) res[i] = 1;
        else { res[i] = 2; maximum_count++;}
      }
      else{
        if(b[i] > b[i+1]) res[i] = -1;
        else { res[i] = -2; minimul_count++;}
      }
    }
    for(int i = 1; i < n-1; i++){
      if(a[i] == a[i-1]) updown[i] = state;
      else { state = res[id]; updown[i] = res[id++];}
    }
  }
};

//Bellman-Ford
template<typename T>
bool find_negative_loop(vector<vector<edge<T>>> &g, int s){
  int n = (int)g.size();
  vector<T> d(n, numeric_limits<T>::max()); d[s] = 0;
  for(int i = 0; i < n; i++) for(int j = 0; j < n; j++) for(edge<T> e : g[j]){
    if(d[j] != numeric_limits<T>::max() && d[e.to] > d[j] + e.cost){
      d[e.to] = d[j] + e.cost;
      if(i == n - 1) return true;
    }
  }
  return false;
}
template<typename T>
bool find_negative_loop_whole(vector<vector<edge<T>>> &g){
  int n = (int)g.size();
  vector<T> d(n);
  for(int i = 0; i < n; i++) for(int j = 0; j < n; j++) for(edge<T> e : g[j]){
    if(d[e.to] > d[j] + e.cost){
      d[e.to] = d[j] + e.cost;
      if(i == n - 1) return true;
    }
  }
  return false;
}
template<typename T>
vector<T> bellman(const vector<vector<edge<T>>> &g, int s){
  int n = (int)g.size();
  vector<T> d(n, numeric_limits<T>::max()); d[s] = 0;
  while(true){
    bool update = false;
    for(int i = 0; i < n; i++) for(edge<T> e : g[i]){
      if(d[i] != numeric_limits<T>::max() && d[e.to] > d[i] + e.cost){
        d[e.to] = d[i] + e.cost;
        update = true;
      }
    }
    if(!update) break;
  }
  return d;
}


//BIGINT (using FFT)
struct BigInt {
  using ll = long long;
  using vll = vector<ll>;
  constexpr static ll base = 1000000000;
  constexpr static ll base_digits = 9;
  vll a;
  ll sign;
  BigInt():sign(1){}
  BigInt(ll v){*this=v;}
  BigInt(const string &s){read(s);}
  void operator=(const BigInt &v){ sign=v.sign; a=v.a;}
  void operator=(ll v){
    sign=1;
    if(v<0) { sign=-1; v=-v;}
    for(;v>0;v=v/base) a.push_back(v%base);
  }
  BigInt operator+(const BigInt &v) const{
    if(sign==v.sign){
      BigInt res=v;
      for(ll i=0,carry=0;i<(ll)max(a.size(),v.a.size())||carry;++i){
        if(i==(ll)res.a.size()) res.a.push_back(0);
        res.a[i]+=carry+(i<(ll)a.size()?a[i]:0);
        carry=res.a[i]>=base;
        if(carry) res.a[i]-=base;
      }
      return res;
    }
    return *this -(-v);
  }
  BigInt operator-(const BigInt &v) const{
    if(sign==v.sign){
      if(abs()>=v.abs()){
        BigInt res=*this;
        for(ll i=0,carry=0;i<(ll)v.a.size()||carry;++i){
          res.a[i]-=carry+(i<(ll)v.a.size()?v.a[i]:0);
          carry=res.a[i]<0;
          if(carry) res.a[i]+=base;
        }
        res.trim();
        return res;
      }
      return -(v-*this);
    }
    return *this+(-v);
  }
  void operator*=(ll v){
    if(v<0) { sign=-sign; v=-v;}
    for(ll i=0,carry=0;i<(ll)a.size()|| carry;++i){
      if(i ==(ll)a.size()) a.push_back(0);
      ll cur=a[i] *(ll)v+carry;
      carry=(ll)(cur/base);
      a[i]=(ll)(cur%base);
    }
    trim();
  }
  BigInt operator*(ll v) const{ BigInt res=*this; res *=v; return res;}
  friend pair<BigInt,BigInt> divmod(const BigInt &a1,const BigInt &b1){
    ll norm=base/(b1.a.back()+1);
    BigInt a=a1.abs()*norm;
    BigInt b=b1.abs()*norm;
    BigInt q,r;
    q.a.resize(a.a.size());
    for(ll i=a.a.size()-1;i>=0;i--){
      r *=base;
      r+=a.a[i];
      ll s1=r.a.size()<=b.a.size() ? 0 : r.a[b.a.size()];
      ll s2=r.a.size()<=b.a.size()-1 ? 0 : r.a[b.a.size()-1];
      ll d=((ll)base*s1+s2)/b.a.back();
      r-=b*d;
      while(r<0) { r+=b; --d;}
      q.a[i]=d;
    }
    q.sign=a1.sign*b1.sign;
    r.sign=a1.sign;
    q.trim();
    r.trim();
    return make_pair(q,r/norm);
  }
  BigInt operator/(const BigInt &v) const{ return divmod(*this,v).first;}
  BigInt operator%(const BigInt &v) const{ return divmod(*this,v).second;}
  void operator/=(ll v){
    if(v<0) { sign=-sign; v=-v;}
    for(ll i=(ll)a.size()-1,rem=0;i>=0;--i){
      ll cur=a[i]+rem *(ll)base;
      a[i]=(ll)(cur/v);
      rem=(ll)(cur%v);
    }
    trim();
  }
  BigInt operator/(ll v) const{ BigInt res=*this; res/=v; return res;}
  ll operator%(ll v) const{
    if(v<0) v=-v;
    ll m=0;
    for(ll i=a.size()-1;i>=0;--i) m=(a[i]+m *(ll)base)%v;
    return m*sign;
  }
  void operator+=(const BigInt &v){ *this=*this+v;}
  void operator-=(const BigInt &v){ *this=*this-v;}
  void operator*=(const BigInt &v){ *this=*this*v;}
  void operator/=(const BigInt &v){ *this=*this/v;}
  bool operator<(const BigInt &v) const{
    if(sign!=v.sign) return sign<v.sign;
    if(a.size()!=v.a.size()) return a.size()*sign<v.a.size()*v.sign;
    for(ll i=a.size()-1;i>=0;i--)
      if(a[i]!=v.a[i]) return a[i]*sign<v.a[i]*sign;
    return false;
  }
  bool operator>(const BigInt &v) const{ return v<*this;}
  bool operator<=(const BigInt &v) const{ return !(v<*this);}
  bool operator>=(const BigInt &v) const{ return !(*this<v);}
  bool operator==(const BigInt &v) const{ return !(*this<v)&&!(v<*this);}
  bool operator!=(const BigInt &v) const{ return *this<v|| v<*this;}
  void trim(){ while(!a.empty()&&!a.back()) a.pop_back(); if(a.empty()) sign=1;}
  bool isZero() const{ return a.empty()||(a.size()==1&&!a[0]);}
  BigInt operator-() const{ BigInt res=*this; res.sign=-sign; return res;}
  BigInt abs() const{ BigInt res=*this; res.sign*=res.sign; return res;}
  ll longValue() const{
    ll res=0;
    for(ll i=a.size()-1;i>=0;i--) res=res*base+a[i];
    return res*sign;
  }
  friend BigInt gcd(const BigInt &a,const BigInt &b){ return b.isZero()?a:gcd(b,a%b);}
  friend BigInt lcm(const BigInt &a,const BigInt &b){ return a/gcd(a,b)*b;}
  void read(const string &s){
    sign=1;
    a.clear();
    ll pos=0;
    while(pos<(ll)s.size()&&(s[pos]=='-'|| s[pos]=='+')){
      if(s[pos]=='-') sign=-sign;
      ++pos;
    }
    for(ll i=s.size()-1;i>=pos;i-=base_digits){
      ll x=0;
      for(ll j=max(pos,i-base_digits+1);j<=i;j++) x=x*10+s[j]-'0';
      a.push_back(x);
    }
    trim();
  }
  friend istream &operator>>(istream &stream,BigInt &v){
    string s;
    stream>>s;
    v.read(s);
    return stream;
  }
  friend ostream &operator<<(ostream &stream,const BigInt &v){
    if(v.sign==-1) stream<<'-';
    stream<<(v.a.empty()?0:v.a.back());
    for(ll i=(ll)v.a.size()-2;i>=0;--i) stream<<setw(base_digits)<<setfill('0')<<v.a[i];
    return stream;
  }
  static vll convert_base(const vll &a,ll old_digits,ll new_digits){
    vll p(max(old_digits,new_digits)+1);
    p[0]=1;
    for(ll i=1;i<(ll)p.size();i++) p[i]=p[i-1]*10;
    vll res;
    ll cur=0;
    ll cur_digits=0;
    for(ll i=0;i<(ll)a.size();i++){
      cur+=a[i]*p[cur_digits];
      cur_digits+=old_digits;
      while(cur_digits>=new_digits){
        res.push_back(signed(cur%p[new_digits]));
        cur/=p[new_digits];
        cur_digits-=new_digits;
      }
    }
    res.push_back((signed)cur);
    while(!res.empty()&&!res.back()) res.pop_back();
    return res;
  }
  static vll karatsubaMultiply(const vll &a,const vll &b){
    ll n=a.size();
    vll res(n+n);
    if(n<=32){
      for(ll i=0;i<n;i++)
        for(ll j=0;j<n;j++)
          res[i+j]+=a[i]*b[j];
      return res;
    }
    ll k=n>>1;
    vll a1(a.begin(),a.begin()+k);
    vll a2(a.begin()+k,a.end());
    vll b1(b.begin(),b.begin()+k);
    vll b2(b.begin()+k,b.end());
    vll a1b1=karatsubaMultiply(a1,b1);
    vll a2b2=karatsubaMultiply(a2,b2);
    for(ll i=0;i<k;i++) a2[i]+=a1[i];
    for(ll i=0;i<k;i++) b2[i]+=b1[i];
    vll r=karatsubaMultiply(a2,b2);
    for(ll i=0;i<(ll)a1b1.size();i++) r[i]-=a1b1[i];
    for(ll i=0;i<(ll)a2b2.size();i++) r[i]-=a2b2[i];
    for(ll i=0;i<(ll)r.size();i++) res[i+k]+=r[i];
    for(ll i=0;i<(ll)a1b1.size();i++) res[i]+=a1b1[i];
    for(ll i=0;i<(ll)a2b2.size();i++) res[i+n]+=a2b2[i];
    return res;
  }
  BigInt operator*(const BigInt &v) const{
    constexpr static ll nbase = 10000;
    constexpr static ll nbase_digits = 4;
    vll a=convert_base(this->a,base_digits,nbase_digits);
    vll b=convert_base(v.a,base_digits,nbase_digits);
    if(a.empty()) a.push_back(0);
    if(b.empty()) b.push_back(0);
    vll c=FFT::multiply(a,b);
    BigInt res;
    res.sign=sign*v.sign;
    for(ll i=0,carry=0;i<(ll)c.size();i++){
      ll cur=c[i]+carry;
      res.a.push_back((ll)(cur%nbase));
      carry=(ll)(cur/nbase);
      if(i+1==(int)c.size()&&carry>0) c.push_back(0);
    }
    res.a=convert_base(res.a,nbase_digits,base_digits);
    res.trim();
    return res;
  }
};


//Binary Search
#define bis(_ac,_wa,_f) [&]{ll ac=_ac,wa=_wa;while(abs(ac-wa)>1){ll wj=(ac+wa)>>1;(_f(wj)?ac:wa)=wj;}return ac;}()


//Binary Indexed Tree (Fenwick Tree)
template<typename T>
struct BIT{
  int _n;
  vector<T> data;
  BIT(): _n(0) {};
  explicit BIT(int n): _n(n), data(n) {};
  void add(int p, T x){
    assert(0 <= p && p < _n);
    p++;
    while(p <= _n){
      data[p-1] += x;
      p += p&-p;
    }
  }
  T sum(int r){
    T s = 0;
    while(r > 0){
      s += data[r-1];
      r -= r&-r;
    }
    return s;
  }
  T sum(int l, int r){
    assert(0 <= l && l <= r && r <= _n);
    return sum(r) - sum(l);
  }
};

//BFS
vector<vector<int>> dist(h, vector<int> (w, INF));
queue<P> q;
auto push = [&](int i, int j, int x){
  if(dist[i][j] != INF) return;
  dist[i][j] = x;
  q.push(P(i, j));
};
auto pop = [&](){
  while(!q.empty()){
    int i = q.front().first;
    int j = q.front().second; q.pop();
    rep(dir, 0, 4){
      int ni = i + di[dir], nj = j + dj[dir];
      if(ni < 0 || ni >= h || nj < 0 || nj >= w) continue;
      if(s[ni][nj] == '#') continue;
    }
  }
};

//BFS/GRID
auto grid_bfs = [&](int si, int sj) -> vector<vector<int>> {
  vector<vector<int>> dist(h, vector<int> (w, INF));
  queue<pair<int, int>> q;
  dist[si][sj] = 0;
  q.push(make_pair(si, sj));
  while (!q.empty()) {
    auto [i, j] = q.front(); q.pop();
    for (int dir = 0; dir < 4; dir++) {
      int ni = i + di[dir], nj = j + dj[dir];
      if (ni < 0 || ni >= h || nj < 0 || nj >= w) continue;
      if (s[ni][nj] == '#') continue;
      if (dist[ni][nj] <= dist[i][j] + 1) continue;
      dist[ni][nj] = dist[i][j] + 1;
      q.push(make_pair(ni, nj));
    }
  }
  return dist;
};

//CHANGE MAX, MIN
template<typename T1,typename T2>bool chmin(T1& x,const T2&y){if(x>y){x=y;return true;}else return false;}
template<typename T1,typename T2>bool chmax(T1& x,const T2&y){if(x<y){x=y;return true;}else return false;}


//Chinese Remainder Theorem
// let (r, m) be return value, the solution is x ≡ r (mod. m)
pair<long long, long long> CRT(const vector<long long>& b, const vector<long long>& m) {
  long long r = 0, M = 1;
  for (int i = 0; i < (int)b.size(); i++) {
    long long p, q;
    long long d = extGCD(M, m[i], p, q); // p is inv of M/d (mod. m[i]/d)
    if ((b[i] - r) % d != 0) return make_pair(0, -1);
    long long tmp = (b[i] - r)/d*p%(m[i]/d);
    r += M*tmp;
    M *= m[i]/d;
  }
  return make_pair((r%M+M)%M, M);
}

// x ≡ r (mod m), r = m1*x1+b1 = m2*x2+b2
pair<long, long> CRT2(long long b1, long long m1, long long b2, long long m2){
  long long p, q;
  long long d = extGCD(m1, m2, p, q); // p is inv of m1/d (mod m2/d)
  if((b2 - b1)%d != 0) return make_pair(0, -1); // no answer
  long long m = m1*(m2/d); // lcm(m1, m2)
  long long tmp = (b2 - b1)/d*p%(m2/d);
  long long r = mod(b1 + m1*tmp, m);
  return make_pair(r, m);
}


//COMBINATION
ll combination(ll n, ll r){
  ll res = 1; r = min(r, n - r);
  rep(i, 0, r){ res *= n - i; res /= i + 1;}
  return res;
}
//COMBINATION_MODver.
mint mod_combination(int n, int r){
  mint res = 1; r = min(r, n - r);
  rep(i, 0, r){ res *= n - i; res /= i + 1;}
  return res;
}
//COMBINATION ENUMERATE
struct modinv{
  int n;
  vector<mint> d;
  modinv(): n(2), d({0, 1}) {}
  mint operator()(int i){ while(n <= i){ d.push_back(-d[mod%n]*(mod/n)); n++;} return d[i];}
  mint operator[](int i) const { return d[i];}
}invs;
struct modfact{
  int n;
  vector<mint> d;
  modfact(): n(2), d({1, 1}) {}
  mint operator()(int i) { while(n <= i){ d.push_back(d.back()*n); n++;} return d[i];}
  mint operator[](int i) const { return d[i];}
}facs;
struct modfactinv{
  int n;
  vector<mint> d;
  modfactinv(): n(2), d({1, 1}) {}
  mint operator()(int i) { while(n <= i){ d.push_back(d.back()*invs(n)); n++;} return d[i];}
  mint operator[](int i) const { return d[i];}
}ifacs;
mint perm(int n, int k) { return facs(n)*ifacs(n-k); }
mint comb(int a, int b){
  if(a < b || b < 0) return 0;
  return facs(a)*ifacs(b)*ifacs(a-b);
}


//COMPRESS
template<typename T>
tuple<unordered_map<T, T>, unordered_map<T, T>> compress(vector<T> a) {
  unordered_map<T, T> mp;  // original_value to compressed_index
  unordered_map<T, T> mpp; // compressed_index to original_value
  sort(a.begin(), a.end());
  a.erase(unique(a.begin(), a.end()), a.end());
  for(int i = 0; i < a.size(); i++) {
    mp[a[i]] = i; mpp[i] = a[i];
  }
  return {mp, mpp};
}

//Coodintate Compression
template<typename T=int>
struct CC {
  bool initialized;
  vector<T> xs;
  CC(): initialized(false) {}
  void add(T x) { xs.push_back(x);}
  void init() {
    sort(xs.begin(), xs.end());
    xs.erase(unique(xs.begin(), xs.end()), xs.end());
    initialized = true;
  }
  int operator()(T x) {
    if (!initialized) init();
    return upper_bound(xs.begin(), xs.end(), x) - xs.begin() - 1;
  }
  T operator[](int i) {
    if (!initialized) init();
    return xs[i];
  }
  int size() {
    if (!initialized) init();
    return xs.size();
  }
};


//Diameter
template<typename T>
void dfs(int p, vector<vector<edge<T>>> &g, vector<T> &dist) {
  if(dist[p] == -1) dist[p] = 0;
  int cur = dist[p];
  for(edge<T> e : g[p]) {
    int to = e.to; T cost = e.cost;
    if(dist[to] != -1) continue;
    dist[to] = cur + cost;
    dfs(to, g, dist);
  }
}
template<typename T>
pair<pair<int, int>, T> diameter(vector<vector<edge<T>>> &g) {
  int n = (int)g.size();
  vector<T> d0(n, -1), d1(n, -1);
  dfs(0, g, d0);
  auto p1 = max_element(d0.begin(), d0.end()) - d0.begin();
  dfs(p1, g, d1);
  auto p2 = max_element(d1.begin(), d1.end()) - d1.begin();
  T d = *max_element(d1.begin(), d1.end());
  return make_pair(make_pair(p1, p2), d);
}


//DIJKSTRA
template<typename T>
vector<long long> dijkstra(const vector<vector<edge<T>>> &g, int s){
  const long long LINF = 1001002003004005006ll;
  int n = (int)g.size();
  vector<long long> dist(n, LINF);
  priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<pair<long long, int>>> q;
  q.emplace(0, s);
  while(!q.empty()){
    int v = q.top().second;
    long long d = q.top().first; q.pop();
    if(dist[v] <= d) continue;
    dist[v] = d;
    for(edge e : g[v]) q.emplace(e.cost+d, e.to);
  }
  return dist;
}


//Dirction
const int di[] = {-1, 0, 1, 0, -1, -1, 1, 1};
const int dj[] = {0, -1, 0, 1, -1, 1, -1, 1};


//DIVISOR
vector<ll> divisor(ll n){
  vector<ll> res;
  for(ll i = 1; i*i <= n; i++){
    if(n%i == 0){
      res.push_back(i);
      if(i != n/i) res.push_back(n/i);
    }
  }
  sort(res.begin(), res.end());
  return res;
}


//Edge
template<typename T=int>
struct edge{
  int from, to, id; T cost;
  edge(int to=0, T cost=1): to(to), cost(cost) {}
  //edge(int to=0, int id=0): to(to), id(id) {}
  //edge(int from=0, int to=0, T cost=1): from(from), to(to), cost(cost) {}
  //edge(int to=0, int from=0, int id=0): to(to), from(from), id(id) {}
  bool operator<(const edge &e) const{return cost < e.cost;}
  bool operator>(const edge &e) const{return cost > e.cost;}
};
template<typename T=int>
using graph = vector<vector<edge<T>>>;


//Eratosthenes
template<typename T>
struct Eratosthenes{
  vector<bool> isprime;
  vector<T> sieves;
  vector<T> minfactor;
  vector<T> mobius;
  Eratosthenes(T n=0):isprime(n+1, true), minfactor(n+1, -1), mobius(n+1, 1){
    isprime[1] = false;
    minfactor[1] = 1;
    for(T i = 2; i <= n; i++){
      if(!isprime[i]) continue;
      minfactor[i] = i;
      mobius[i] = -1;
      for(T j = i*2; j <= n; j += i){
        isprime[j] = false;
        if(minfactor[j] == -1) minfactor[j] = i;
        if((j/i)%i) mobius[j] = -mobius[j];
        else mobius[j] = 0;
      }
    }
    for(T i = 2; i <= n; i++) if(isprime[i]) sieves.emplace_back(i);
  }
  vector<pair<T, T>> factorize(T n){
    vector<pair<T, T>> res;
    while(n > 1){
      int p = minfactor[n];
      int exp = 0;
      while(minfactor[n] == p){
        n /= p;
        exp++;
      }
      res.emplace_back(p, exp);
    }
    return res;
  }
  vector<T> divisors(T n){
    vector<T> res({1});
    auto pf = factorize(n);
    for(auto p : pf){
      int s = (int)res.size();
      for(int i = 0; i < s; i++){
        T v = 1;
        for(int j = 0; j < p.second; j++){
          v *= p.first;
          res.push_back(res[i]*v);
        }
      }
    }
    return res;
  }
};



//Euler_phi
int euler_phi(int n){
  int res = n;
  for(int i = 2; i*i < n; i++){
    if(n%i == 0){
      res = res/i*(i - 1);
      for(; n%i == 0; n /= i);
    }
  }
  if(n != 1) res = res/n*(n - 1);
  return res;
}


//EXTGCD  a*x + b*y = extGCD
ll extGCD(ll a, ll b, ll &x, ll &y){
  if(b == 0) { x = 1; y = 0; return a;}
  ll d = extGCD(b, a%b, y, x);
  y -= a/b*x;
  return d;
}


//FFT
namespace FFT{
  struct num{
    double x,y;
    num(){x=y=0;}
    num(double x,double y):x(x),y(y){}
  };
  inline num operator+(num a,num b){ return num(a.x+b.x,a.y+b.y);}
  inline num operator-(num a,num b){ return num(a.x-b.x,a.y-b.y);}
  inline num operator*(num a,num b){ return num(a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x);}
  inline num conj(num a){ return num(a.x,-a.y);}
  int base=1;
  vector<num> rts={{0,0},{1,0}};
  vector<int> rev={0,1};
  const double PI=acosl(-1.0);
  void ensure_base(int nbase){
    if(nbase<=base) return;
    rev.resize(1<<nbase);
    for(int i=0;i<(1<<nbase);i++) rev[i]=(rev[i>>1]>>1)+((i&1)<<(nbase-1));
    rts.resize(1<<nbase);
    while(base<nbase){
      double angle=2*PI/(1<<(base+1));
      for(int i=1<<(base-1);i<(1<<base);i++){
        rts[i<<1]=rts[i];
        double angle_i=angle*(2*i+1-(1<<base));
        rts[(i<<1)+1]=num(cos(angle_i),sin(angle_i));
      }
      base++;
    }
  }
  void fft(vector<num> &a,int n=-1){
    if(n == -1) n = (int)a.size();
    assert((n&(n-1))==0);
    int zeros=__builtin_ctz(n);
    ensure_base(zeros);
    int shift=base-zeros;
    for(int i=0;i<n;i++) if(i<(rev[i]>>shift)) swap(a[i],a[rev[i]>>shift]);
    for(int k=1;k<n;k<<=1){
      for(int i=0;i<n;i+=2*k){
        for(int j=0;j<k;j++){
          num z=a[i+j+k]*rts[j+k];
          a[i+j+k]=a[i+j]-z;
          a[i+j]=a[i+j]+z;
        }
      }
    }
  }
  vector<num> fa;
  template<typename T>
  vector<long long> multiply(const vector<T> &a,const vector<T> &b){
    int need=(int)a.size()+b.size()-1;
    int nbase=0;
    while((1<<nbase)<need) nbase++;
    ensure_base(nbase);
    int sz=1<<nbase;
    if(sz>(int)fa.size()) fa.resize(sz);
    for(int i=0;i<sz;i++){
      int x=(i<(int)a.size()?a[i]:0);
      int y=(i<(int)b.size()?b[i]:0);
      fa[i]=num(x,y);
    }
    fft(fa,sz);
    num r(0,-0.25/sz);
    for(int i=0;i<=(sz>>1);i++){
      int j=(sz-i)&(sz-1);
      num z=(fa[j]*fa[j]-conj(fa[i]*fa[i]))*r;
      if(i!=j)
        fa[j]=(fa[i]*fa[i]-conj(fa[j]*fa[j]))*r;
      fa[i]=z;
    }
    fft(fa,sz);
    vector<long long> res(need);
    for(int i=0;i<need;i++) res[i]=fa[i].x+0.5;
    return res;
  }
};


//FLOOR_SUM ∑[i=0, n−1](floor(a×i+b)/m)
template<typename T>
T floor_sum(T n, T m, T a, T b){
  T ans = 0;
  if(a >= m){
    ans += (n - 1)*n*(a/m)/2;
    a %= m;
  }
  if(b >= m){
    ans += n*(b/m);
    b %= m;
  }
  T y_max = (a*n + b)/m, x_max = (y_max*m - b);
  if(y_max == 0) return ans;
  ans += (n - (x_max + a - 1)/a)*y_max;
  ans += floor_sum(y_max, a, m, (a - x_max%a)%a);
  return ans;
}


//FPS (Formal Power Series)
struct fps : vector<mint> {
#define d (*this)
#define s int(vector<mint>::size())
  template<class...Args> fps(Args...args): vector<mint>(args...) {}
  fps(initializer_list<mint> a): vector<mint>(a.begin(),a.end()) {}
  void rsz(int n) { if (s < n) resize(n);}
  fps& low_(int n) { resize(n); return d;}
  fps low(int n) const { return fps(d).low_(n);}
  mint& operator[](int i) { rsz(i+1); return vector<mint>::operator[](i);}
  mint operator[](int i) const { return i<s ? vector<mint>::operator[](i) : 0;}
  mint operator()(mint x) const {
    mint r;
    for (int i = s-1; i >= 0; --i) r = r*x+d[i];
    return r;
  }
  fps operator-() const { fps r(d); rep(i, 0, s) r[i] = -r[i]; return r;}
  fps& operator+=(const fps& a) { rsz((int)a.size()); rep(i,0,a.size()) d[i] += a[i]; return d;}
  fps& operator-=(const fps& a) { rsz((int)a.size()); rep(i,0,a.size()) d[i] -= a[i]; return d;}
  fps& operator*=(const fps& a) { return d = ntt(d, a);}
  fps& operator*=(mint a) { rep(i, 0, s) d[i] *= a; return d;}
  fps& operator/=(mint a) { rep(i, 0, s) d[i] /= a; return d;}
  fps operator+(const fps& a) const { return fps(d) += a;}
  fps operator-(const fps& a) const { return fps(d) -= a;}
  fps operator*(const fps& a) const { return fps(d) *= a;}
  fps operator*(mint a) const { return fps(d) *= a;}
  fps operator/(mint a) const { return fps(d) /= a;}
  fps operator~() const {
    fps r({d[0].inv()});
    for (int i = 1; i < s; i <<= 1) r = r*mint(2) - (r*r*low(i<<1)).low(i<<1);
    return r.low_(s);
  }
  fps& operator/=(const fps& a) { int w = s; d *= ~a; return d.low_(w);}
  fps operator/(const fps& a) const { return fps(d) /= a;}
  fps integ() const {
    fps r;
    rep(i, 0, s) r[i+1] = d[i]/(i+1);
    return r;
  }
#undef s
#undef d
};


//FRACTION
template<typename T=ll>
struct frac{
  T a, b;
  frac(T _a=0, T _b=1): a(_a), b(_b) {
    if(b == 0) { a = 1; return;}
    if(b < 0) a = -a, b = -b;
    T g = gcd(abs(a), b);
    a /= g; b /= g;
  }
  bool operator<(const frac& x) const{ return a*x.b < x.a*b;}
  bool operator==(const frac& x) const{ return a == x.a && b == x.b;}
};


//Fraction
struct frac{
  __int128_t num, den;
  void simplify(){
    __int128_t d = gcd(num, den);
    num /= (den >= 0 ? d : -d);
    den /= (den >= 0 ? d : -d);
  }
  frac(): num(0), den(1) {}
  frac(const frac &a) { num = a.num; den = a.den; }
  frac(__int128_t n): num(n), den(1) {}
  // false を添えると約分されない
  frac(__int128_t numerator, __int128_t denominator, bool with_simplify = true) : num(numerator), den(denominator) { if(with_simplify) simplify(); }
  friend const frac operator+(const frac &a){ return a; }
  friend const frac operator-(const frac &a){ return {-a.num, a.den, false}; }
  friend const frac operator+(const frac &a, const frac &b){
    return {a.num * b.den + b.num * a.den, a.den * b.den};
  }
  friend const frac operator-(const frac &a, const frac &b){
    return {a.num * b.den - b.num * a.den, a.den * b.den};
  }
  friend const frac operator*(const frac &a, const frac &b){
    __int128_t gcd_tmp1 = gcd(a.num, b.den), gcd_tmp2 = gcd(b.num, a.den);
    return {(a.num / gcd_tmp1) * (b.num / gcd_tmp2), (a.den / gcd_tmp2) * (b.den / gcd_tmp1), false};
  }
  friend const frac operator/(const frac &a, const frac &b){
    __int128_t gcd_tmp1 = gcd(a.num, b.num), gcd_tmp2 = gcd(b.den, a.den);
    return {(b.num < 0 ? -1 : 1) * (a.num / gcd_tmp1) * (b.den / gcd_tmp2), (b.num < 0 ? -1 : 1) * (a.den / gcd_tmp2) * (b.num / gcd_tmp1), false};
  }
  friend bool operator==(const frac &a, const frac &b){ return a.num == b.num && a.den == b.den; }
  friend bool operator!=(const frac &a, const frac &b){ return a.num != b.num || a.den != b.den; }
  friend bool operator>(const frac &a, const frac &b) {
    return a.num * b.den > b.num * a.den;
  }
  friend bool operator>=(const frac &a, const frac &b) {
    return a.num * b.den >= b.num * a.den;
  }
  friend bool operator<(const frac &a, const frac &b) {
    return a.num * b.den < b.num * a.den;
  }
  friend bool operator<=(const frac &a, const frac &b){
    return a.num * b.den <= b.num * a.den;
  }
  frac &operator=(const frac &a){
    num = a.num;
    den = a.den;
    return *this;
  }
  frac &operator+=(const frac &a){
    num = num * a.den + a.num * den;
    den *= a.den;
    simplify();
    return *this;
  }
  frac &operator-=(const frac &a){
    num = num * a.den - a.num * den;
    den *= a.den;
    simplify();
    return *this;
  }
  frac &operator*=(const frac &a){
    __int128_t gcd_tmp1 = gcd(num, a.den), gcd_tmp2 = gcd(a.num, den);
    num = (num / gcd_tmp1) * (a.num / gcd_tmp2);
    den = (den / gcd_tmp2) * (a.den / gcd_tmp1);
    return *this;
  }
  frac &operator/=(const frac &a){
    __int128_t gcd_tmp1 = gcd(num, a.num), gcd_tmp2 = gcd(a.den, den);
    num = (a.num < 0 ? -1 : 1) * (num / gcd_tmp1) * (a.den / gcd_tmp2);
    den = (a.num < 0 ? -1 : 1) * (den / gcd_tmp2) * (a.num / gcd_tmp1);
    return *this;
  }
  long long numerator() const { return num; }
  long long denomitnator() const { return den; }
  long long floor() const { return num / den; }
  long long ceil() const { return (num + den - 1) / den; }
  double real() const { return (double)num / den; }
  frac abs() const { return {(num >= 0 ? num : -num), den, false}; }
  frac inverse() const {
    assert(num != 0);
    return {(num >= 0 ? den : -den), (num >= 0 ? num : -num), false};
  }
};
istream& operator>>(istream &is, frac &a) {
  string buf;
  is >> buf;
  a.num = a.den = 0;
  int i = (buf[0] == '-');
  for(; i < buf.size() && buf[i] != '/'; i++) a.num = a.num * 10 + buf[i] - '0';
  if(i == buf.size()) a.den = 1;
  else for(i = i + 1; i < buf.size(); i++) a.den = a.den * 10 + buf[i] - '0';
  if(buf[0] == '-') a.num *= -1;
  a.simplify();
  return is;
}
ostream& operator<<(ostream &os, const frac &a) {
  if(a.den == 0) os << (a.num >= 0 ? "inf" : "-inf");
  else if(a.den == 1) os << (long long)a.num;
  else os << (long long)a.num << '/' << (long long)a.den;
  return os;
}


//GCD,LCM
ll gcd(ll a, ll b) { return b ? gcd(b, a%b) : a;}
ll lcm(ll a, ll b) { return a/gcd(a, b)*b;}


//GraphLink
template<typename T>
struct graphLink{
  vector<int> ord, low, parent, cmp;
  vector<vector<edge<T>>> g, h;
  vector<pair<int, int>> bridges;
  int cnt, v;
  graphLink(vector<vector<edge<T>>> &s, int root=0){
    int n = (int)s.size();
    ord.resize(n, -1); low.resize(n, 0); parent.resize(n, -1); cmp.resize(n, -1);
    cnt = 0; v = n; g = s;
    dfs(root);
  }
  bool is_bridge(int x, int y){
    if(ord[x] > ord[y]) swap(x, y);
    return ord[x] < low[y];
  }
  void dfs(int cur, int prev=-1){
    low[cur] = cnt;
    ord[cur] = cnt++;
    for(auto x : g[cur]){
      if(x.to == prev) continue;
      if(ord[x.to] < 0){
        parent[x.to] = cur;
        dfs(x.to, cur);
        low[cur] = min(low[cur], low[x.to]);
      }
      else{
        low[cur] = min(low[cur], ord[x.to]);
      }
      if(is_bridge(cur, x.to)){
        int a = min(cur, x.to);
        int b = max(cur, x.to);
        bridges.emplace_back(make_pair(a, b));
      }
    }
  }
  set<int> artPoint(int root=0){
    set<int> st;
    int num = 0;
    for(int i = 0; i < v; i++){
      if(parent[i] < 0) continue;
      if(parent[i] == root) num++;
      else if(ord[parent[i]] <= low[i]) st.insert(parent[i]);
    }
    if(num >= 2) st.insert(0);
    return st;
  }
};


//HL-decomposition (using SegmentTree, edge)
using HL_S = long long;
HL_S hl_op(HL_S a, HL_S b){ return a + b;}
HL_S hl_e(){ return 0;}
template<typename T=long long>
struct HL {
  int n, root;
  vector<vector<int>> to;
  SEG<HL_S, hl_op, hl_e> t;
  HL(int n=0): n(n), to(n) {}
  HL(const vector<vector<edge<T>>>& s): n((int)s.size()), to(n) {
    for(int i = 0; i < n; i++) for(auto j : s[i]) to[i].push_back(j.to);
    init();
  }
  void addEdge(int a, int b) {
    to[a].push_back(b);
    to[b].push_back(a);
  }
  vector<int> tsz, pa, dep;
  int tfs(int v, int p=-1) {
    pa[v] = p; tsz[v] = 1;
    for (int i = 0; i < (int)to[v].size(); i++) {
      int& u = to[v][i];
      if (u == p) swap(u, to[v].back());
      if (u == p) break;
      dep[u] = dep[v] + 1;
      tsz[v] += tfs(u, v);
      if (tsz[u] > tsz[to[v][0]]) swap(u, to[v][0]);
    }
    if (~p) to[v].pop_back();
    return tsz[v];
  }
  vector<int> in, out, nxt, kv;
  void dfs(int v) {
    in[v] = (int)size(kv); kv.push_back(v);
    for (int i = 0; i < (int)to[v].size(); i++) {
      int u = to[v][i];
      nxt[u] = i ? u : nxt[v];
      dfs(u);
    }
    out[v] = (int)size(kv);
  }
  void init(int v=0) {
    tsz = pa = dep = vector<int>(n);
    tfs(root = v);
    in = out = nxt = vector<int>(n);
    nxt[root] = root;
    dfs(root);
    t = SEG<HL_S, hl_op, hl_e>(n);
  }
  int up(int v, int l) {
    while (~v) {
      int u = nxt[v], nl = dep[v] - dep[u] + 1;
      if (l < nl) return kv[in[v]-l];
      l -= nl; v = pa[u];
    }
    return -1;
  }
  int lca(int a, int b) {
    while (nxt[a] != nxt[b]) {
      if (in[a] < in[b]) swap(a, b);
      a = pa[nxt[a]];
    }
    if (in[a] < in[b]) swap(a, b);
    return b;
  }
  int len(int a, int b) { return dep[a]+dep[b]-dep[lca(a, b)]*2;}
  void initSeg(vector<HL_S>& a) { t = SEG<HL_S, hl_op, hl_e>(a);}
  void add(int v, HL_S x) { t.set(in[v], x);}
  HL_S sum(int p, int v) {
    HL_S res = 0;
    int ip = ~p ? in[p] : -1;
    while (~v) {
      int u = nxt[v];
      if (in[u] <= ip) {
        res += t.prod(ip+1, in[v]+1);
        break;
      }
      res += t.prod(in[u], in[v]+1);
      v = pa[u];
    }
    return res;
  }
  HL_S sums(int a, int b) { int c = lca(a, b); return sum(c, a)+sum(c, b);}
  HL_S sumr(int v) { return t.prod(in[v], out[v]);}
};


//ImplicitTreap
// モノイドの型
using S = int;
// S × S → S を計算する関数
S op(S a, S b) { return a + b; }
// モノイドの単位元
S e() { return 0; }
// 作用素 (!=の定義を必ずする！！)
using F = int;
// F × S(区間sizeをszにいれる)
S mapping(F f, S b, int sz) {
  return f + b;
}
// F × F
F composition(F a, F b) { return a + b; }
// 作用素の単位元
F id() { return 0; }

template <class S, S (*op)(S, S), S (*e)(), class F, S (*mapping)(F, S, int), F (*composition)(F, F), F (*id)()>
struct ImplicitTreap {
  struct xorshift {
    using u32 = uint32_t;
    u32 x = 123456789, y = 362436069, z = 521288629, w = 88675123;
    xorshift(u32 seed = 0) { z ^= seed; }
    u32 operator()() {
      u32 t = x ^ (x << 11);
      x = y; y = z; z = w;
      return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
    }
  };
  xorshift rnd;
  struct Node {
    S val, acc;
    F lazy;
    int pri, cnt;
    bool rev;
    Node *l, *r;
    Node(S val_, int pri_) : val(val_), acc(e()), lazy(id()), pri(pri_), cnt(1), rev(false), l(nullptr), r(nullptr) {}
  }* root = nullptr;
  using Tree = Node*;
  int count(Tree t) { return t ? t->cnt : 0;}
  S acc(Tree t) { return t ? t->acc : e();}
  void proc(Tree& t) {
    if(t) {
      t->cnt = 1 + count(t->l) + count(t->r);
      t->acc = op(acc(t->l), op(t->val, acc(t->r)));
    }
  }
  void eval(Tree& t) {
    if(t and t->rev) {
      t->rev = false;
      std::swap(t->l, t->r);
      if(t->l) t->l->rev ^= 1;
      if(t->r) t->r->rev ^= 1;
    }
    if(t and t->lazy != id()) {
      if(t->l) {
        t->l->lazy = composition(t->lazy, t->l->lazy);
        t->l->acc = mapping(t->lazy, t->l->acc, count(t->l));
      }
      if(t->r) {
        t->r->lazy = composition(t->lazy, t->r->lazy);
        t->r->acc = mapping(t->lazy, t->r->acc, count(t->r));
      }
      t->val = mapping(t->lazy, t->val, 1);
      t->lazy = id();
    }
    proc(t);
  }
  void split(Tree t, int key, Tree& l, Tree& r) {
    if(!t) {
      l = r = nullptr;
      return;
    }
    eval(t);
    int implicit_key = count(t->l) + 1;
    if(key < implicit_key) {
      split(t->l, key, l, t->l), r = t;
    } else {
      split(t->r, key - implicit_key, t->r, r), l = t;
    }
    proc(t);
  }
  void merge(Tree& t, Tree l, Tree r) {
    eval(l);
    eval(r);
    if(!l or !r)
      t = l ? l : r;
    else if(l->pri > r->pri)
      merge(l->r, l->r, r), t = l;
    else
      merge(r->l, l, r->l), t = r;
    proc(t);
  }
  void insert(Tree& t, int key, Tree item) {
    Tree t1, t2;
    split(t, key, t1, t2);
    merge(t1, t1, item);
    merge(t, t1, t2);
  }
  void erase(Tree& t, int key) {
    Tree t1, t2, t3;
    split(t, key + 1, t1, t2);
    split(t1, key, t1, t3);
    merge(t, t1, t2);
  }
  void update(Tree t, int l, int r, F x) {
    Tree t1, t2, t3;
    split(t, l, t1, t2);
    split(t2, r - l, t2, t3);
    t2->lazy = composition(x, t2->lazy);
    t2->acc = mapping(x, t2->acc, count(t2));
    merge(t2, t2, t3);
    merge(t, t1, t2);
  }
  S query(Tree t, int l, int r) {
    Tree t1, t2, t3;
    split(t, l, t1, t2);
    split(t2, r - l, t2, t3);
    S ret = t2->acc;
    merge(t2, t2, t3);
    merge(t, t1, t2);
    return ret;
  }
  int find(Tree t, S x, int offset, bool left = true) {
    if(op(t->acc, x) == x) {
      return -1;
    } else {
      if(left) {
        if(t->l and op(t->l->acc, x) != x)
          return find(t->l, x, offset, left);
        else
          return (op(t->val, x) != x) ? offset + count(t->l) : find(t->r, x, offset + count(t->l) + 1, left);
        
      } else {
        if(t->r and op(t->r->acc, x) != x)
          return find(t->r, x, offset + count(t->l) + 1, left);
        else
          return (op(t->val, x) != x) ? offset + count(t->l) : find(t->l, x, offset, left);
      }
    }
  }
  void reverse(Tree t, int l, int r) {
    if(l > r) return;
    Tree t1, t2, t3;
    split(t, l, t1, t2);
    split(t2, r - l, t2, t3);
    t2->rev ^= 1;
    merge(t2, t2, t3);
    merge(t, t1, t2);
  }
  // [l, r)の先頭がmになるようにシフトさせる。std::rotateと同じ仕様
  void rotate(Tree t, int l, int m, int r) {
    reverse(t, l, r);
    reverse(t, l, l + r - m);
    reverse(t, l + r - m, r);
  }
  
public:
  ImplicitTreap() {}
  ImplicitTreap(std::vector<S> as) {
    std::reverse(as.begin(), as.end());
    for(S& a : as) { insert(0, a); }
  }
  int size() { return count(root); }
  void insert(int pos, S x) { insert(root, pos, new Node(x, rnd()));}
  void update(int l, int r, F x) { update(root, l, r, x);}
  S query(int l, int r) { return query(root, l, r);}
  // 二分探索。[l, r)内のkでMonoid::op(tr[k], x) != xとなる最左/最右のもの。存在しない場合は-1
  // たとえばMinMonoidの場合、x未満の最左/最右の要素の位置を返す
  int find(int l, int r, S x, bool left = true) {
    Tree t1, t2, t3;
    split(root, l, t1, t2);
    split(t2, r - l, t2, t3);
    int ret = find(t2, x, l, left);
    merge(t2, t2, t3);
    merge(root, t1, t2);
    return ret;
  }
  void erase(int pos) { erase(root, pos);}
  void reverse(int l, int r) { reverse(root, l, r);}
  void rotate(int l, int m, int r) { rotate(root, l, m, r);}
  S operator[](int pos) {
    Tree t1, t2, t3;
    split(root, pos + 1, t1, t2);
    split(t1, pos, t1, t3);
    S ret = t3->acc;
    merge(t1, t1, t3);
    merge(root, t1, t2);
    return ret;
  }
};


//IS_PRIME
bool is_prime(long long n) {
  if (n <= 1) return false;
  if (n == 2 || n == 3) return true;
  if (n%2 == 0) return false;
  auto pow_mod = [&](__int128_t a, __int128_t n, __int128_t m) -> __int128_t {
    __int128_t res = 1%m;
    a %= m;
    for (; n > 0; n >>= 1){
      if (n&1) res = (res*a)%m;
      a = (a*a)%m;
    }
    return res;
  };
  vector<long long> a = {2, 325, 9375, 28178, 450775,
    9780504, 1795265022};
  long long s = 0, d = n - 1;
  for (; d%2 == 0; d >>= 1) s++;
  for (auto i : a) {
    if (i%n == 0) return true;
    long long t, x = pow_mod(i, d, n);
    if (x != 1) {
      for (t = 0; t < s; t++) {
        if (x == n - 1) break;
        x = __int128_t(x)*x%n;
      }
      if (t == s) return false;
    }
  }
  return true;
}


//__int128
std::ostream &operator<<(std::ostream &dest, __int128_t value) {
  std::ostream::sentry s(dest);
  if (s) {
    __uint128_t tmp = value < 0 ? -value : value;
    char buffer[128];
    char *d = std::end(buffer);
    do {
      --d;
      *d = "0123456789"[tmp % 10];
      tmp /= 10;
    } while (tmp != 0);
    if (value < 0) {
      --d;
      *d = '-';
    }
    int len = std::end(buffer) - d;
    if (dest.rdbuf()->sputn(d, len) != len) {
      dest.setstate(std::ios_base::badbit);
    }
  }
  return dest;
}

__int128 parse(string &s) {
  __int128 ret = 0;
  for (int i = 0; i < s.length(); i++)
    if ('0' <= s[i] && s[i] <= '9')
      ret = 10 * ret + s[i] - '0';
  return ret;
}


//KMP (Morris-Pratt)
template<typename T>
struct MP {
  int n;
  T t;
  vector<int> a;
  MP() {}
  MP(const T& t): t(t) {
    n = t.size();
    a = vector<int>(n+1);
    a[0] = -1;
    int j = -1;
    for (int i = 0; i < n; ++i) {
      while (j != -1 && t[j] != t[i]) j = a[j];
      j++;
      a[i+1] = j;
    }
  }
  int operator[](int i) { return a[i];}
  vector<int> findAll(const T& s) {
    vector<int> res;
    int j = 0;
    for (int i = 0; i < s.size(); ++i) {
      while (j != -1 && t[j] != s[i]) j = a[j];
      j++;
      if (j == n) {
        res.push_back(i-j+1);
        j = a[j];
      }
    }
    return res;
  }
};


//Kruskal
template<typename T>
T kruskal(vector<vector<edge<T>>>& g){
  int n = (int)g.size();
  vector<edge<T>> es;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < g[i].size(); j++){
      es.push_back(g[i][j]);
    }
  }
  sort(es.begin(), es.end());
  UnionFind<int> uf(n);
  T res = 0;
  for(int i = 0; i < es.size(); i++){
    edge<T> e = es[i];
    if(!uf.same(e.from, e.to)){
      uf.unite(e.from, e.to);
      res += e.cost;
    }
  }
  return res;
}


//Matrix
template<typename T>
struct Matrix {
  int h, w;
  vector<vector<T>> d;
  Matrix(int h, int w, T val=0): h(h), w(w), d(h, vector<T>(w, val)) {}
  Matrix(vector<vector<T>> a): h((int)a.size()), w((int)a[0].size()), d(a) {}
  const vector<T>& operator[](int i) const { return d[i];}
  vector<T>& operator[](int i) { return d[i];}
  Matrix& operator+=(const Matrix& a) {
    assert(h == a.h && w == a.w);
    for(int i=0;i<h;i++)for(int j=0;j<w;j++) d[i][j] += a[i][j];
    return *this;
  }
  Matrix& operator-=(const Matrix& a) {
    assert(h == a.h && w == a.w);
    for(int i=0;i<h;i++)for(int j=0;j<w;j++) d[i][j] -= a[i][j];
    return *this;
  }
  Matrix& operator*=(const Matrix& a) {
    assert(w == a.h);
    Matrix r(h, a.w);
    for(int k=0;k<w;k++)for(int i=0;i<h;i++)for(int j=0;j<a.w;j++) r[i][j] += d[i][k]*a[k][j];
    w = a.w; for(int i=0;i<h;i++) d[i].resize(w);
    for(int i=0;i<h;i++)for(int j=0;j<w;j++) d[i][j] = r[i][j];
    return *this;
  }
  Matrix operator+(const Matrix& a) const { return Matrix(*this) += a;}
  Matrix operator-(const Matrix& a) const { return Matrix(*this) -= a;}
  Matrix operator*(const Matrix& a) const { return Matrix(*this) *= a;}
  bool operator==(const Matrix& a) {
    assert(h == a.h && w == a.w);
    for(int i=0;i<h;i++)for(int j=0;j<w;j++) if (d[i][j] != a[i][j]) return false;
    return true;
  }
  Matrix& operator+=(const T& a) { for(int i=0;i<h;i++)for(int j=0;j<w;j++) d[i][j] += a; return *this;}
  Matrix& operator-=(const T& a) { for(int i=0;i<h;i++)for(int j=0;j<w;j++) d[i][j] -= a; return *this;}
  Matrix& operator*=(const T& a) { for(int i=0;i<h;i++)for(int j=0;j<w;j++) d[i][j] *= a; return *this;}
  Matrix& operator/=(const T& a) { for(int i=0;i<h;i++)for(int j=0;j<w;j++) d[i][j] /= a; return *this;}
  Matrix& operator%=(const T& a) { for(int i=0;i<h;i++)for(int j=0;j<w;j++) d[i][j] %= a; return *this;}
  Matrix operator+(const T& a) const { return Matrix(*this) += a;}
  Matrix operator-(const T& a) const { return Matrix(*this) -= a;}
  Matrix operator*(const T& a) const { return Matrix(*this) *= a;}
  Matrix operator/(const T& a) const { return Matrix(*this) /= a;}
  Matrix operator%(const T& a) const { return Matrix(*this) %= a;}
  Matrix& unit() {
    assert(h == w);
    for (int i = 0; i < h; i++) d[i][i] = 1;
    return *this;
  }
  Matrix pow(long long t) const {
    assert(h == w);
    if (!t) return Matrix(h, h).unit();
    if (t == 1) return *this;
    Matrix r = pow(t>>1);
    r = r*r;
    if (t&1) r = r*(*this);
    return r;
  }
  Matrix& rot(int deg) {
    assert(1 <= deg && deg <= 3);
    Matrix r(h, w);
    if (deg&1) {
      if (deg == 1) for(int i=0;i<h;i++)for(int j=0;j<w;j++) r[j][h-i-1] = d[i][j];
      if (deg == 3) for(int i=0;i<h;i++)for(int j=0;j<w;j++) r[w-j-1][i] = d[i][j];
      swap(h, w); d.resize(h); for(int i=0;i<h;i++) d[i].resize(w);
    }
    else for(int i=0;i<h;i++)for(int j=0;j<w;j++) r[h-i-1][w-j-1] = d[i][j];
    for(int i=0;i<h;i++)for(int j=0;j<w;j++) d[i][j] = r[i][j];
    return *this;
  }
  /* mint only
  mint det() {
    assert(h == w);
    mint res = 1;
    for (int k = 0; k < h; k++) {
      for (int i = k; i < h; i++) {
        if (d[i][k] == 0) continue;
        if (i != k) {
          swap(d[i], d[k]);
          res = -res;
        }
      }
      if (d[k][k] == 0) return 0;
      res *= d[k][k];
      mint inv = mint(1)/d[k][k];
      for (int j = 0; j < h; j++) d[k][j] *= inv;
      for (int i = k+1; i < h; i++) {
        mint c = d[i][k];
        for (int j = k; j < h; j++) d[i][j] -= d[k][j]*c;
      }
    }
    return res;
  }
   */
  void show() { for(int i=0;i<h;i++)for(int j=0;j<w;j++)cout<<d[i][j]<<" \n"[j+1 == w];}
};


//Maximum Flow
template<class Cap>
struct mf_graph {
  int _n;
  struct _edge {
    int to, rev;
    Cap cap;
  };
  vector<pair<int, int>> pos;
  vector<vector<_edge>> g;
  mf_graph(): _n(0) {}
  mf_graph(int n): _n(n), g(n) {};
  int add_edge(int from, int to, Cap cap) {
    assert(0 <= from && from < _n);
    assert(0 <= to && to < _n);
    assert(0 <= cap);
    int m = int(pos.size());
    pos.push_back({from, int(g[from].size())});
    int from_id = int(g[from].size());
    int to_id = int(g[to].size());
    if (from == to) to_id++;
    g[from].push_back(_edge{to, to_id, cap});
    g[to].push_back(_edge{from, from_id, 0});
    return m;
  }
  struct edge {
    int from, to;
    Cap cap, flow;
  };
  edge get_edge(int i) {
    int m = int(pos.size());
    assert(0 <= i && i < m);
    auto _e = g[pos[i].first][pos[i].second];
    auto _re = g[_e.to][_e.rev];
    return edge{pos[i].first, _e.to, _e.cap + _re.cap, _re.cap};
  }
  vector<edge> edges() {
    int m = int(pos.size());
    vector<edge> result;
    for (int i = 0; i < m; i++) result.push_back(get_edge(i));
    return result;
  }
  void change_edge(int i, Cap new_cap, Cap new_flow) {
    int m = int(pos.size());
    assert(0 <= i && i < m);
    assert(0 <= new_flow && new_flow <= new_cap);
    auto& _e = g[pos[i].first][pos[i].second];
    auto& _re = g[_e.to][_e.rev];
    _e.cap = new_cap - new_flow;
    _re.cap = new_flow;
  }
  
  Cap flow(int s, int t) { return flow(s, t, std::numeric_limits<Cap>::max());}
  Cap flow(int s, int t, Cap flow_limit) {
    assert(0 <= s && s < _n);
    assert(0 <= t && t < _n);
    assert(s != t);
    vector<int> level(_n), iter(_n);
    queue<int> que;
    auto bfs = [&]() {
      fill(level.begin(), level.end(), -1);
      level[s] = 0;
      queue<int>().swap(que);
      que.push(s);
      while (!que.empty()) {
        int v = que.front();
        que.pop();
        for (auto e : g[v]) {
          if (e.cap == 0 || level[e.to] >= 0) continue;
          level[e.to] = level[v] + 1;
          if (e.to == t) return;
          que.push(e.to);
        }
      }
    };
    auto dfs = [&](auto self, int v, Cap up) {
      if (v == s) return up;
      Cap res = 0;
      int level_v = level[v];
      for (int& i = iter[v]; i < int(g[v].size()); i++) {
        _edge& e = g[v][i];
        if (level_v <= level[e.to] || g[e.to][e.rev].cap == 0) continue;
        Cap d = self(self, e.to, std::min(up - res, g[e.to][e.rev].cap));
        if (d <= 0) continue;
        g[v][i].cap += d;
        g[e.to][e.rev].cap -= d;
        res += d;
        if (res == up) break;
      }
      return res;
    };
    Cap flow = 0;
    while (flow < flow_limit) {
      bfs();
      if (level[t] == -1) break;
      fill(iter.begin(), iter.end(), 0);
      while (flow < flow_limit) {
        Cap f = dfs(dfs, t, flow_limit - flow);
        if (!f) break;
        flow += f;
      }
    }
    return flow;
  }
  vector<bool> min_cut(int s) {
    vector<bool> visited(_n);
    queue<int> que;
    que.push(s);
    while (!que.empty()) {
      int p = que.front(); que.pop();
      visited[p] = true;
      for (auto e : g[p]) {
        if (e.cap && !visited[e.to]) {
          visited[e.to] = true;
          que.push(e.to);
        }
      }
    }
    return visited;
  }
};


//MEX
template<typename T=int>
struct MEX {
  set<pair<int, int>> st;
  MEX() { init(); }
  MEX(vector<T>& a) { init(); for (T i : a) insert(i);}
  void init() { st.emplace(INT_MIN, INT_MIN); st.emplace(INT_MAX, INT_MAX); }
  bool contains(int x) const {
    auto [l, r] = *prev(st.lower_bound(make_pair(x+1, x+1)));
    return l <= x && x <= r;
  }
  bool insert(int x) {
    auto nit = st.lower_bound(make_pair(x+1, x+1));
    auto it = prev(nit);
    auto [l, r] = *it;
    auto [nl, nr] = *nit;
    if (l <= x && x <= r) return false;
    if (r == x-1) {
      if (nl == x+1) { st.erase(it); st.erase(nit); st.emplace(l, nr); }
      else { st.erase(it); st.emplace(l, x); }
    } else {
      if (nl == x+1) { st.erase(nit); st.emplace(x, nr); }
      else { st.emplace(x, x); }
    }
    return true;
  }
  int mex(int x=0) const {
    auto [l, r] = *prev(st.lower_bound(make_pair(x+1, x+1)));
    return (l <= x && x <= r) ? r+1 : x;
  }
};


//MINT
struct mint {
  unsigned x;
  mint(): x(0) {}
  mint(ll x):x((x%mod+mod)%mod) {}
  mint operator-() const { return mint(0) - *this;}
  mint operator~() const { return mint(1) / *this;}
  mint& operator+=(const mint& a) { if((x+=a.x)>=mod) x-=mod; return *this;}
  mint& operator-=(const mint& a) { if((x+=mod-a.x)>=mod) x-=mod; return *this;}
  mint& operator*=(const mint& a) { x=(unsigned long long)x*a.x%mod; return *this;}
  mint& operator/=(const mint& a) { x=(unsigned long long)x*a.pow(mod-2).x%mod; return *this;}
  mint operator+(const mint& a) const { return mint(*this) += a;}
  mint operator-(const mint& a) const { return mint(*this) -= a;}
  mint operator*(const mint& a) const { return mint(*this) *= a;}
  mint operator/(const mint& a) const { return mint(*this) /= a;}
  mint pow(ll t) const {
    if (!t) return 1;
    mint res = pow(t>>1);
    res *= res;
    return (t&1)?res*x:res;
  }
  bool operator<(const mint& a) const { return x < a.x;}
  bool operator==(const mint& a) const { return x == a.x;}
  bool operator!=(const mint& a) const { return x != a.x;}
  mint inv() const { return pow(mod-2);}
};
mint ex(mint x, ll t) { return x.pow(t);}
istream& operator>>(istream& i, mint& a) { unsigned long long t; i>>t; a=mint(t); return i;}
ostream& operator<<(ostream& o, const mint& a) { return o<<a.x;}


//MOD
inline ll mod(ll a, ll m) { return (a%m+m)%m;}


//Mo (ABC293-G)
struct D {
  int n; vector<int>& a;
  vector<int> cnt; ll d;
  D(vector<int>& a): n((int)a.size()), a(a), cnt(200005), d(0) {}
  void push_back(int i) { add(i);}
  void push_front(int i) { add(i);}
  void pop_back(int i) { del(i);}
  void pop_front(int i) { del(i);}
  void add(int i, int x=1) {
    auto c3 = [&](ll n) { return n*(n-1)*(n-2)/6;};
    d -= c3(cnt[a[i]]);
    cnt[a[i]] += x;
    d += c3(cnt[a[i]]);
  }
  void del(int i) { add(i,-1);}
  ll get() { return d;}
};

template<typename T=long long>
vector<T> Mo(vector<pair<int, int>> lr, D& d) {
  int q = (int)lr.size(); vector<T> res(q);
  int W = d.n/(sqrt(q)+1)+1;
  vector<int> is(q); iota(is.begin(), is.end(), 0);
  vector<int> li(q); for (int i = 0; i < q; i++) li[i] = lr[i].first/W;
  sort(is.begin(), is.end(), [&](int i, int j) {
    if (li[i] != li[j]) return li[i] < li[j];
    return (li[i]&1) ? (lr[i].second > lr[j].second) : lr[i].second < lr[j].second;
  });
  int l = 0, r = 0;
  for (int i : is) {
    auto [nl, nr] = lr[i];
    while (r < nr) d.push_back(r++);
    while (l > nl) d.push_front(--l);
    while (l < nl) d.pop_front(l++);
    while (r > nr) d.pop_back(--r);
    res[i] = d.get();
  }
  return res;
}


//NTT
struct NTT {
  int ceil_pow2(int n) { int x = 0; while ((1U<<x) < (uint)(n)) x++; return x;}
  int bsf(uint n) { return __builtin_ctz(n);}
  const static int _g = 3; // primitive root !!!
  void butterfly(vector<mint>& a) {
    static mint g = _g; int n = (int)a.size(), h = ceil_pow2(n);
    static bool first = true; static mint sum_e[30];
    if (first) {
      first = false; mint es[30], ies[30]; int cnt2 = bsf(mod - 1);
      mint e = g.pow((mod - 1) >> cnt2), ie = ~e;
      for (int i = cnt2; i >= 2; i--) {
        es[i - 2] = e; ies[i - 2] = ie; e *= e; ie *= ie;
      }
      mint now = 1;
      rep(i, 0, cnt2-1) { sum_e[i] = es[i] * now; now *= ies[i];}
    }
    for (int ph = 1; ph <= h; ph++) {
      int w = 1<<(ph-1), p = 1<<(h-ph); mint now = 1;
      rep(s, 0, w) {
        int offset = s << (h - ph + 1);
        for (int i = 0; i < p; i++) {
          auto l = a[i + offset], r = a[i + offset + p] * now;
          a[i + offset] = l + r; a[i + offset + p] = l - r;
        }
        now *= sum_e[bsf(~(uint)(s))];
      }
    }
  }
  void butterfly_inv(vector<mint>& a) {
    static mint g = _g; int n = (int)a.size(), h = ceil_pow2(n);
    static bool first = true; static mint sum_ie[30];
    if (first) {
      first = false; mint es[30], ies[30];
      int cnt2 = bsf(mod - 1); mint e = g.pow((mod - 1) >> cnt2), ie = ~e;
      for (int i = cnt2; i >= 2; i--) {
        es[i - 2] = e; ies[i - 2] = ie; e *= e; ie *= ie;
      }
      mint now = 1; rep(i, 0, cnt2-1) { sum_ie[i] = ies[i] * now; now *= es[i];}
    }
    for (int ph = h; ph >= 1; ph--) {
      int w = 1 << (ph - 1), p = 1 << (h - ph); mint inow = 1;
      rep(s, 0, w) {
        int offset = s << (h - ph + 1);
        for (int i = 0; i < p; i++) {
          auto l = a[i + offset], r = a[i + offset + p];
          a[i + offset] = l + r;
          a[i + offset + p] = (unsigned long long)(mod + l.x - r.x) * inow.x;
        }
        inow *= sum_ie[bsf(~(uint)(s))];
      }
    }
  }
  vector<mint> operator()(vector<mint> a, vector<mint> b) {
    int n = (int)a.size(), m = (int)b.size();
    if (!n || !m) return {};
    if (min(n, m) <= 60) {
      if (n < m) { swap(n, m); swap(a, b);}
      vector<mint> ans(n+m-1); rep(i,0,n)rep(j,0,m) ans[i + j] += a[i] * b[j];
      return ans;
    }
    int z = 1 << ceil_pow2(n+m-1);
    a.resize(z); butterfly(a); b.resize(z); butterfly(b);
    rep(i, 0, z) a[i] *= b[i];
    butterfly_inv(a); a.resize(n + m - 1); mint iz = ~(mint(z));
    rep(i, 0, n+m-1) a[i] *= iz;
    return a;
  }
} ntt;


//NTT
template<int mod=1012924417>
struct NTT {
  vector<int> rev, rts;
  int base, max_base, root;
  NTT(): base(1), rev{0, 1}, rts{0, 1} {
    assert(mod >= 3 && mod&1);
    auto tmp = mod-1;
    max_base = 0;
    while (tmp%2 == 0) { tmp >>= 1; max_base++; }
    root = 2;
    while (mod_pow(root, (mod-1)>>1) == 1) root++;
    assert(mod_pow(root, mod-1) == 1);
    root = mod_pow(root, (mod-1)>>max_base);
  }
  inline int mod_pow(int x, int n) {
    int res = 1;
    while (n > 0) {
      if (n&1) res = mul(res, x);
      x = mul(x, x);
      n >>= 1;
    }
    return res;
  }
  inline int inverse(int x) { return mod_pow(x, mod-2); }
  inline unsigned add(unsigned x, unsigned y) {
    x += y;
    if (x >= mod) x -= mod;
    return x;
  }
  inline unsigned mul(unsigned a, unsigned b) {
    return 1ull*a*b%(unsigned long long)mod;
  }
  void ensure_base(int nbase) {
    if (nbase <= base) return;
    rev.resize(1<<nbase);
    rts.resize(1<<nbase);
    for (int i = 0; i < (1<<nbase); i++) rev[i] = (rev[i>>1]>>1)+((i&1)<<(nbase-1));
    assert(nbase <= max_base);
    while (base < nbase) {
      int z = mod_pow(root, 1<<(max_base-1-base));
      for (int i = 1<<(base-1); i < (1<<base); i++) {
        rts[i<<1] = rts[i];
        rts[(i<<1)+1] = mul(rts[i], z);
      }
      base++;
    }
  }
  void ntt(vector<int>& a) {
    const int n = (int)a.size();
    assert((n&(n-1)) == 0);
    int zeros = __builtin_ctz(n);
    ensure_base(zeros);
    int shift = base - zeros;
    for (int i = 0; i < n; i++) if (i < (rev[i]>>shift)) swap(a[i], a[rev[i]>>shift]);
    for (int k = 1; k < n; k <<= 1) {
      for (int i = 0; i < n; i += 2*k) {
        for (int j = 0; j < k; j++) {
          int z = mul(a[i+j+k], rts[j+k]);
          a[i+j+k] = add(a[i+j], mod-z);
          a[i+j] = add(a[i+j], z);
        }
      }
    }
  }
  vector<int> multiply(vector<int> a, vector<int> b) {
    int need = (int)a.size() + (int)b.size() - 1;
    int nbase = 1;
    while ((1<<nbase) < need) nbase++;
    ensure_base(nbase);
    int sz = 1<<nbase;
    a.resize(sz, 0); b.resize(sz, 0);
    ntt(a); ntt(b);
    int inv_sz = inverse(sz);
    for (int i = 0; i < sz; i++) a[i] = mul(a[i], mul(b[i], inv_sz));
    reverse(a.begin()+1, a.end());
    ntt(a);
    a.resize(need);
    return a;
  }
};


//Notation systemp of base N
ll baseN_to_long(string s, int N){
  ll res = 0;
  for(int i = 0; i < s.size(); i++){
    if('A' <= s[i] && s[i] <= 'Z') res = res*N + s[i] - 'A' + 10;
    else res = res*N + s[i] - '0';
  }
  return res;
}

string long_to_baseN(ll n, int N){
  if(n == 0) return "0";
  string s = "0123456789ABCDEF";
  string res;
  while(n > 0){ res = s[n%N] + res; n /= N;}
  return res;
}

//Lazy_segtree

// RMQ, RAQ
using S = ll;
using F = ll;
S op(S a, S b){ return max(a, b);}
S e(){ return -INF;}
S mapping(F f, S x){ return f+x;}
F composition(F f, F g){ return f+g;}
F id(){ return 0;}

// RSQ, RAQ
struct S{ ll val; int size;};
using F = ll;
S op(S a, S b){ return {a.val + b.val, a.size + b.size};}
S e(){ return {0, 1};}
S mapping(F f, S x){ return {x.val + f*x.size, x.size};}
F composition(F f, F g){ return f + g;}
F id() { return 0;}

template<class S, S (*op)(S, S), S (*e)(), class F, S (*mapping)(F, S), F (*composition)(F, F), F (*id)()>
struct LSEG{
  int n, size, log = 0;
  vector<S> d;
  vector<F> lz;
  LSEG(): LSEG(0) {}
  explicit LSEG(int n): LSEG(vector<S>(n, e())) {}
  explicit LSEG(const vector<S>& v): n((int)v.size()){
    while((1<<log) < n) log++;
    size = 1<<log;
    d = vector<S>(2*size, e());
    lz = vector<F>(size, id());
    for(int i = 0; i < n; i++) d[size+i] = v[i];
    for(int i = size-1; i >= 1; i--) update(i);
  }
  void update(int k){ d[k] = op(d[2*k], d[2*k+1]);}
  void all_apply(int k, F f){
    d[k] = mapping(f, d[k]);
    if(k < size) lz[k] = composition(f, lz[k]);
  }
  void push(int k){
    all_apply(2*k, lz[k]);
    all_apply(2*k+1, lz[k]);
    lz[k] = id();
  }
  void set(int p, S x){
    assert(0 <= p && p < n);
    p += size;
    for(int i = log; i >= 1; i--) push(p>>i);
    d[p] = x;
    for(int i = 1; i <= log; i++) update(p>>i);
  }
  S get(int p){
    assert(0 <= p && p < n);
    p += size;
    for(int i = log; i >= 1; i--) push(p>>i);
    return d[p];
  }
  S prod(int l, int r){
    assert(0 <= l && l <= r && r <= n);
    if(l == r) return e();
    l += size;
    r += size;
    for(int i = log; i >= 1; i--){
      if(((l>>i)<<i) != l) push(l>>i);
      if(((r>>i)<<i) != r) push((r-1)>>i);
    }
    S sml = e(), smr = e();
    while(l < r){
      if(l&1) sml = op(sml, d[l++]);
      if(r&1) smr = op(d[--r], smr);
      l >>= 1; r >>= 1;
    }
    return op(sml, smr);
  }
  S all_prod() { return d[1];}
  void apply(int p, F f){
    assert(0 <= p && p < n);
    p += size;
    for(int i = log; i >= 1; i--) push(p>>i);
    d[p] = mapping(f, d[p]);
    for(int i = 1; i <= log; i++) update(p>>i);
  }
  void apply(int l, int r, F f){
    assert(0 <= l && l <= r && r <= n);
    if(l == r) return;
    l += size; r += size;
    for(int i = log; i >= 1; i--){
      if(((l>>i)<<i) != l) push(l>>i);
      if(((r>>i)<<i) != r) push((r-1)>>i);
    }
    {
      int l2 = l, r2 = r;
      while(l < r){
        if(l&1) all_apply(l++, f);
        if(r&1) all_apply(--r, f);
        l >>= 1; r >>= 1;
      }
      l = l2; r = r2;
    }
    for(int i = 1; i <= log; i++){
      if(((l>>i)<<i) != l) update(l>>i);
      if(((r>>i)<<i) != r) update((r-1)>>i);
    }
  }
  template<bool (*g)(S)>
  int max_right(int l){return max_right(l, [](S x){ return g(x);});}
  template<class G>
  int max_right(int l, G g){
    assert(0 <= l && l <= n);
    assert(g(e()));
    if(l == n) return n;
    l += size;
    for(int i = log; i >= 1; i--) push(l>>i);
    S sm = e();
    do{
      while((l&1) == 0) l >>= 1;
      if(!g(op(sm, d[l]))){
        while(l < size){
          push(l);
          l = 2*l;
          if(g(op(sm, d[l]))){
            sm = op(sm, d[l]);
            l++;
          }
        }
        return l - size;
      }
      sm = op(sm, d[l]);
      l++;
    }while((l&-l) != l);
    return n;
  }
  template<bool (*g)(S)>
  int min_left(int r){ return min_left(r, [](S x){ return g(x);});}
  template<class G>
  int min_left(int r, G g){
    assert(0 <= r && r <= n);
    assert(g(e()));
    if(r == 0) return 0;
    r += size;
    for(int i = log; i >= 1; i--) push((r-1)>>i);
    S sm = e();
    do{
      r--;
      while(r > 1 && (r&1)) r >>= 1;
      if(!g(op(d[r], sm))){
        while(r < size){
          push(r);
          r = 2*r + 1;
          if(g(op(d[r], sm))){
            sm = op(d[r], sm);
            r--;
          }
        }
        return r + 1 - size;
      }
      sm = op(d[r], sm);
    }while((r&-r) != r);
    return 0;
  }
};


//LCA
template<typename T>
struct treeLCA{
  int const MAX_LOG_V;
  vector<vector<edge<T>>> g;
  int root, vn;
  vector<vector<int>> parent;
  vector<int> depth;
  treeLCA(vector<vector<edge<T>>> &_g, int _r=0):
  MAX_LOG_V(35), g(_g), root(_r), vn((int)g.size()),
  parent(MAX_LOG_V, vector<int>(vn, 0)), depth(vn, -1)
  {depth[root] = 0; init(vn);}
  void dfs(int v, int p, int d){
    parent[0][v] = p;
    depth[v] = d;
    for(int i=0; i<g[v].size(); i++){
      if(depth[ g[v][i].to ] >= 0) continue;
      if(g[v][i].to != p) dfs(g[v][i].to, v, d+1);
    }
  }
  void init(int V) {
    dfs(root, -1, 0);
    for(int k=0; k+1 < MAX_LOG_V; k++){
      for(int v=0; v < V; v++) {
        if(parent[k][v] < 0) parent[k+1][v] = -1;
        else parent[k+1][v] = parent[k][parent[k][v]];
      }
    }
  }
  int lca(int u, int v) {
    if(depth[u] > depth[v]) swap(u, v);
    for(int k=0; k < MAX_LOG_V; k++){
      if((depth[v] - depth[u]) >> k & 1){
        v = parent[k][v];
      }
    }
    if(u == v) return u;
    for(int k=MAX_LOG_V - 1; k>=0; k--){
      if(parent[k][u] != parent[k][v]){
        u = parent[k][u];
        v = parent[k][v];
      }
    }
    return parent[0][u];
  }
  int dist(int u, int v){
    int anc = lca(u, v);
    return depth[u] + depth[v] - 2*depth[anc];
  }
};


//LCS
string LCS(string s, string t) {
  int n = (int)s.size(), m = (int)t.size();
  vector<vector<int>> dp(n+1, vector<int> (m+1));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      dp[i+1][j+1] = max(dp[i+1][j], dp[i][j+1]);
      if (s[i] == t[j]) dp[i+1][j+1] = dp[i][j] + 1;
    }
  }
  string res;
  while (n + m) {
    if (n && dp[n][m] == dp[n-1][m]) n--;
    else if (m && dp[n][m] == dp[n][m-1]) m--;
    else { n--; m--; res += s[n];}
  }
  reverse(res.begin(), res.end());
  return res;
}


//Levenshtein_distance
template<typename T=int>
T Levenshtein_distance(string s, string t, T INSERT_COST=1, T DELETE_COST=1, T CHANGE_COST=1) {
  int n = (int)s.size(), m = (int)t.size();
  vector<vector<int>> dp(n+1, vector<int> (m+1));
  for (int i = 0; i <= n; i++) dp[i][0] = i*INSERT_COST;
  for (int j = 0; j <= m; j++) dp[0][j] = j*INSERT_COST;
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= m; j++) {
      T D = dp[i-1][j] + DELETE_COST;
      T I = dp[i][j-1] + INSERT_COST;
      T C = dp[i-1][j-1] + (s[i-1]==t[j-1]?0:CHANGE_COST);
      dp[i][j] = min({D, I, C});
    }
  }
  return dp[n][m];
}


//LIS
// dp[i] := minimum element at length i
template<typename T>
vector<long long> LIS(vector<T>& a){
  const long long LINF = 1001002003004005006ll;
  int n = (int)a.size();
  vector<long long> dp(n+1, LINF);
  for (int i = 0; i < n; i++) *lower_bound(dp.begin(), dp.end(), a[i]) = a[i];
  return dp;
}


//Online_dynamic_connectivity
const int NODE_SIZE = 303030 * 10;
struct euler_tour_tree {
  using value_type = long long;
  using size_type = std::size_t;
  using node_index = std::int_least32_t;
  using vertex_index = std::int_least32_t;
  struct node;
  static struct node n[NODE_SIZE];
  static node_index ni;
  struct node {
    vertex_index s, d;
    node_index c[3];
    int sz;
    int flag;
    value_type val;
    value_type Sigma;
    node(): sz(1) {}
    inline node& operator[](size_type d) { return n[c[d]]; }
  };
  node_index new_edge(int s, int d, bool hi) {
    int i = ni++;
    int ri = ni++;
    n[i].s = n[ri].d = s;
    n[i].d = n[ri].s = d;
    n[i].sz = n[ri].sz = 0;
    n[i].flag = hi;
    return i;
  }
  static void fix(node_index i) {
    n[i].sz = (n[i].s == n[i].d) ? 1 : 0;
    if(n[i].c[0]) n[i].sz += n[i][0].sz;
    if(n[i].c[1]) n[i].sz += n[i][1].sz;
    n[i].flag &= 0b0101;
    n[i].flag |= n[i].flag << 1;
    if(n[i].c[0]) n[i].flag |= n[i][0].flag & 0b1010;
    if(n[i].c[1]) n[i].flag |= n[i][1].flag & 0b1010;
    n[i].Sigma = n[i].val;
    if(n[i].c[0]) n[i].Sigma += n[i][0].Sigma;
    if(n[i].c[1]) n[i].Sigma += n[i][1].Sigma;
  }
  static int child_dir(node_index i) {
    if(n[i].c[2]) {
      if(n[i][2].c[0] == i) { return 0; }
      else if(n[i][2].c[1] == i) { return 1; }
    }
    return 2;
  }
  static void rotate(node_index x, size_type dir) {
    node_index p = n[x].c[2];
    int x_dir = child_dir(x);
    node_index y = n[x].c[dir^1];
    if(n[y].c[dir]) n[y][dir].c[2] = x;
    n[x].c[dir^1] = n[y].c[dir];
    n[n[x].c[2] = y].c[dir] = x;
    n[y].c[2] = p;
    if(x_dir < 2) n[p].c[x_dir] = y;
    if(n[x].c[dir^1]) fix(n[x].c[dir^1]);
    fix(x);
  }
  static void splay(node_index i) {
    int i_dir;
    int j_dir;
    while((i_dir = child_dir(i)) < 2) {
      node_index j = n[i].c[2];
      if((j_dir = child_dir(j)) < 2) {
        node_index k = n[j].c[2];
        if(i_dir == j_dir) { rotate(k, j_dir^1); rotate(j, i_dir^1);}
          else { rotate(j, i_dir^1); rotate(k, j_dir^1);}
      }
      else rotate(j, i_dir ^ 1);
    }
    fix(i);
  }
  static node_index merge_back(node_index l, node_index r) {
    if(!l) return r;
    if(!r) return l;
    while(n[l].c[1]) l = n[l].c[1];
    splay(l);
    n[n[r].c[2] = l].c[1] = r;
    fix(l);
    return l;
  }
  static std::pair<node_index, node_index> split(node_index i) {
    splay(i);
    node_index l = n[i].c[0];
    n[i].c[0] = n[l].c[2] = 0;
    fix(i);
    return { l, i };
  }
  static void reroot(node_index v) {
    auto p = split(v);
    merge_back(p.second, p.first);
    splay(v);
  }
  static bool same_root(node_index i, node_index j) {
    if(i) splay(i);
    if(j) splay(j);
    while(n[i].c[2]) i = n[i].c[2];
    while(n[j].c[2]) j = n[j].c[2];
    return i == j;
  }
  node_index n_start;
  std::unordered_map<long long, node_index> emp;
  euler_tour_tree() {}
  euler_tour_tree(int N): n_start(ni) {
    ni += N;
    for(int i = 0; i < N; i++) {
      n[i + n_start].s = n[i + n_start].d = i;
    }
  }
  bool edge_exist(vertex_index x, vertex_index y) {
    if(x > y) std::swap(x, y);
    return emp.count(((long long)x << 32) | (long long)y);
  }
  void link(vertex_index x, vertex_index y, bool hi) {
    if(x > y) std::swap(x, y);
    int ei = new_edge(x, y, hi);
    emp[((long long)x << 32) | (long long)y] = ei;
    x += n_start;
    y += n_start;
    reroot(x);
    reroot(y);
    n[n[x].c[2] = ei].c[0] = x;
    n[n[y].c[2] = ei].c[1] = y;
    fix(ei);
    merge_back(ei, ei + 1);
  }
  void cut(vertex_index x, vertex_index y) {
    if(x > y) std::swap(x, y);
    auto iter = emp.find(((long long)x << 32) | (long long)y);
    int ei = iter->second;
    int rei = ei + 1;
    emp.erase(iter);
    auto p = split(ei);
    if(p.first && same_root(p.first, rei)) {
      auto q = split(rei);
      node_index left = q.first;
      node_index center = n[q.second].c[1];
      node_index right = n[p.second].c[1];
      n[center].c[2] = 0;
      n[right].c[2] = 0;
      merge_back(left, right);
    }
    else {
      splay(ei);
      ei = n[ei].c[1];
      n[ei].c[2] = 0;
      auto q = split(rei);
      splay(p.first);
      node_index left = p.first;
      node_index center = q.first;
      node_index right = n[q.second].c[1];
      n[right].c[2] = 0;
      merge_back(left, right);
    }
  }

  bool same_tree(vertex_index x, vertex_index y) {
    return same_root(x + n_start, y + n_start);
  }
  int tree_size(vertex_index x) {
    x += n_start;
    splay(x);
    return n[x].sz;
  }
  void subedge_set(vertex_index x, bool val) {
    x += n_start;
    splay(x);
    if(val) n[x].flag |= (0b0100);
    else n[x].flag &= ~(0b0100);
    fix(x);
  }
  void add_val(vertex_index x, value_type val) {
    x += n_start;
    splay(x);
    n[x].val += val;
    fix(x);
  }
  value_type tree_sum(vertex_index x) {
    x += n_start;
    splay(x);
    return n[x].Sigma;
  }
  template<class Func>
  void hilevel_edges(vertex_index v, Func f) {
    node_index i = v + n_start;
    splay(i);
    while(i && (n[i].flag & 0b0010)) {
      while(1) {
        if(n[i].flag & 0b0001) {
          f(n[i].s, n[i].d);
          splay(i);
          n[i].flag &= ~(0b0001);
          fix(i);
          break;
        }
        else if(n[i].c[0] && (n[i][0].flag & 0b0010)) i = n[i].c[0];
        else i = n[i].c[1];
      }
    }
  }
  template<class Func>
  int subedges(vertex_index v, Func f) {
    node_index i = v + n_start;
    splay(i);
    while(i && (n[i].flag & 0b1000)) {
      while(1) {
        if(n[i].flag & 0b0100) {
          if(f(n[i].s)) {
            return 1;
          }
          splay(i);
          break;
        }
        else if(n[i].c[0] && (n[i][0].flag & 0b1000)) i = n[i].c[0];
        else i = n[i].c[1];
      }
    }
    return 0;
  }
};

euler_tour_tree::node_index euler_tour_tree::ni = 1;
euler_tour_tree::node euler_tour_tree::n[NODE_SIZE];

struct online_dynamic_connectivity {
  int N;
  std::vector<euler_tour_tree> ett;
  std::vector<std::vector<std::unordered_set<int>>> E;
  online_dynamic_connectivity(int N): N(N) {
    ett.emplace_back(N);
    E.emplace_back(N);
  }
  void link(int x, int y) {
    if(ett[0].same_tree(x, y)) {
      if(E[0][x].size() == 0) ett[0].subedge_set(x, 1);
      if(E[0][y].size() == 0) ett[0].subedge_set(y, 1);
      E[0][x].insert(y);
      E[0][y].insert(x);
    }
    else {
      ett[0].link(x, y, true);
    }
  }
  void replace(int x, int y, int level) {
    for(int k = 0; k < level; k++) {
      ett[k].cut(x, y);
    }
    for(int k = level; k-- > 0;) {
      if(ett[k].tree_size(x) > ett[k].tree_size(y)) std::swap(x, y);
      ett[k].hilevel_edges(x, [&](int s, int d) { ett[k + 1].link(s, d, true); });
      int res = ett[k].subedges(x, [&](int s) {
        for(auto iter = E[k][s].begin(); iter != E[k][s].end(); ) {
          int d = *iter;
          iter = E[k][s].erase(iter);
          E[k][d].erase(s);
          if(E[k][s].size() == 0) ett[k].subedge_set(s, 0);
          if(E[k][d].size() == 0) ett[k].subedge_set(d, 0);
          if(ett[k].same_tree(s, d)) {
            if(E[k + 1][s].size() == 0) ett[k + 1].subedge_set(s, 1);
            if(E[k + 1][d].size() == 0) ett[k + 1].subedge_set(d, 1);
            E[k + 1][s].insert(d);
            E[k + 1][d].insert(s);
          }
          else {
            for(int kk = k + 1; kk-- > 0;) {
              ett[kk].link(s, d, kk == k);
            }
            return 1;
          }
        }
        return 0;
        });
      if(res) return;
    }
  }
  void cut(int x, int y) {
    for(int k = 0; k < ett.size(); k++) {
      if(E[k][x].count(y)) {
        E[k][x].erase(y);
        E[k][y].erase(x);
        if(E[k][x].size() == 0) ett[k].subedge_set(x, 0);
        if(E[k][y].size() == 0) ett[k].subedge_set(y, 0);
        return;
      }
    }
    for(int k = (int)ett.size(); k-- > 0;) {
      if(ett[k].edge_exist(x, y)) {
        if(k + 1 == ett.size()) {
          ett.emplace_back(N);
          E.emplace_back(N);
        }
        replace(x, y, k + 1);
      }
    }
  }
  void add_val(int x, long long val) { ett[0].add_val(x, val);}
  int size(int x) { return ett[0].tree_size(x);}
  long long sum(int x) { return ett[0].tree_sum(x);}
  bool same(int x, int y) { return ett[0].same_tree(x, y);}
};


//Perm
struct Perm : vector<int> {
#define n (int)(size())
#define p (*this)
  Perm(int _n): vector<int>(_n) { iota(begin(), end(), 0);}
  template<class...Args> Perm(Args...args): vector<int>(args...) {}
  Perm(initializer_list<int> a): vector<int>(a.begin(),a.end()) {}
  Perm operator+(const Perm& a) const {
    Perm r(n);
    for (int i = 0; i < n; ++i) r[i] = p[a[i]];
    return r;
  }
  Perm& operator+=(const Perm& a) { return *this = (*this)+a;}
  Perm operator-() const {
    Perm r(n);
    for (int i = 0; i < n; ++i) r[p[i]] = i;
    return r;
  }
  Perm operator-(const Perm& a) const { return *this + -a;}
  Perm& operator-=(const Perm& a) { return *this += -a;}
  // next permutation
  bool operator++() { return next_permutation(begin(),end());}
#undef n
#undef p
};


//PERMUTATION_SUM
ll permutation_sum(ll l, ll r, ll n){ return (l+r)*n/2;}
mint mod_permutation_sum(ll l, ll r, ll n){ return mint(l+r)*n/2;}
ll permutation_sum2(ll a, ll d, ll n){ return (2*a + (n-1)*d)*n/2;}
mint mod_permutation2(ll a, ll d, ll n){ return mint(2*a + (n-1)*d)*n/2;}
ll permutation_sum3(ll l, ll r, ll d=1){ return (l+r)*(r-l+1)/2;}
mint mod_permutation_sum3(ll l, ll r, ll d=1){ return mint(l+r)*(r-l+1)/2;}


//POWER_MODver. N^k % MOD
ll mod_pow(ll n, ll k){
  ll res = 1;
  for(; k > 0; k >>= 1){
    if(k&1) res = (res*n)%mod;
    n = (n*n)%mod;
  }
  return res;
}


//Prim
template<typename T>
T prim(vector<vector<edge<T>>> &g){
  int n = (int)g.size(); T ans = 0;
  vector<bool> used(n);
  priority_queue<edge<T>, vector<edge<T>>, greater<edge<T>>> q;
  q.push(edge<T>(0, 0));
  while(!q.empty()){
    edge<T> e = q.top(); q.pop();
    if(used[e.to]) continue;
    used[e.to] = true;
    ans += e.cost;
    for(int i = 0; i < g[e.to].size(); i++) q.emplace(g[e.to][i]);
  }
  return ans;
}


//PRIME_FACTOR
vector<pair<ll, ll>> prime_factor(ll n){
  vector<pair<ll, ll>> pf;
  // if(n == 1) pf.push_back(P(1, 1));
  for(ll i = 2; i*i <= n; i++){
    if(n%i == 0){
      int cnt = 0;
      while(n%i == 0){
        n /= i;
        cnt++;
      }
      pf.push_back(pair<ll, ll>(i, cnt));
    }
  }
  if(n != 1) pf.push_back(pair<ll, ll>(n, 1));
  return pf;
}


//PRIME_FACTOR_fastver.(Pollard’s rho algorithm, *randomized)
vector<pair<long long, long long>> prime_factor(long long n) {
  auto pollard = [&](long long n) -> long long {
    if (n%2 == 0) return 2;
    if (is_prime(n)) return n;
    auto f = [&](long long x) -> long long {
      return (__int128_t(x)*x + 1)%n;
    };
    long long step = 0;
    while (true) {
      long long x = ++step, y = f(x);
      while (true) {
        long long p = gcd(y - x + n, n);
        if (p == 0 || p == n) break;
        if (p != 1) return p;
        x = f(x); y = f(f(y));
      }
    }
  };
  auto prime_factorize = [&](auto f, long long n) -> vector<long long> {
    if (n == 1) return {};
    long long p = pollard(n);
    if (p == n) return {p};
    vector<long long> left = f(f, p), right = f(f, n/p);
    left.insert(left.end(), right.begin(), right.end());
    sort(left.begin(), left.end());
    return left;
  };
  vector<pair<long long, long long>> pf;
  auto res = prime_factorize(prime_factorize, n);
  res.erase(unique(res.begin(), res.end()), res.end());
  for (auto i : res) {
    long long cnt = 0;
    while(n%i == 0){
      n /= i;
      cnt++;
    }
    pf.push_back(make_pair(i, cnt));
  }
  return pf;
}


//Range Maximun Query
template<typename T>
struct RMQ{
  int n;
  vector<ll> dat;
  RMQ(){}
  RMQ(int n_) {init(n_);}
  void init(int n_){
    n = 1;
    while(n < n_) n *= 2;
    dat.resize(n*2 - 1, -LINF);
  }
  void update(T k, T a){
    k += n - 1;
    dat[k] = a;
    while(k > 0){
      k = (k-1)/2;
      dat[k] = max(dat[k*2+1], dat[k*2+2]);
    }
  }
  T query(T a, T b, T k, T l, T r){
    if(r <= a || b <= l) return -LINF;
    if(a <= l && r <= b) return dat[k];
    else{
      T vl = query(a, b, k*2 + 1, l, (l + r)/2);
      T vr = query(a, b, k*2 + 2, (l + r)/2, r);
      return max(vl, vr);
    }
  }
  T query(T a, T b){
    return query(a, b, 0, 0, n);
  }
};


//Range Minimum Query
template<typename T>
struct RMQ{
  int n;
  vector<ll> dat;
  RMQ(){}
  RMQ(int n_) {init(n_);}
  void init(int n_){
    n = 1;
    while(n < n_) n *= 2;
    dat.resize(2*n - 1, LINF);
  }
  void update(T k, T a){
    k += n - 1;
    dat[k] = a;
    while(k > 0){
      k = (k-1)/2;
      dat[k] = min(dat[k*2+1], dat[k*2+2]);
    }
  }
  T query(T a, T b, T k, T l, T r){
    if(r <= a || b <= l) return LINF;
    if(a <= l && r <= b) return dat[k];
    else{
      T vl = query(a, b, k*2 + 1, l, (l + r)/2);
      T vr = query(a, b, k*2 + 2, (l + r)/2, r);
      return min(vl, vr);
    }
  }
  T query(T a, T b){
    return query(a, b, 0, 0, n);
  }
};


//Range Minimum Query -- pair.ver
template<typename T>
struct RMQ{
  int n;
  vector<pair<ll, int>> dat;
  RMQ(){}
  RMQ(int n_) {init(n_);}
  void init(int n_){
    n = 1;
    while(n < n_) n *= 2;
    dat.resize(n*2 - 1);
    for(int i = 0; i < n*2 - 1; i++){
      dat[i].first = LINF;
      dat[i].second = i - n + 1;
    }
  }
  void update(T k, T a){
    k += n - 1;
    dat[k].first = a;
    while(k > 0){
      k = (k-1)/2;
      dat[k] = dat[k*2+1].first < dat[k*2+2].first ? dat[k*2+1] : dat[k*2+2];
    }
  }
  pair<ll, ll> query(T a, T b, T k, T l, T r){
    if(r <= a || b <= l) return make_pair(LINF, -1);
    if(a <= l && r <= b) return dat[k];
    else{
      pair<ll, int> vl = query(a, b, k*2 + 1, l, (l + r)/2);
      pair<ll, int> vr = query(a, b, k*2 + 2, (l + r)/2, r);
      return vl.first < vr.first ? vl : vr;
    }
  }
  T query(T a, T b){
    pair<ll, int> q = query(a, b, 0, 0, n);
    // update(q.second, LINF);
    return q.first;
  }
};


//RecSum
template<typename T>
struct RecSum{
  int h, w;
  vector<vector<T>> d;
  RecSum(const vector<vector<T>> &v): h(v.size()), w(v[0].size()), d(h+1, {0}) {
    for(int i=0;i<h;i++)for(int j=0;j<w;j++) d[i+1].push_back(v[i][j]+d[i+1][j]);
    for(int j=0;j<w;j++)for(int i=0;i<h;i++) d[i+1][j+1] += d[i][j+1];
  }
  // return [x1, x2) * [y1, y2)
  T get(int x1, int x2, int y1, int y2) {
    return d[y2][x2]-d[y1][x2]-d[y2][x1]+d[y1][x1];
  }
};


//Rerooting
template<typename T>
struct Rerooting{
  struct DP{
    ll mx; // Maximum distance among all vertices
    ll ans; // Sum of the minimum distances to all vertices
    int t;
    DP(ll mx=-1, ll ans=0, int t=0): mx(mx), ans(ans), t(t) {}
    DP operator+(const DP& a) const{
      return DP(max(mx, a.mx), ans+a.ans, t+a.t);
    }
    DP addRoot() const{
      return DP(mx+1, ans+t, t+1);
    }
  };
  int n;
  vector<vector<int>> to;
  vector<vector<DP>> dp;
  vector<DP> ans;
  Rerooting(int n=0): n(n), to(n), dp(n), ans(n) {}
  Rerooting(vector<vector<edge<T>>> &s): n((int)s.size()), to(n), dp(n), ans(n) {
    for(int i = 0; i < n; i++) for(edge j : s[i]) addEdge(i, j.to);
    init();
  }
  void addEdge(int a, int b){ to[a].push_back(b);}
  void init(){ dfs(0); bfs(0);}
  DP dfs(int v, int p=-1){
    DP dpSum;
    dp[v] = vector<DP>(to[v].size());
    for(int i = 0; i < to[v].size(); i++){
      int u = to[v][i];
      if(u == p) continue;
      dp[v][i] = dfs(u, v);
      dpSum = dpSum + dp[v][i];
    }
    return dpSum.addRoot();
  }
  void bfs(int v, const DP& dpP=DP(), int p=-1){
    int deg = (int)to[v].size();
    vector<DP> dpSumL(deg+1), dpSumR(deg+1);
    for(int i = 0; i < deg; i++) if(to[v][i] == p) dp[v][i] = dpP;
    for(int i = 0; i < deg; i++) dpSumL[i+1] = dpSumL[i] + dp[v][i];
    for(int i = deg-1; i >= 0; i--) dpSumR[i] = dpSumR[i+1] + dp[v][i];
    ans[v] = dpSumL[deg].addRoot();
    for(int i = 0; i < deg; i++){
      int u = to[v][i];
      if(u == p) continue;
      DP d = dpSumL[i] + dpSumR[i+1];
      bfs(u, d.addRoot(), v);
    }
  }
};


//RollingHash
struct RollingHash {
  static const int base1 = 1007, base2 = 2009;
  static const int mod1 = 1000000007, mod2 = 1000000009;
  vector<long long> hash1, hash2, power1, power2;
  RollingHash(const string &s) {
    int n = (int)s.size();
    hash1.assign(n+1, 0); hash2.assign(n+1, 0);
    power1.assign(n+1, 1); power2.assign(n+1, 1);
    for (int i = 0; i < n; i++) {
      hash1[i+1] = (hash1[i]*base1 + s[i])%mod1;
      hash2[i+1] = (hash2[i]*base2 + s[i])%mod2;
      power1[i+1] = (power1[i]*base1)%mod1;
      power2[i+1] = (power2[i]*base2)%mod2;
    }
  }
  // get hash value of s[l:r]
  inline long long get(int l, int r) const {
    long long res1 = hash1[r] - hash1[l]*power1[r-l]%mod1;
    long long res2 = hash2[r] - hash2[l]*power2[r-l]%mod2;
    if (res1 < 0) res1 += mod1;
    if (res2 < 0) res2 += mod2;
    return res1*mod2 + res2;
  }
  
  // get hash value of s
  inline long long get() const {
    return hash1.back()*mod2 + hash2.back();
  }
  
  // get lcp of s[a:] and s[b:]
  inline int getLCP(int a, int b) const {
    int len = min((int)hash1.size()-a, (int)hash1.size()-b);
    int low = 0, high = len;
    while (high - low > 1) {
      int mid = (low + high)>>1;
      if (get(a, a+mid) != get(b, b+mid)) high = mid;
      else low = mid;
    }
    return low;
  }
  
  // get lcp of s[a:] and t[b:]
  inline int getLCP(const RollingHash &t, int a, int b) const {
    int len = min((int)hash1.size()-a, (int)hash1.size()-b);
    int low = 0, high = len;
    while (high - low > 1) {
      int mid = (low + high)>>1;
      if (get(a, a+mid) != t.get(b, b+mid)) high = mid;
      else low = mid;
    }
    return low;
  }
};


//SCC_GRAPH
template<typename T=int>
struct scc_graph {
  int n;
  vector<pair<int, int>> edges;
  struct csr {
    vector<int> start, elist;
    csr(int n, const vector<pair<int, int>>& edges): start(n+1), elist(edges.size()) {
      for (auto e : edges) start[e.first + 1]++;
      for (int i = 1; i <= n; i++) start[i] += start[i - 1];
      auto counter = start;
      for (auto e : edges) elist[counter[e.first]++] = e.second;
    }
  };
  scc_graph(int _n): n(_n){}
  scc_graph(vector<vector<edge<T>>> &s): n((int)s.size()) {
    for (int i = 0; i < n; i++) for (int j = 0; j < s[i].size(); j++) {
      add_edge(i, s[i][j].to);
    }
  }
  
  int num_vertices() { return n; }
  
  void add_edge(int from, int to) { edges.push_back({from, to}); }
  
  // @return pair of (# of scc, scc id)
  pair<int, vector<int>> scc_ids() {
    auto g = csr(n, edges);
    int now_ord = 0, group_num = 0;
    vector<int> visited, low(n), ord(n, -1), ids(n);
    visited.reserve(n);
    auto dfs = [&](auto self, int v) -> void {
      low[v] = ord[v] = now_ord++;
      visited.push_back(v);
      for (int i = g.start[v]; i < g.start[v + 1]; i++) {
        auto to = g.elist[i];
        if (ord[to] == -1) {
          self(self, to);
          low[v] = min(low[v], low[to]);
        } else {
          low[v] = min(low[v], ord[to]);
        }
      }
      if (low[v] == ord[v]) {
        while (true) {
          int u = visited.back();
          visited.pop_back();
          ord[u] = n;
          ids[u] = group_num;
          if (u == v) break;
        }
        group_num++;
      }
    };
    for (int i = 0; i < n; i++) if (ord[i] == -1) dfs(dfs, i);
    for (auto& x : ids) x = group_num - 1 - x;
    return {group_num, ids};
  }
  
  vector<vector<int>> scc() {
    auto ids = scc_ids();
    int group_num = ids.first;
    vector<int> counts(group_num);
    for (auto x : ids.second) counts[x]++;
    vector<vector<int>> groups(ids.first);
    for (int i = 0; i < group_num; i++) groups[i].reserve(counts[i]);
    for (int i = 0; i < n; i++) groups[ids.second[i]].push_back(i);
    return groups;
  }
};


//SegmentTree
template<class S, S (*op)(S, S), S (*e)()>
struct SEG{
  int n, size, log = 0;
  vector<S> d;
  SEG(): SEG(0) {}
  explicit SEG(int n): SEG(vector<S>(n, e())) {}
  explicit SEG(const vector<S>& v): n((int)v.size()){
    while((1<<log) < n) log++;
    size = 1<<log;
    d = vector<S>(2*size, e());
    for(int i = 0; i < n; i++) d[size+i] = v[i];
    for(int i = size-1; i >= 1; i--) update(i);
  }
  void update(int k) { d[k] = op(d[2*k], d[2*k+1]);}
  void set(int p, S x){
    assert(0 <= p && p < n);
    p += size;
    d[p] = x;
    for(int i = 1; i <= log; i++) update(p>>i);
  }
  S get(int p) const{
    assert(0 <= p && p < n);
    return d[p+size];
  }
  S prod(int l, int r) const{
    assert(0 <= l && l <= r && r <= n);
    S sml = e(), smr = e();
    l += size;
    r += size;
    while(l < r){
      if(l&1) sml = op(sml, d[l++]);
      if(r&1) smr = op(d[--r], smr);
      l >>= 1; r >>= 1;
    }
    return op(sml, smr);
  }
  S all_prod() const{ return d[1];}
  template<bool (*f)(S)>
  int max_right(int l) const{ return max_right(l, [](S x) {return f(x);});}
  template<class F>
  int max_right(int l, F f) const{
    assert(0 <= l && l <= n);
    assert(f(e()));
    if(l == n) return n;
    l += size;
    S sm = e();
    do{
      while((l&1) == 0) l >>= 1;
      if(!f(op(sm, d[l]))){
        while(l < size){
          l = 2*l;
          if(f(op(sm, d[l]))){
            sm = op(sm, d[l]);
            l++;
          }
        }
        return l - size;
      }
      sm = op(sm, d[l]);
      l++;
    }while((l&-l) != l);
    return n;
  }
  template<bool (*f)(S)>
  int min_left(int r) const{ return min_left(r, [](S x){ return f(x);});}
  template<class F>
  int min_left(int r, F f) const{
    assert(0 <= r && r <= n);
    assert(f(e()));
    if(r == 0) return 0;
    r += size;
    S sm = e();
    do{
      r--;
      while(r > 1 && (r&1)) r >>= 1;
      if(!f(op(d[r], sm))){
        while(r < size){
          r = 2*r + 1;
          if(f(op(d[r], sm))){
            sm = op(d[r], sm);
            r--;
          }
        }
        return r + 1 - size;
      }
      sm = op(d[r], sm);
    }while((r&-r) != r);
    return 0;
  }
};



//Topological sort - Tarjan.ver
template<typename T>
void dfs_t(vector<vector<edge<T>>>& g, int v, vector<bool>& used, vector<int>& ans){
  if(used[v]) return;
  used[v] = true;
  for(edge<T> u : g[v]) dfs_t(g, u.to, used, ans);
  ans.emplace_back(v);
}
template<typename T>
vector<int> tsort(vector<vector<edge<T>>>& g){
  int n = (int)g.size();
  vector<bool> used(n);
  vector<int> res;
  for(int v = 0; v < n; v++) dfs_t(g, v, used, res);
  reverse(res.begin(), res.end());
  return res;
}


//Topological sort - Kahn.ver (Closed circuit detectable)
template<typename T>
vector<int> tsort(vector<vector<edge<T>>>& g){
  vector<int> res;
  int n = (int)g.size();
  vector<int> ind(n);
  queue<int> q;
  for(int i = 0; i < n; i++) for(edge<T> e : g[i]) ind[e.to]++;
  for(int i = 0; i < n; i++) if(ind[i] == 0) q.push(i);
  while(!q.empty()){
    int now = q.front(); q.pop();
    res.emplace_back(now);
    for(edge<T> e : g[now]){
      ind[e.to]--;
      if(ind[e.to] == 0) q.push(e.to);
    }
  }
  return res;
}


//TreeHeight
template<typename T>
vector<T> treeHeight(vector<vector<edge<T>>> &g){
  int n = (int)g.size();
  vector<T> v1, v2, v3, ret(n);
  int p1, p2;
  v1 = dijkstra<T>(g, 0);
  p1 = max_element(v1.begin(), v1.end()) - v1.begin();
  v2 = dijkstra<T>(g, p1);
  p2 = max_element(v2.begin(), v2.end()) - v2.begin();
  v3 = dijkstra<T>(g, p2);
  for(int i=0; i<n; i++) ret[i] = max(v2[i], v3[i]);
  return ret;
}


//Trie
template<size_t X>
struct Trie {
  struct Node {
    array<int, X> nxt;
    int idx, cnt;
    Node(): idx(-1), cnt(0) { fill(nxt.begin(), nxt.end(), -1);}
  };
  using F = function<int(char)>;
  vector<Node> vs;
  F conv;
  Trie(F conv): conv(conv) { vs.emplace_back();}
  Trie(char start): Trie([=](char a) { return a-start;}){}
  inline int &next(int i, int j) { return vs[i].nxt[j];}
  void add(const string &s) {
    int pos = 0;
    for (char c : s) {
      int k = conv(c);
      if (~next(pos, k)) { pos=next(pos,k); continue;}
      int npos = (int)vs.size();
      next(pos, k) = npos;
      vs.emplace_back();
      pos = npos;
    }
    pos = 0;
    for (char c : s) {
      int k = conv(c);
      pos = next(pos, k);
      vs[pos].cnt++;
    }
  }
  int calc(const string &s) {
    int pos = 0, dep = 0, res = 0;
    for (char c : s) {
      int k = conv(c);
      pos = next(pos, k);
      dep++;
      if (vs[pos].cnt >= 2) res = max(res, dep);
    }
    return res;
  }
  int find(const string &s) {
    int pos = 0;
    for (char c : s) {
      int k = conv(c);
      if (next(pos, k) < 0) return -1;
      pos = next(pos, k);
    }
    return pos;
  }
  int move(int pos, char c) {
    assert(pos < (int)vs.size());
    return pos < 0 ? -1 : next(pos, conv(c));
  }
  int size() { return (int)vs.size();}
  int idx(int pos) { return pos < 0 ? -1 : vs[pos].idx;}
  vector<int> idxs(int pos) { return pos < 0 ? vector<int>() : vs[pos].idxs;}
};


//UNIONFIND
template<typename T>
struct UnionFind {
  int num;
  vector<T> d;
  UnionFind(int n=0): num(n), d(n,-1) {}
  T find(int x) {
    if (d[x] < 0) return x;
    return d[x] = find(d[x]);
  }
  bool unite(int x, int y) {
    x = find(x); y = find(y);
    if (x == y) return false;
    if (d[x] > d[y]) swap(x,y);
    d[x] += d[y];
    d[y] = x;
    num--;
    return true;
  }
  bool same(int x, int y) { return find(x) == find(y);}
  T size(int x) { return -d[find(x)];}
  int count() { return num;}
};


//Vector
struct Vector {
  double x, y;
  const double eps = 1e-9;
  Vector(double x=0, double y=0): x(x), y(y) {}
  Vector& operator+=(const Vector& v) { x += v.x; y += v.y; return *this;}
  Vector operator+(const Vector& v) const { return Vector(*this) += v;}
  Vector& operator-=(const Vector& v) { x -= v.x; y -= v.y; return *this;}
  Vector operator-(const Vector& v) const { return Vector(*this) -= v;}
  Vector& operator*=(double s) { x *= s; y *= s; return *this;}
  Vector operator*(double s) const { return Vector(*this) *= s;}
  Vector& operator/=(double s) { x /= s; y /= s; return *this;}
  Vector operator/(double s) const { return Vector(*this) /= s;}
  double dot(const Vector& v) const { return x*v.x + y*v.y;}
  double cross(const Vector& v) const { return x*v.y - v.x*y;}
  double norm2() const { return x*x + y*y;}
  double norm() const { return sqrt(norm2());}
  Vector rotate90() const { return Vector(y, -x);}
  int ort() const { // orthant
    if (abs(x) < eps && abs(y) < eps) return 0;
    if (y > 0) return x>0 ? 1 : 2;
    else return x>0 ? 4 : 3;
  }
  bool operator<(const Vector& v) const {
    int o = ort(), vo = v.ort();
    if (o != vo) return o < vo;
    return cross(v) > 0;
  }
};
istream& operator>>(istream& is, Vector& v) { is >> v.x >> v.y; return is;}
ostream& operator<<(ostream& os, const Vector& v) { os<<"("<<v.x<<","<<v.y<<")"; return os;}



//WAVELETMATRIX
struct FullyIndexableDictionary{
  int len, blk;
  vector<unsigned> bit;
  vector<int> sum;
  FullyIndexableDictionary() {}
  FullyIndexableDictionary(int len): len(len), blk((len+31)>>5), bit(blk, 0), sum(blk, 0) {}
  void set(int k){ bit[k>>5] |= 1u<<(k&31);}
  void build(){
    sum[0] = 0;
    for(int i=1;i<blk;i++) sum[i]=sum[i-1]+__builtin_popcount(bit[i-1]);
  }
  bool operator[](int k) const{ return bool((bit[k>>5]>>(k&31))&1);}
  int rank(int k){ return sum[k>>5]+__builtin_popcount(bit[k>>5]&((1u<<(k&31))-1));}
  int rank(bool v, int k){ return (v?rank(k):k-rank(k));}
  int select(bool v, int k){
    if(k<0 || rank(v, len) <= k) return -1;
    int l = 0, r = len;
    while(r - l > 1){
      int m = (l+r)>>1;
      if(rank(v, m) >= k+1) r = m;
      else l = m;
    }
    return r - 1;
  }
  int select(bool v, int i, int l){ return select(v, i+rank(v, l));}
};

template<class T, int MAXLOG=30>
struct WaveletMatrix{
  int len;
  FullyIndexableDictionary mat[MAXLOG];
  int zs[MAXLOG], buff1[MAXLOG], buff2[MAXLOG];
  vector<vector<T>> acc;
  static const T npos=-1;
  WaveletMatrix(vector<T> data){
    len = (int)data.size();
    acc = vector<vector<T>>(MAXLOG, vector<T>(len+1));
    vector<T> ls(len), rs(len);
    for(int dep=0;dep<MAXLOG;dep++){
      mat[dep] = FullyIndexableDictionary(len+1);
      int p = 0, q = 0;
      for(int i=0;i<len;i++){
        bool k = (data[i]>>(MAXLOG-(dep+1)))&1;
        if(k) rs[q++] = data[i], mat[dep].set(i);
        else ls[p++] = data[i];
      }
      zs[dep] = p;
      mat[dep].build();
      for(int i=0;i<len;i++){
        if(!mat[dep][i]) acc[dep][i+1]=data[i];
        acc[dep][i+1]+=acc[dep][i];
      }
      swap(ls, data);
      for(int i=0;i<q;i++) data[p+i] = rs[i];
    }
  }
  
  T access(int k){
    T res = 0;
    for(int dep = 0; dep < MAXLOG; dep++){
      bool bit = mat[dep][k];
      res = (res<<1)|bit;
      k = mat[dep].rank(bit, k) + zs[dep]*dep;
    }
    return res;
  }
  
  // return the number of v in [0, k)
  int rank(T v, int k){
    int l = 0, r = k;
    for(int dep = 0; dep < MAXLOG; dep++){
      buff1[dep] = l; buff2[dep] = r;
      bool bit = (v>>(MAXLOG-(dep+1)))&1;
      l = mat[dep].rank(bit, l) + zs[dep]*bit;
      r = mat[dep].rank(bit, r) + zs[dep]*bit;
    }
    return r - l;
  }
  
  // return the position of k-th v
  int select(T v, int k){
    rank(v, len);
    for(int dep = MAXLOG-1; dep >= 0; dep--){
      bool bit = (v>>(MAXLOG-(dep+1)))&1;
      k = mat[dep].select(bit, k, buff1[dep]);
      if(k >= buff2[dep] || k < 0) return -1;
      k -= buff1[dep];
    }
    return k;
  }
  int select(T v, int k, int l){ return select(v, k+rank(v, l));}
  
  // return k-th largest value in [l, r)
  T quantile(int l, int r, int k){
    if(r - l <= k || k < 0) return -1;
    T res = 0;
    for(int dep = 0; dep < MAXLOG; dep++){
      int p = mat[dep].rank(1, l);
      int q = mat[dep].rank(1, r);
      if(q - p > k){
        l = p + zs[dep];
        r = q + zs[dep];
        res |= T(1)<<(MAXLOG-(dep+1));
      }
      else{
        k -= (q - p);
        l -= p;
        r -= q;
      }
    }
    return res;
  }
  T rquantile(int l, int r, int k){ return quantile(l, r, r-l-k-1);}
  
  int freq_dfs(int d, int l, int r, T val, T a, T b){
    if(l == r) return 0;
    if(d == MAXLOG) return (a <= val && val < b) ? r-l : 0;
    T nv = T(1)<<(MAXLOG-d-1)|val;
    T nnv = ((T(1)<<(MAXLOG-d-1))-1)|nv;
    if(nnv < a || b <= val) return 0;
    if(a <= val && nnv < b) return r - l;
    int lc = mat[d].rank(1, l), rc = mat[d].rank(1, r);
    return freq_dfs(d+1, l-lc, r-rc, val, a, b) + freq_dfs(d+1, lc+zs[d], rc+zs[d], nv, a, b);
  }
  
  // return number of points in [left, right) * [lower, upper)
  int range_freq(int left, int right, T lower, T upper){ return freq_dfs(0, left, right, 0, lower, upper);}
  
  // Sum of the first k+1 intervals [l, r) in ascending order
  T kthLowerSum(int l, int r, int k){
    assert(r - l > k);
    assert(l < r);
    long long kth = 0, res = 0;
    for(int dep = 0; dep < MAXLOG; dep++){
      long long p = mat[dep].rank(0, l);
      long long q = mat[dep].rank(0, r);
      long long bit = (k < q - p) ? 0 : 1;
      if(bit){
        res += acc[dep][r] - acc[dep][l];
        k -= (q - p);
        l += zs[dep] - p;
        r += zs[dep] - q;
      }
      else{
        l = p; r = q;
      }
      kth <<= 1;
      kth |= bit;
    }
    res += kth*(k+1);
    return res;
  }
  
  pair<int, int> ll(int l, int r, T v){
    int res = 0;
    for(int dep = 0; dep < MAXLOG; dep++){
      buff1[dep] = l; buff2[dep] = r;
      bool bit = (v>>(MAXLOG-(dep+1)))&1;
      if(bit) res += r - l + mat[dep].rank(bit, l) - mat[dep].rank(bit, r);
      l = mat[dep].rank(bit, l) + zs[dep]*bit;
      r = mat[dep].rank(bit, r) + zs[dep]*bit;
    }
    return make_pair(res, r - l);
  }
  int lt(int l, int r, T v){
    auto p = ll(l, r, v);
    return p.first;
  }
  int le(int l, int r, T v){
    auto p = ll(l, r, v);
    return p.first + p.second;
  }
  T succ(int l, int r, T v){
    int k = le(l, r, v);
    return k == r - l ? npos : rquantile(l, r, k);
  }
  T pred(int l, int r, T v){
    int k = lt(l, r, v);
    return k ? rquantile(l, r, k-1) : npos;
  }
};


//Warshall_floyd
bool negative_loop = false;
template<typename T>
vector<vector<long long>> warshall(vector<vector<edge<T>>> &g){
  const long long LINF = 1001002003004005006ll;
  int n = (int)g.size();
  vector<vector<long long>> d(n, vector<long long> (n, LINF));
  for(int i = 0; i < n; i++) d[i][i] = 0;
  for(int i = 0; i < n; i++) for(edge<T> e : g[i]) d[i][e.to] = e.cost;
  for(int k = 0; k < n; k++) for(int i = 0; i < n; i++) for(int j = 0; j < n; j++){
    if(d[i][k] != LINF && d[k][j] != LINF){
      d[i][j] = min(d[i][j], d[i][k] + d[k][j]);
    }
  }
  for(int i = 0; i < n; i++) if(d[i][i] < 0) negative_loop = true;
  return d;
}


//WEIGHT_UNIONFIND
template<typename T>
struct WeightUnionFind {
  vector<int> rs, ps;
  vector<T> ws;
  WeightUnionFind(int n): rs(n, 1), ps(n), ws(n, T(0)){
    iota(ps.begin(), ps.end(), 0);
  }
  int find(int x){
    if(x == ps[x]) return x;
    int t = find(ps[x]);
    ws[x] += ws[ps[x]];
    return ps[x] = t;
  }
  T weight(int x){ find(x); return ws[x];}
  bool same(int x, int y){ return find(x) == find(y);}
  void unite(int x, int y, T w){
    w += weight(x);
    w -= weight(y);
    x = find(x); y = find(y);
    if(x == y) return;
    if(rs[x] < rs[y]) swap(x, y), w = -w;
    rs[x] += rs[y];
    ps[y] = x;
    ws[y] = w;
  }
  T diff(int x, int y){ return weight(y) - weight(x);}
};


//Zobrist Hash
template <typename T>
vector<uint64_t> hashing(T a) {
  int n = (int)a.size();
  vector<uint64_t> hash(n+1);
  set<long long> st;
  struct custom_hash {
    static uint64_t splitmix64(uint64_t x) {
      x += 0x9e3779b97f4a7c15;
      x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
      x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
      return x ^ (x >> 31);
    }
    size_t operator() (uint64_t x) const {
      static const uint64_t FIXED_RANDOM =
          chrono::steady_clock::now().time_since_epoch().count();
      return splitmix64(x + FIXED_RANDOM);
    }
  }rng;
  for (int i = 0; i < n; i++) {
    if (st.count(a[i])) { hash[i+1] = hash[i];}
    else { st.insert(a[i]); hash[i+1] = hash[i]^rng(a[i]);}
  }
  return hash;
}


//Z-algorithm
struct Z {
  int n;
  string s;
  vector<int> z;
  Z() {}
  Z(const string& s): s(s) { init();}
  void init() {
    n = (int)s.size();
    z = vector<int>(n);
    z[0] = n;
    for (int i = 1, j = 0; i < n;) {
      while (i+j < n && s[i+j] == s[j]) j++;
      z[i] = j;
      if (j) {
        int k = 1;
        while (i+k < n && k+z[k] < j) {
          z[i+k] = z[k];
          k++;
        }
        i += k; j -= k;
      } else i++;
    }
  }
  int operator[](int i) const { return z[i];}
};





//CLOCK
clock_t begin, end;
begin = clock();
end = clock();
if(double(end - begin)/CLOCKS_PER_SEC > LIMIT){}

//snuke template
#include <bits/stdc++.h>
#define fi first
#define se second
#define rep(i,n) for(int i = 0; i < (n); ++i)
#define rrep(i,n) for(int i = 1; i <= (n); ++i)
#define drep(i,n) for(int i = (n)-1; i >= 0; --i)
#define srep(i,s,t) for (int i = s; i < t; ++i)
#define rng(a) a.begin(),a.end()
#define rrng(a) a.rbegin(),a.rend()
#define isin(x,l,r) ((l) <= (x) && (x) < (r))
#define pb push_back
#define eb emplace_back
#define sz(x) (int)(x).size()
#define pcnt __builtin_popcountll
#define uni(x) x.erase(unique(rng(x)),x.end())
#define snuke srand((unsigned)clock()+(unsigned)time(NULL));
#define show(x) cerr<<#x<<" = "<<x<<endl;
#define PQ(T) priority_queue<T,v(T),greater<T> >
#define bn(x) ((1<<x)-1)
#define dup(x,y) (((x)+(y)-1)/(y))
#define newline puts("")
#define v(T) vector<T>
#define vv(T) v(v(T))
using namespace std;
typedef long long int ll;
typedef unsigned uint;
typedef unsigned long long ull;
typedef pair<int,int> P;
typedef tuple<int,int,int> T;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<ll> vl;
typedef vector<P> vp;
typedef vector<T> vt;
int getInt(){int x;scanf("%d",&x);return x;}
template<typename T>istream& operator>>(istream&i,v(T)&v){rep(j,sz(v))i>>v[j];return i;}
template<typename T>string join(const v(T)&v){stringstream s;rep(i,sz(v))s<<' '<<v[i];return s.str().substr(1);}
template<typename T>ostream& operator<<(ostream&o,const v(T)&v){if(sz(v))o<<join(v);return o;}
template<typename T1,typename T2>istream& operator>>(istream&i,pair<T1,T2>&v){return i>>v.fi>>v.se;}
template<typename T1,typename T2>ostream& operator<<(ostream&o,const pair<T1,T2>&v){return o<<v.fi<<","<<v.se;}
template<typename T>bool mins(T& x,const T&y){if(x>y){x=y;return true;}else return false;}
template<typename T>bool maxs(T& x,const T&y){if(x<y){x=y;return true;}else return false;}
template<typename T>ll suma(const v(T)&a){ll res(0);for(auto&&x:a)res+=x;return res;}
const double eps = 1e-10;
const ll LINF = 1001002003004005006ll;
const int INF = 1001001001;
#define dame { puts("-1"); return 0;}
#define yn {puts("Yes");}else{puts("No");}

#include <iostream> // cout, endl, cin
#include <string> // string, to_string, stoi
#include <vector> // vector
#include <algorithm> // min, max, swap, sort, reverse, lower_bound, upper_bound
#include <utility> // pair, make_pair
#include <tuple> // tuple, make_tuple
#include <cstdint> // int64_t, int*_t
#include <cstdio> // printf
#include <map> // map
#include <queue> // queue, priority_queue
#include <set> // set
#include <stack> // stack
#include <deque> // deque
#include <unordered_map> // unordered_map
#include <unordered_set> // unordered_set
#include <bitset> // bitset
#include <cctype> // isupper, islower, isdigit, toupper, tolower
#include <climits>
#include <iomanip> // setprecision
#include <complex> // complex
#include <math.h>
#include <functional>
#include <cassert>
#define rep(i, a, n) for(int i = a; i < (n); i++)
using namespace std;
using ll = long long;
using P = pair<ll, ll>;
const int INF = 1001001001;
const ll LINF = 1001002003004005006ll;
const int mod = 1000000007;
//const int mod = 998244353;

#include <bits/stdc++.h>
#define rep(i, a, n) for(int i = a; i < (n); i++)
#define drep(i, a, n) for(int i = (n)-1; i >= a; i--)
using namespace std;
using ll = long long;
using P = pair<int, int>;
const int INF = 1001001001;
const ll LINF = 1001002003004005006ll;
//const int mod = 1000000007;
//const int mod = 998244353;

