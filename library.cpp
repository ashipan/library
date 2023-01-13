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
const int di[] = {-1, 0, 1, 0, -1, -1, 1, 1};
const int dj[] = {0, -1, 0, 1, -1, 1, -1, 1};

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
// x ≡ r (mod m)
// r = m1*x1+b1
// r = m2*x2+b2
pair<ll, ll> ChineseRem(ll b1, ll m1, ll b2, ll m2){
  ll p, q;
  ll d = extGCD(m1, m2, p, q); // p is inv of m1/d (mod m2/d)
  if((b2 - b1)%d != 0) return make_pair(0, -1); // no answer
  ll m = m1*(m2/d); // lcm(m1, m2)
  ll tmp = (b2 - b1)/d*p%(m2/d);
  ll r = mod(b1 + m1*tmp, m);
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
  modinv(): n(2), d({0, 1}) {};
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
vector<T> dijkstra(const vector<vector<edge<T>>> &g, int s){
  int n = (int)g.size();
  vector<T> dist(n, numeric_limits<T>::max());
  priority_queue<pair<T, int>, vector<pair<T, int>>, greater<pair<T, int>>> q;
  q.emplace(0, s);
  while(!q.empty()){
    int v = q.top().second;
    int d = q.top().first; q.pop();
    if(dist[v] <= d) continue;
    dist[v] = d;
    for(edge e : g[v]) q.emplace(e.cost+d, e.to);
  }
  return dist;
}


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


//GCD,LCM
ll gcd(ll a, ll b) { return b ? gcd(b, a%b) : a;}
ll lcm(ll a, ll b) { return a/gcd(a, b)*b;}


//GraphSCC
template<typename T=int>
struct GraphSCC{
public:
  const int n;
  vector<bool> used;
  vector<int> vs, cmp;
  vector<vector<int>> g, rg, h;
  GraphSCC(vector<vector<edge<T>>> &s): n((int)s.size()), used(n), cmp(n), g(n), rg(n) {
    for(int i = 0; i < n; i++) for(int j = 0; j < s[i].size(); j++){
      g[i].emplace_back(s[i][j].to);
      rg[s[i][j].to].emplace_back(i);
    }
  }
  void SCC_dfs_one(int cur){
    used[cur] = true;
    for(int i : g[cur]) if(!used[i]) SCC_dfs_one(i);
    vs.emplace_back(cur);
  }
  void SCC_dfs_two(vector<int> &vec, int cur, int k){
    cmp[cur] = k;
    used[cur] = true;
    vec.push_back(cur);
    for(int i = 0; i < rg[cur].size(); i++){
      if(!used[rg[cur][i]]) SCC_dfs_two(vec, rg[cur][i], k);
    }
  }
  pair<vector<int>, int> scc(){
    for(int i = 0; i < n; i++) if(!used[i]) SCC_dfs_one(i);
    fill(used.begin(), used.end(), false);
    reverse(vs.begin(), vs.end());
    int k = 0;
    vector<vector<int>> s;
    for(int i = 0; i < vs.size(); i++){
      if(!used[vs[i]]){
        s.push_back(vector<int>());
        SCC_dfs_two(s.back(), vs[i], k++);
      }
    }
    h.resize(k);
    fill(used.begin(), used.end(), false);
    for(int i = 0; i < k; i++){
      for(int j = 0; j < s[i].size(); j++){
        int v = s[i][j];
        for(int x = 0; x < g[v].size(); x++){
          int u = g[v][x];
          if(used[cmp[u]] || cmp[v] == cmp[u]) continue;
          used[cmp[u]] = true;
          h[cmp[v]].push_back(cmp[u]);
        }
      }
      for(int j = 0; j < h[i].size(); j++) used[h[i][j]] = false;
    }
    return make_pair(cmp, k);
  }
};


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
};
mint ex(mint x, ll t) { return x.pow(t);}
istream& operator>>(istream& i, mint& a) { unsigned long long t; i>>t; a=mint(t); return i;}
ostream& operator<<(ostream& o, const mint& a) { return o<<a.x;}


//MOD
inline ll mod(ll a, ll m) { return (a%m+m)%m;}


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


//Lazy Segment Tree(range_max)
template<typename T>
struct SEG{
  int n;
  vector<T> seg, lazy;
  SEG(){}
  SEG(int n_) { init(n_);}
  void init(int n_){
    n = 1;
    while(n < n_) n *= 2;
    seg = lazy = vector<T> (n*2, -INT_MAX);
  }
  void push(int k){
    if(k < n){
      lazy[k*2] = max(lazy[k*2], lazy[k]);
      lazy[k*2+1] = max(lazy[k*2+1], lazy[k]);
    }
    seg[k] = max(seg[k], lazy[k]);
    lazy[k] = 0;
  }
  void update(int a, int b, T x, int k, int l, int r){
    push(k);
    if(r <= a || b <= l) return;
    if(a <= l && r <= b){
      lazy[k] = x;
      push(k);
      return;
    }
    update(a, b, x, k*2, l, (l + r)/2);
    update(a, b, x, k*2 + 1, (l + r)/2, r);
    seg[k] = max(seg[k*2], seg[k*2+1]);
  }
  void update(int l, int r, T x){ update(l, r, x, 1, 0, n);}
  int range_max(int a, int b, T k, int l, int r){
    push(k);
    if(r <= a || b <= l) return 0;
    if(a <= l && r <= b) return seg[k];
    int vl = range_max(a, b, k*2, l, (l + r)/2);
    int vr = range_max(a, b, k*2 + 1, (l + r)/2, r);
    return max(vl, vr);
  }
  int range_max(int l, int r){ return range_max(l, r, 1, 0, n);}
};


//Lazy Segment Tree(range_min)
template<typename T>
struct SEG{
  int n;
  vector<T> seg, lazy;
  SEG(){}
  SEG(int n_) { init(n_);}
  void init(int n_){
    n = 1;
    while(n < n_) n *= 2;
    seg = lazy = vector<T> (n*2, INT_MAX);
  }
  void push(int k){
    if(k < n){
      lazy[k*2] = min(lazy[k*2], lazy[k]);
      lazy[k*2+1] = min(lazy[k*2+1], lazy[k]);
    }
    seg[k] = min(seg[k], lazy[k]);
    lazy[k] = INT_MAX;
  }
  void update(int a, int b, T x, int k, int l, int r){
    push(k);
    if(r <= a || b <= l) return;
    if(a <= l && r <= b){
      lazy[k] = x;
      push(k);
      return;
    }
    update(a, b, x, k*2, l, (l + r)/2);
    update(a, b, x, k*2 + 1, (l + r)/2, r);
    seg[k] = min(seg[k*2], seg[k*2+1]);
  }
  void update(int l, int r, T x){ update(l, r, x, 1, 0, n);}
  int range_min(int a, int b, T k, int l, int r){
    push(k);
    if(r <= a || b <= l) return 0;
    if(a <= l && r <= b) return seg[k];
    int vl = range_min(a, b, k*2, l, (l + r)/2);
    int vr = range_min(a, b, k*2 + 1, (l + r)/2, r);
    return min(vl, vr);
  }
  int range_min(int l, int r){ return range_min(l, r, 1, 0, n);}
};


// Lazy_segtree

/* RMQ, RAQ
 using S = ll;
 using F = ll;
 S op(S a, S b){ return max(a, b);}
 S e(){ return -INF;}
 S mapping(F f, S x){ return f+x;}
 F composition(F f, F g){ return f+g;}
 F id(){ return 0;}
 */

/* RSQ, RAQ
 struct S{
 ll val;
 int size;
 };
 using F = ll;
 S op(S a, S b){ return {a.val + b.val, a.size + b.size};}
 S e(){ return {0, 1};}
 S mapping(F f, S x){ return {x.val + f*x.size, x.size};}
 F composition(F f, F g){ return f + g;}
 F id() { return 0;}
 */

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

//LIS
template<typename T>
vector<T> LIS(vector<T>& a){
  int n = (int)a.size();
  vector<T> dp(n+1, numeric_limits<T>::max());
  for(int i = 0; i < n; i++) *lower_bound(dp.begin(), dp.end(), a[i]) = a[i];
  return dp;
}

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
  //if(n == 1) pf.push_back(P(1, 1));
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
      if(p - q > k){
        l = p + zs[dep];
        r = p + zs[dep];
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
vector<vector<T>> warshall(vector<vector<edge<T>>> &g){
  int n = (int)g.size();
  vector<vector<T>> d(n, vector<T> (n, numeric_limits<T>::max()));
  for(int i = 0; i < n; i++) d[i][i] = 0;
  for(int i = 0; i < n; i++) for(edge<T> e : g[i]) d[i][e.to] = e.cost;
  for(int k = 0; k < n; k++) for(int i = 0; i < n; i++) for(int j = 0; j < n; j++){
    if(d[i][k] != numeric_limits<T>::max() && d[k][j] != numeric_limits<T>::max()){
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

