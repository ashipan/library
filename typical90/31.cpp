#include <bits/stdc++.h>
#define rep(i, a, n) for(int i = a; i < (n); i++)
using namespace std;
using ll = long long;
using P = pair<ll, ll>;
const int INF = 1001001001;
const ll LINF = 1001002003004005006ll;
//const int mod = 1000000007;
//const int mod = 998244353;

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

int main()
{
  int n;
  cin >> n;
  vector<int> a(n), b(n);
  rep(i, 0, n) cin >> a[i];
  rep(i, 0, n) cin >> b[i];
  vector<vector<int>> dp(51, vector<int> (1505, -1));
  auto f = [&](auto f, int w, int b) -> int {
    if (dp[w][b] != -1) return dp[w][b];
    MEX<int> mex;
    if (w == 0 && b <= 1) return 0;
    if (w >= 1) mex.insert(f(f, w-1, w+b));
    rep(k, 1, b/2+1) mex.insert(f(f, w, b-k));
    return dp[w][b] = mex.mex();
  };
  ll flag = 0;
  rep(i, 0, n) flag ^= f(f, a[i], b[i]);
  if (flag) cout << "First" << endl;
  else cout << "Second" << endl;
  return 0;
}

