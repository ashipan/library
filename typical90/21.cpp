#include <bits/stdc++.h>
#define rep(i, a, n) for(int i = a; i < (n); i++)
using namespace std;
using ll = long long;
using P = pair<ll, ll>;
const int INF = 1001001001;
const ll LINF = 1001002003004005006ll;
//const int mod = 1000000007;
//const int mod = 998244353;

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

int main()
{
  int n, m;
  cin >> n >> m;
  graph<int> to(n);
  rep(i, 0, m) {
    int a, b;
    cin >> a >> b;
    a--; b--;
    to[a].emplace_back(b);
  }
  scc_graph<int> g(to);
  ll ans = 0;
  for (auto s : g.scc()) {
    if (s.size() >= 2) ans += (s.size())*(s.size()-1)/2;
  }
  cout << ans << endl;
  return 0;
}
