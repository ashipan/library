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
#include <iomanip> // setprecision
#include <complex> // complex
#include <math.h>
#include <functional>
#include <cassert>
#define rep(i, a, n) for(int i = a; i < (n); i++)
using namespace std;
using ll = long long;
using P = pair<int, int>;
const int INF = 1001001001;
const ll LINF = 1001002003004005006ll;
const int mod = 1000000007;
//const int mod = 998244353;

struct edge{
    int to, cost;
    edge(int to=0, int cost=0): to(to), cost(cost) {}
};

vector<vector<edge>> to;
vector<int> dist;


int main()
{
    int n, m;
    cin >> n >> m;
    to = vector<vector<edge>> (n);
    dist = vector<int> (n, INF);
    rep(i, 0, m){
        int a, b, c;
        cin >> a >> b >> c;
        a--; b--;
        to[a].emplace_back(b, c);
    }
    priority_queue<P, vector<P>, greater<P>> q;
    rep(i, 0, n){
        fill(dist.begin(), dist.end(), INF);
        for(edge u : to[i]) q.push(P(u.cost, u.to));
        while(!q.empty()){
            int v = q.top().second;
            int d = q.top().first; q.pop();
            if(dist[v] <= d) continue;
            dist[v] = d;
            for(edge u : to[v]){
                int nv = u.to, nd = u.cost;
                q.push(P(d + nd, nv));
            }
        }
        if(dist[i] == INF) cout << -1 << endl;
        else cout << dist[i] << endl;
    }
    return 0;
}

