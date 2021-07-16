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
using P = pair<ll, ll>;
const int INF = 1001001001;
const ll LINF = 1001002003004005006ll;
const int mod = 1000000007;
//const int mod = 998244353;

vector<vector<int>> to, from;
vector<int> tmp;
vector<bool> used;
ll cnt = 0;

void dfs(int v){
    used[v] = true;
    for(int u : to[v]){
        if(used[u]) continue;
        dfs(u);
    }
    tmp.push_back(v);
}

void dfs2(int v){
    used[v] = true;
    for(int u : from[v]){
        if(used[u]) continue;
        dfs2(u);
    }
    cnt++;
}

int main()
{
    int n, m;
    cin >> n >> m;
    to.resize(n);
    used.resize(n);
    from.resize(n);
    rep(i, 0, m){
        int a, b;
        cin >> a >> b;
        a--; b--;
        to[a].push_back(b);
        from[b].push_back(a);
    }
    rep(i, 0, n){
        if(used[i]) continue;
        dfs(i);
    }
    ll ans = 0;
    reverse(tmp.begin(), tmp.end());
    fill(used.begin(), used.end(), false);
    for(int i : tmp){
        if(used[i]) continue;
        cnt = 0;
        dfs2(i);
        ans += cnt*(cnt - 1)/2;
    }
    cout << ans << endl;
    return 0;
}
