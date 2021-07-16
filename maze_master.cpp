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

const int di[] = {-1, 0, 1, 0, -1, -1, 1, 1};
const int dj[] = {0, -1, 0, 1, -1, 1, -1, 1};

int main()
{
    int h, w;
    cin >> h >> w;
    vector<string> s(h);
    rep(i, 0, h) cin >> s[i];
    vector<vector<int>> dist;
    queue<pair<int, int>> q;
    auto push = [&](int i, int j, int x){
        dist[i][j] = x;
        q.push(make_pair(i, j));
    };
    auto pop = [&](){
        while(!q.empty()){
            int i = q.front().first;
            int j = q.front().second; q.pop();
            for(int dir = 0; dir < 4; dir++){
                int ni = i + di[dir], nj = j + dj[dir];
                if(ni < 0 || ni >= h || nj < 0 || nj >= w) continue;
                if(s[ni][nj] == '#') continue;
                if(dist[ni][nj] != INF) continue;
                push(ni, nj, dist[i][j] + 1);
            }
        }
    };
    
    int ans = 0;
    rep(i, 0, h){
        rep(j, 0, w){
            if(s[i][j] == '#') continue;
            dist = vector<vector<int>> (h, vector<int> (w, INF));
            push(i, j, 0);
            pop();
            rep(ni, 0, h){
                rep(nj, 0, w){
                    if(dist[ni][nj] == INF) continue;
                    ans = max(ans, dist[ni][nj]);
                }
            }
        }
    }
    cout << ans << endl;
    return 0;
}
