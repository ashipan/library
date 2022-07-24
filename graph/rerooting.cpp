// Rerooting
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
