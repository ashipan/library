// graphLink
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