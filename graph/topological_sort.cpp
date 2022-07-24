// topological sort - Tarjan.ver
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

// topological sort - Kahn.ver (Closed circuit detectable)
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