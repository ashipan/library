// prim
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