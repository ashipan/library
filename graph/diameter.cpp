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