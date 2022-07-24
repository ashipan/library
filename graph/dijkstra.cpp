//DIJKSTRA
template<typename T>
vector<T> dijkstra(const vector<vector<edge<T>>> &g, int s){
    int n = (int)g.size();
    vector<T> dist(n, numeric_limits<T>::max());
    priority_queue<pair<T, int>, vector<pair<T, int>>, greater<pair<T, int>>> q;
    q.emplace(0, s);
    while(!q.empty()){
        int v = q.top().second;
        int d = q.top().first; q.pop();
        if(dist[v] <= d) continue;
        dist[v] = d;
        for(edge e : g[v]) q.emplace(e.cost+d, e.to);
    }
    return dist;
}