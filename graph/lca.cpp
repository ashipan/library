// LCA
template<typename T>
struct treeLCA{
    int const MAX_LOG_V;
    vector<vector<edge<T>>> g;
    int root, vn;
    vector<vector<int>> parent;
    vector<int> depth;
    treeLCA(vector<vector<edge<T>>> &_g, int _r=0):
        MAX_LOG_V(35), g(_g), root(_r), vn((int)g.size()),
        parent(MAX_LOG_V, vector<int>(vn, 0)), depth(vn, -1)
        {depth[root] = 0; init(vn);}
    void dfs(int v, int p, int d){
        parent[0][v] = p;
        depth[v] = d;
        for(int i=0; i<g[v].size(); i++){
            if(depth[ g[v][i].to ] >= 0) continue;
            if(g[v][i].to != p) dfs(g[v][i].to, v, d+1);
        }
    }
    void init(int V) {
        dfs(root, -1, 0);
        for(int k=0; k+1 < MAX_LOG_V; k++){
            for(int v=0; v < V; v++) {
                if(parent[k][v] < 0) parent[k+1][v] = -1;
                else parent[k+1][v] = parent[k][parent[k][v]];
            }
        }
    }
    int lca(int u, int v) {
        if(depth[u] > depth[v]) swap(u, v);
        for(int k=0; k < MAX_LOG_V; k++){
            if((depth[v] - depth[u]) >> k & 1){
                v = parent[k][v];
            }
        }
        if(u == v) return u;
        for(int k=MAX_LOG_V - 1; k>=0; k--){
            if(parent[k][u] != parent[k][v]){
                u = parent[k][u];
                v = parent[k][v];
            }
        }
        return parent[0][u];
    }
    int dist(int u, int v){
        int anc = lca(u, v);
        return depth[u] + depth[v] - 2*depth[anc];
    }
};