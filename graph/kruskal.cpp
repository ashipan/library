//Kruskal
template<typename T>
T kruskal(vector<vector<edge<T>>>& g){
    int n = (int)g.size();
    vector<edge<T>> es;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < g[i].size(); j++){
            es.push_back(g[i][j]);
        }
    }
    sort(es.begin(), es.end());
    UnionFind<int> uf(n);
    T res = 0;
    for(int i = 0; i < es.size(); i++){
        edge<T> e = es[i];
        if(!uf.same(e.from, e.to)){
            uf.unite(e.from, e.to);
            res += e.cost;
        }
    }
    return res;
}