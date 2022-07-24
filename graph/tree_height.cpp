// TreeHeight
template<typename T>
vector<T> treeHeight(vector<vector<edge<T>>> &g){
    int n = (int)g.size();
    vector<T> v1, v2, v3, ret(n);
    int p1, p2;
    v1 = dijkstra<T>(g, 0);
    p1 = max_element(v1.begin(), v1.end()) - v1.begin();
    v2 = dijkstra<T>(g, p1);
    p2 = max_element(v2.begin(), v2.end()) - v2.begin();
    v3 = dijkstra<T>(g, p2);
    for(int i=0; i<n; i++) ret[i] = max(v2[i], v3[i]);
    return ret;
}