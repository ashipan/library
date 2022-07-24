// Warshall_floyd
bool negative_loop = false;
template<typename T>
vector<vector<T>> warshall(vector<vector<edge<T>>> &g){
    int n = (int)g.size();
    vector<vector<T>> d(n, vector<T> (n, numeric_limits<T>::max()));
    for(int i = 0; i < n; i++) d[i][i] = 0;
    for(int i = 0; i < n; i++) for(edge<T> e : g[i]) d[i][e.to] = e.cost;
    for(int k = 0; k < n; k++) for(int i = 0; i < n; i++) for(int j = 0; j < n; j++){
        if(d[i][k] != numeric_limits<T>::max() && d[k][j] != numeric_limits<T>::max()){
            d[i][j] = min(d[i][j], d[i][k] + d[k][j]);
        }
    }
    for(int i = 0; i < n; i++) if(d[i][i] < 0) negative_loop = true;
    return d;
}