//Bellman-Ford
template<typename T>
bool find_negative_loop(vector<vector<edge<T>>> &g, int s){
    int n = (int)g.size();
    vector<T> d(n, numeric_limits<T>::max()); d[s] = 0;
    for(int i = 0; i < n; i++) for(int j = 0; j < n; j++) for(edge<T> e : g[j]){
        if(d[j] != numeric_limits<T>::max() && d[e.to] > d[j] + e.cost){
            d[e.to] = d[j] + e.cost;
            if(i == n - 1) return true;
        }
    }
    return false;
}
template<typename T>
bool find_negative_loop_whole(vector<vector<edge<T>>> &g){
    int n = (int)g.size();
    vector<T> d(n);
    for(int i = 0; i < n; i++) for(int j = 0; j < n; j++) for(edge<T> e : g[j]){
        if(d[e.to] > d[j] + e.cost){
            d[e.to] = d[j] + e.cost;
            if(i == n - 1) return true;
        }
    }
    return false;
}
template<typename T>
vector<T> bellman(const vector<vector<edge<T>>> &g, int s){
    int n = (int)g.size();
    vector<T> d(n, numeric_limits<T>::max()); d[s] = 0;
    while(true){
        bool update = false;
        for(int i = 0; i < n; i++) for(edge<T> e : g[i]){
            if(d[i] != numeric_limits<T>::max() && d[e.to] > d[i] + e.cost){
                d[e.to] = d[i] + e.cost;
                update = true;
            }
        }
        if(!update) break;
    }
    return d;
}