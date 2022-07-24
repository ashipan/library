// GraphSCC
template<typename T=int>
struct GraphSCC{
public:
    const int n;
    vector<bool> used;
    vector<int> vs, cmp;
    vector<vector<int>> g, rg, h;
    GraphSCC(vector<vector<edge<T>>> &s): n((int)s.size()), used(n), cmp(n), g(n), rg(n) {
        for(int i = 0; i < n; i++) for(int j = 0; j < s[i].size(); j++){
            g[i].emplace_back(s[i][j].to);
            rg[s[i][j].to].emplace_back(i);
        }
    }
    void SCC_dfs_one(int cur){
        used[cur] = true;
        for(int i : g[cur]) if(!used[i]) SCC_dfs_one(i);
        vs.emplace_back(cur);
    }
    void SCC_dfs_two(vector<int> &vec, int cur, int k){
        cmp[cur] = k;
        used[cur] = true;
        vec.push_back(cur);
        for(int i = 0; i < rg[cur].size(); i++){
            if(!used[rg[cur][i]]) SCC_dfs_two(vec, rg[cur][i], k);
        }
    }
    pair<vector<int>, int> scc(){
        for(int i = 0; i < n; i++) if(!used[i]) SCC_dfs_one(i);
        fill(used.begin(), used.end(), false);
        reverse(vs.begin(), vs.end());
        int k = 0;
        vector<vector<int>> s;
        for(int i = 0; i < vs.size(); i++){
            if(!used[vs[i]]){
                s.push_back(vector<int>());
                SCC_dfs_two(s.back(), vs[i], k++);
            }
        }
        h.resize(k);
        fill(used.begin(), used.end(), false);
        for(int i = 0; i < k; i++){
            for(int j = 0; j < s[i].size(); j++){
                int v = s[i][j];
                for(int x = 0; x < g[v].size(); x++){
                    int u = g[v][x];
                    if(used[cmp[u]] || cmp[v] == cmp[u]) continue;
                    used[cmp[u]] = true;
                    h[cmp[v]].push_back(cmp[u]);
                }
            }
            for(int j = 0; j < h[i].size(); j++) used[h[i][j]] = false;
        }
        return make_pair(cmp, k);
    }
};
