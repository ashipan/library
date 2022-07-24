//BFS
vector<vector<int>> dist(h, vector<int> (w, INF));
queue<P> q;
auto push = [&](int i, int j, int x){
    if(dist[i][j] != INF) return;
    dist[i][j] = x;
    q.push(P(i, j));
};
auto pop = [&](){
    while(!q.empty()){
        int i = q.front().first;
        int j = q.front().second; q.pop();
        rep(dir, 0, 4){
            int ni = i + di[dir], nj = j + dj[dir];
            if(ni < 0 || ni >= h || nj < 0 || nj >= w) continue;
            if(s[ni][nj] == '#') continue;
        }
    }
};
