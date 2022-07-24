//UNIONFIND
template<typename T>
struct UnionFind {
    int num;
    vector<T> d;
    UnionFind(int n=0): num(n), d(n,-1) {}
    T find(int x) {
        if (d[x] < 0) return x;
        return d[x] = find(d[x]);
    }
    bool unite(int x, int y) {
        x = find(x); y = find(y);
        if (x == y) return false;
        if (d[x] > d[y]) swap(x,y);
        d[x] += d[y];
        d[y] = x;
        num--;
        return true;
    }
    bool same(int x, int y) { return find(x) == find(y);}
    T size(int x) { return -d[find(x)];}
    int count() { return num;}
};