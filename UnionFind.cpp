//UNIONFIND
template<typename T>
struct UnionFind {
    int size_;
    vector<T> d;
    UnionFind(int n=0): d(n,-1), size_(n) {}
    T find(int x) {
        if (d[x] < 0) return x;
        return d[x] = find(d[x]);
    }
    bool unite(int x, int y) {
        x = find(x); y = find(y);
        if (x == y) return false;
        size_--;
        if (d[x] > d[y]) swap(x,y);
        d[x] += d[y];
        d[y] = x;
        return true;
    }
    bool same(int x, int y) { return find(x) == find(y);}
    T size(int x) { return -d[find(x)];}
    int size() { return size_;}
};
