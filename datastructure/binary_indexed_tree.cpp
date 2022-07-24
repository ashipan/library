//Binary Indexed Tree (Fenwick Tree)
template<typename T>
struct BIT{
    int _n;
    vector<T> data;
    BIT(): _n(0) {};
    explicit BIT(int n): _n(n), data(n) {};
    void add(int p, T x){
        assert(0 <= p && p < _n);
        p++;
        while(p <= _n){
            data[p-1] += x;
            p += p&-p;
        }
    }
    T sum(int r){
        T s = 0;
        while(r > 0){
            s += data[r-1];
            r -= r&-r;
        }
        return s;
    }
    T sum(int l, int r){
        assert(0 <= l && l <= r && r <= _n);
        return sum(r) - sum(l);
    }
};