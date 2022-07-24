// SegmentTree
using S = ll;
ll op(ll a, ll b){ return min(a, b);}
ll e(){ return ll(1e9);}
int target = (int)(1e9);
bool f(int v){ return v < target;}

template<class S, S (*op)(S, S), S (*e)()>
struct SEG{
    int n, size, log = 0;
    vector<S> d;
    SEG(): SEG(0) {}
    explicit SEG(int n): SEG(vector<S>(n, e())) {}
    explicit SEG(const vector<S>& v): n((int)v.size()){
        while((1<<log) < n) log++;
        size = 1<<log;
        d = vector<S>(2*size, e());
        for(int i = 0; i < n; i++) d[size+i] = v[i];
        for(int i = size-1; i >= 1; i--) update(i);
    }
    void update(int k) { d[k] = op(d[2*k], d[2*k+1]);}
    void set(int p, S x){
        assert(0 <= p && p < n);
        p += size;
        d[p] = x;
        for(int i = 1; i <= log; i++) update(p>>i);
    }
    S get(int p) const{
        assert(0 <= p && p < n);
        return d[p+size];
    }
    S prod(int l, int r) const{
        assert(0 <= l && l <= r && r <= n);
        S sml = e(), smr = e();
        l += size;
        r += size;
        while(l < r){
            if(l&1) sml = op(sml, d[l++]);
            if(r&1) smr = op(d[--r], smr);
            l >>= 1; r >>= 1;
        }
        return op(sml, smr);
    }
    S all_prod() const{ return d[1];}
    template<bool (*f)(S)>
    int max_right(int l) const{ return max_right(l, [](S x) {return f(x);});}
    template<class F>
    int max_right(int l, F f) const{
        assert(0 <= l && l <= n);
        assert(f(e()));
        if(l == n) return n;
        l += size;
        S sm = e();
        do{
            while((l&1) == 0) l >>= 1;
            if(!f(op(sm, d[l]))){
                while(l < size){
                    l = 2*l;
                    if(f(op(sm, d[l]))){
                        sm = op(sm, d[l]);
                        l++;
                    }
                }
                return l - size;
            }
            sm = op(sm, d[l]);
            l++;
        }while((l&-l) != l);
        return n;
    }
    template<bool (*f)(S)>
    int min_left(int r) const{ return min_left(r, [](S x){ return f(x);});}
    template<class F>
    int min_left(int r, F f) const{
        assert(0 <= r && r <= n);
        assert(f(e()));
        if(r == 0) return 0;
        r += size;
        S sm = e();
        do{
            r--;
            while(r > 1 && (r&1)) r >>= 1;
            if(!f(op(d[r], sm))){
                while(r < size){
                    r = 2*r + 1;
                    if(f(op(d[r], sm))){
                        sm = op(d[r], sm);
                        r--;
                    }
                }
                return r + 1 - size;
            }
            sm = op(d[r], sm);
        }while((r&-r) != r);
        return 0;
    }
};