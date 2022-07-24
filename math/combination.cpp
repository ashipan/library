//COMBINATION
ll combination(ll n, ll r){
    ll res = 1; r = min(r, n - r);
    rep(i, 0, r){ res *= n - i; res /= i + 1;}
    return res;
}
//COMBINATION_MODver.
mint mod_combination(int n, int r){
    mint res = 1; r = min(r, n - r);
    rep(i, 0, r){ res *= n - i; res /= i + 1;}
    return res;
}
//COMBINATION ENUMERATE
struct modinv{
    int n;
    vector<mint> d;
    modinv(): n(2), d({0, 1}) {};
    mint operator()(int i){ while(n <= i){ d.push_back(-d[mod%n]*(mod/n)); n++;} return d[i];}
    mint operator[](int i) const { return d[i];}
}invs;
struct modfact{
    int n;
    vector<mint> d;
    modfact(): n(2), d({1, 1}) {}
    mint operator()(int i) { while(n <= i){ d.push_back(d.back()*n); n++;} return d[i];}
    mint operator[](int i) const { return d[i];}
}facs;
struct modfactinv{
    int n;
    vector<mint> d;
    modfactinv(): n(2), d({1, 1}) {}
    mint operator()(int i) { while(n <= i){ d.push_back(d.back()*invs(n)); n++;} return d[i];}
    mint operator[](int i) const { return d[i];}
}ifacs;
mint comb(int a, int b){
    if(a < b || b < 0) return 0;
    return facs(a)*ifacs(b)*ifacs(a-b);
}