// POWER_MODver. N^k % MOD
ll mod_pow(ll n, ll k){
    ll res = 1;
    for(; k > 0; k >>= 1){
        if(k&1) res = (res*n)%mod;
        n = (n*n)%mod;
    }
    return res;
}
