//PRIME_FACTOR
vector<pair<ll, ll>> prime_factor(ll n){
    vector<pair<ll, ll>> pf;
    //if(n == 1) pf.push_back(P(1, 1));
    for(ll i = 2; i*i <= n; i++){
        if(n%i == 0){
            int cnt = 0;
            while(n%i == 0){
                n /= i;
                cnt++;
            }
            pf.push_back(pair<ll, ll>(i, cnt));
        }
    }
    if(n != 1) pf.push_back(pair<ll, ll>(n, 1));
    return pf;
}