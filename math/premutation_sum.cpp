// PERMUTATION_SUM
ll permutation_sum(ll l, ll r, ll n){ return (l+r)*n/2;}
mint mod_permutation_sum(ll l, ll r, ll n){ return mint(l+r)*n/2;}
ll permutation_sum2(ll a, ll d, ll n){ return (2*a + (n-1)*d)*n/2;}
mint mod_permutation2(ll a, ll d, ll n){ return mint(2*a + (n-1)*d)*n/2;}
ll permutation_sum3(ll l, ll r, ll d=1){ return (l+r)*(r-l+1)/2;}
mint mod_permutation_sum3(ll l, ll r, ll d=1){ return mint(l+r)*(r-l+1)/2;}