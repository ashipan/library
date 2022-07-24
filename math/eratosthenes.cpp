//Eratosthenes
template<typename T>
struct Eratosthenes{
    vector<bool> isprime;
    vector<T> sieves;
    vector<T> minfactor;
    vector<T> mobius;
    Eratosthenes(T n=0):isprime(n+1, true), minfactor(n+1, -1), mobius(n+1, 1){
        isprime[1] = false;
        minfactor[1] = 1;
        for(T i = 2; i <= n; i++){
            if(!isprime[i]) continue;
            minfactor[i] = i;
            mobius[i] = -1;
            for(T j = i*2; j <= n; j += i){
                isprime[j] = false;
                if(minfactor[j] == -1) minfactor[j] = i;
                if((j/i)%i) mobius[j] = -mobius[j];
                else mobius[j] = 0;
            }
        }
        for(T i = 2; i <= n; i++) if(isprime[i]) sieves.emplace_back(i);
    }
    vector<pair<T, T>> factorize(T n){
        vector<pair<T, T>> res;
        while(n > 1){
            int p = minfactor[n];
            int exp = 0;
            while(minfactor[n] == p){
                n /= p;
                exp++;
            }
            res.emplace_back(p, exp);
        }
        return res;
    }
    vector<T> divisors(T n){
        vector<T> res({1});
        auto pf = factorize(n);
        for(auto p : pf){
            int s = (int)res.size();
            for(int i = 0; i < s; i++){
                T v = 1;
                for(int j = 0; j < p.second; j++){
                    v *= p.first;
                    res.push_back(res[i]*v);
                }
            }
        }
        return res;
    }
};