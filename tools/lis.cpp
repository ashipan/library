// LIS
template<typename T>
vector<T> LIS(vector<T>& a){
    int n = (int)a.size();
    vector<T> dp(n+1, numeric_limits<T>::max());
    for(int i = 0; i < n; i++) *lower_bound(dp.begin(), dp.end(), a[i]) = a[i];
    return dp;
}