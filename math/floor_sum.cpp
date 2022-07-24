//FLOOR_SUM ∑[i=0, n−1](floor(a×i+b)/m)
template<typename T>
T floor_sum(T n, T m, T a, T b){
    T ans = 0;
    if(a >= m){
        ans += (n - 1)*n*(a/m)/2;
        a %= m;
    }
    if(b >= m){
        ans += n*(b/m);
        b %= m;
    }
    T y_max = (a*n + b)/m, x_max = (y_max*m - b);
    if(y_max == 0) return ans;
    ans += (n - (x_max + a - 1)/a)*y_max;
    ans += floor_sum(y_max, a, m, (a - x_max%a)%a);
    return ans;
}