//FRACTION
template<typename T=ll>
struct frac{
    T a, b;
    frac(T _a=0, T _b=1): a(_a), b(_b) {
        if(b == 0) { a = 1; return;}
        if(b < 0) a = -a, b = -b;
        T g = gcd(abs(a), b);
        a /= g; b /= g;
    }
    bool operator<(const frac& x) const{ return a*x.b < x.a*b;}
    bool operator==(const frac& x) const{ return a == x.a && b == x.b;}
};