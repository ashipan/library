//Chinese Remainder Theorem
// x ≡ r (mod m)
// r = m1*x1+b1
// r = m2*x2+b2
pair<ll, ll> ChineseRem(ll b1, ll m1, ll b2, ll m2){
    ll p, q;
    ll d = extGCD(m1, m2, p, q); // p is inv of m1/d (mod m2/d)
    if((b2 - b1)%d != 0) return make_pair(0, -1); // no answer
    ll m = m1*(m2/d); // lcm(m1, m2)
    ll tmp = (b2 - b1)/d*p%(m2/d);
    ll r = mod(b1 + m1*tmp, m);
    return make_pair(r, m);
}