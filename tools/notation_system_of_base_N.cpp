//Notation systemp of base N
ll baseN_to_long(string s, int N){
    ll res = 0;
    for(int i = 0; i < s.size(); i++){
        if('A' <= s[i] && s[i] <= 'Z') res = res*N + s[i] - 'A' + 10;
        else res = res*N + s[i] - '0';
    }
    return res;
}

string long_to_baseN(ll n, int N){
    if(n == 0) return "0";
    string s = "0123456789ABCDEF";
    string res;
    while(n > 0){ res = s[n%N] + res; n /= N;}
    return res;
}