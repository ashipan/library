//edge
template<typename T=int>
struct edge{
    int from, to, id; T cost;
    edge(int to=0, T cost=1): to(to), cost(cost) {}
    //edge(int to=0, int id=0): to(to), id(id) {}
    //edge(int from=0, int to=0, T cost=1): from(from), to(to), cost(cost) {}
    //edge(int to=0, int from=0, int id=0): to(to), from(from), id(id) {}
    bool operator<(const edge &e) const{return cost < e.cost;}
    bool operator>(const edge &e) const{return cost > e.cost;}
};
template<typename T=int>
using graph = vector<vector<edge<T>>>;