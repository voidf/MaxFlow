#include <memory>
#include <cmath>
#include <vector>
#include <queue>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <numeric>
#include <algorithm>
using namespace std;
#define vec vector
#define mp unordered_map

template <typename T, typename U>
struct cached_edge
{
    cached_edge(T to, U c, T rev) : to(to), rev(rev), cap(c), rcap(c) {}
    T to, rev;
    U cap, rcap;
};

template <class T, class U>
struct Dinic
{
    vec<T> cur;
    // vec<T> nxt;
    vec<T> level;
    // vec<T> to;
    // vec<U> flow;
    vec<vec<cached_edge<T, U>>> &rnet;

    T _src, _dst;
    const U INF = ~0u >> 1;

    inline void prepare(const T s, const T t)
    {
        _src = s;
        _dst = t;
        for (auto &u : rnet)
            for (auto &edge : u)
                edge.cap = edge.rcap = (edge.cap + edge.rcap) / 2;
    }

    Dinic(vec<vec<cached_edge<T, U>>> &graph) : rnet(graph) { level.assign(rnet.size(), 0); }

    inline bool bfs(const T s, const T t)
    {
        level.assign(rnet.size(), 0);
        level[t] = 1;
        queue<T> Q({t});
        while (Q.size() && !level[s])
        {
            int p = Q.front();
            Q.pop();
            for (const auto &e : rnet[p])
                if (e.rcap && !level[e.to])
                {
                    level[e.to] = level[p] + 1;
                    Q.push(e.to);
                }
        }
        return !!level[s];
    }

    U dfs(const T p, const U maxflow)
    {
        if (p == _dst)
            return maxflow;
        U nowflow = maxflow, tflow;
        for (T &cid = cur[p]; cid < rnet[p].size(); ++cid)
        {
            auto &e = rnet[p][cid];
            if (e.cap && level[p] == level[e.to] + 1)
            {
                tflow = dfs(e.to, min(e.cap, nowflow));
                if (tflow)
                {
                    rnet[e.to][e.rev].rcap -= tflow;
                    e.cap -= tflow;
                    rnet[e.to][e.rev].cap += tflow;
                    e.rcap += tflow;
                    nowflow -= tflow; // nowflow实际上就是本点的超额流
                }
                if (!nowflow)
                    break;
            }
        }
        if (nowflow == maxflow)
            level[p] = 0;
        return maxflow - nowflow;
    }

    U dinic()
    {
        U maxflow = 0, tflow;
        while (bfs(_src, _dst))
        {
            cur.assign(rnet.size(), 0);
            while (tflow = dfs(_src, INF))
                maxflow += tflow;
        }
        return maxflow;
    }
    // ===
    auto stcut() const
    {
        vec<int> D(rnet.size());
        vec<int> Q(rnet.size());
        Q[0] = _src;
        D[_dst] = 1;
        int top = 0;
        int cap = 1;
        while (top < cap)
        {
            int v = Q[top++];
            for (auto &e : rnet[v])
            {
                if (!D[e.to] && e.cap > 0)
                {
                    D[e.to] = D[v] + 1;
                    Q[cap++] = e.to;
                }
            }
        }
        return D;
    }

    auto cutset(T *d) const
    {

        vec<vec<T>> cset;
        for (T u = 0; u < rnet.size(); ++u)
            for (auto &edge : rnet[u])
                if (u > edge.to && (d[u] == 0) != (0 == d[edge.to]))
                    cset.push_back({u, edge.to});
        return cset;
    }
};

template <class U>
struct Solver
{
    vec<vec<cached_edge<int, U>>> G;

    void add_edge(int u, int v, U w)
    {
        G[u].emplace_back(v, w, G[v].size());
        G[v].emplace_back(u, w, G[u].size() - 1);
    }

    unique_ptr<int[]> node, tmp1, d;
    vec<vec<U>> ans;
    vec<vec<int>> ansedge;
    vec<vec<vec<int>>> cutsets;
    int _n;
    Solver(int n) : node(make_unique<int[]>(n)), tmp1(make_unique<int[]>(n)), d(make_unique<int[]>(n)),
                    //    ans(n, vec<int>(n, ~0u >> 1)),
                    _n(n),
                    G(n)
    {
        init();
    }
    void init()
    {
        iota(node.get(), node.get() + _n, 0);
        ans.assign(_n, vec<U>(_n, ~0u >> 1));
        ansedge.assign(_n, vec<int>(_n, ~0u >> 1));
    }

    void work(int l, int r, Dinic<int, U> &din)
    {
        if (l == r)
            return;
        int S = node[l], T = node[l + 1];
        U t = 0;

        fill(d.get(), d.get() + _n, 0);
        int c = 0;
        din.prepare(S, T);
        t = din.dinic();
        for (auto i : din.level)
            d[c++] = !i;
        ans[T][S] = ans[S][T] = t;
        int pos = ansedge[T][S] = ansedge[S][T] = cutsets.size();
        // cerr << "DEBUG: [before] cutsets size:" << cutsets.size() << endl;
        cutsets.emplace_back(din.cutset(d.get()));
        // cerr << "flow: " << t << ", pos:" << pos << ", size:" << cutsets.back().size() << endl;
        // cerr << "DEBUG: [after] cutsets size:" << cutsets.size() << endl;

        int cnt1 = 0, cnt2 = r - l;
        for (int i = l; i <= r; ++i)
            if (d[node[i]])
                tmp1[cnt1++] = node[i];
            else
                tmp1[cnt2--] = node[i];
        copy(tmp1.get(), tmp1.get() + r - l + 1, node.get() + l);
        cnt2 = r - l - cnt2;
        work(l, l + cnt1 - 1, din);
        work(l + cnt1, r, din);
        for (int i = 0; i < cnt1; ++i)
            for (int j = 0; j < cnt2; ++j)
            {
                int ii = node[i + l], jj = node[j + cnt1 + l];
                U minans = ans[ii][S];
                int mincset = ansedge[ii][S];
                if (minans > ans[S][T])
                {
                    minans = ans[S][T];
                    mincset = pos;
                }
                if (minans > ans[T][jj])
                {
                    minans = ans[T][jj];
                    mincset = ansedge[T][jj];
                }
                ans[jj][ii] = ans[ii][jj] = minans;
                ansedge[jj][ii] = ansedge[ii][jj] = mincset;
                // cerr << "[" << jj << ", " << ii << "] minans:" << minans << ", mincset:" << mincset << ", size:" << cutsets[mincset].size() << endl;
            }
    }

    void solve()
    {
        Dinic<int, U> din(G);
        work(0, _n - 1, din);
    }
};

signed main()
{
    using U = int;
    ios::sync_with_stdio(0);
    int n, m;
    cin >> n >> m;
    ++n;
    auto sol = Solver<U>(n);
    while (m--)
    {
        int u, v;
        U w;
        cin >> u >> v >> w;
        sol.add_edge(u, v, w);
    }
    sol.solve();
    cin >> n;
    while (n--)
    {
        int x, y;
        cin >> x >> y;
        cout << sol.ans[x][y] << '\n';
    }

    return 0;
}