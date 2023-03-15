#include <memory>
#include <cmath>
#include <vector>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <numeric>
#include <algorithm>
using namespace std;
#define vec vector
#define mp unordered_map

template <class T>
struct HLPP
{
    const T INF = ~0u >> 1;
    struct edge
    {
        int to, rev;
        T f;
    };
    vec<edge> edges; // edges[MAXM << 1] write on push()

    int maxn, s, t;

    vec<int> first_edge; // [MAXN + 1]  readonly
    vec<int> _cur_edge;  // [MAXN];     write on addEdge()
    vec<int> nxt;        // [MAXN];     write on pushLst()
    vec<int> lst;        // [MAXN];     write on pushLst(), globalRelabel(), calc()
    vec<T> excess;       //[MAXN];      write on push(), calc()
    vec<int> arc;        // [MAXN];     write on discharge(), calc()

    vec<int> gapNxt; //[MAXN << 1];     write on globalRelabel(), updHeight(), discharge()
    vec<int> gapPrv; //[MAXN << 1];     write on globalRelabel(), updHeight(), discharge()

    vec<int> height; //[MAXN];          write on updHeight(), globalRelabel(), discharge()
    int highest, highestGap, work;
    vec<int> q; //[MAXM << 1];          write on globalRelabel()
    // vec<int> degs;

    HLPP(vec<int> &degrees, int MAXN, int MAXM)
    {
        edges.assign(MAXM << 1, {});
        first_edge.assign(MAXN + 1, 0);
        _cur_edge.assign(MAXN + 1, 0);
        nxt.assign(MAXN, 0);
        lst.assign(MAXN, 0);
        excess.assign(MAXN, 0);
        arc.assign(MAXN, 0);
        gapNxt.assign(MAXN << 1, 0);
        gapPrv.assign(MAXN << 1, 0);
        height.assign(MAXN, 0);
        q.assign(MAXM << 1, 0);

        maxn = degrees.size();
        // assert(maxn <= height.size());
        int cnt(0);
        for (int i(0); i < maxn; ++i)
        {
            first_edge[i] = cnt;
            cnt += degrees[i];
        }
        first_edge[maxn] = cnt;
        copy(first_edge.begin(), first_edge.begin() + maxn + 1, _cur_edge.begin());
    }

    inline void addEdge(int from, int to, int f, bool isDirected = true)
    {
        edges[_cur_edge[from]++] = {to, _cur_edge[to], f};
        edges[_cur_edge[to]++] = {from, _cur_edge[from] - 1, isDirected ? 0 : f};
    }

    inline void pushLst(int h, int v)
    {
        nxt[v] = lst[h];
        lst[h] = v;
    }

    inline void updHeight(int v, int nh)
    {
        if (height[v] != maxn)
        {
            gapNxt[gapPrv[v]] = gapNxt[v];
            gapPrv[gapNxt[v]] = gapPrv[v];
        }

        height[v] = nh;
        if (nh == maxn)
            return;

        highestGap = max(highestGap, nh);
        if (excess[v] > 0)
        {
            highest = max(highest, nh);
            pushLst(nh, v);
        }

        nh += maxn;
        gapNxt[v] = gapNxt[nh];
        gapPrv[v] = nh;
        gapNxt[nh] = v;
        gapPrv[gapNxt[v]] = v;
    }

    inline void globalRelabel(bool from_t = true)
    {
        work = 0;
        fill(height.begin(), height.begin() + maxn, maxn);
        fill(lst.begin(), lst.begin() + maxn, -1);
        iota(gapNxt.begin(), gapNxt.begin() + maxn, 0);
        iota(gapPrv.begin(), gapPrv.begin() + maxn, 0);

        int src = (from_t ? t : s);

        height[src] = 0;
        q[0] = src;
        for (int i(0), sz(1); i < sz; ++i)
        {
            int v = q[i];
            for (int ie = first_edge[v]; ie < first_edge[v + 1]; ++ie)
            {
                auto &e = edges[ie];
                if (height[e.to] == maxn && edges[e.rev].f > 0)
                    q[sz++] = e.to, updHeight(e.to, height[v] + 1);
            }
            // if (from_t == false && v != t)
            highest = highestGap = height[v];
        }
    }

    inline void push(int v, edge &e)
    {
        T df = min(excess[v], e.f);
        if (df > 0)
        {
            if (!excess[e.to])
                pushLst(height[e.to], e.to);
            e.f -= df, edges[e.rev].f += df;
            excess[v] -= df, excess[e.to] += df;
        }
    }

    inline void discharge(int v)
    {
        int nh = maxn;

        for (int i(arc[v]); i < first_edge[v + 1]; ++i)
        {
            auto &e = edges[i];
            if (e.f > 0)
            {
                if (height[v] == height[e.to] + 1)
                {
                    push(v, e);
                    if (excess[v] <= 0)
                    {
                        arc[v] = i;
                        return;
                    }
                }
                else
                    nh = min(nh, height[e.to] + 1);
            }
        }

        for (int i(first_edge[v]); i < arc[v]; ++i)
        {
            auto &e = edges[i];
            if (e.f > 0)
            {
                if (height[v] == height[e.to] + 1)
                {
                    push(v, e);
                    if (excess[v] <= 0)
                    {
                        arc[v] = i;
                        return;
                    }
                }
                else
                    nh = min(nh, height[e.to] + 1);
            }
        }

        ++work;

        if (gapNxt[gapNxt[height[v] + maxn]] != height[v] + maxn)
        {
            updHeight(v, nh);
        }
        else
        {
            int oldH = height[v];
            for (int h(oldH); h < highestGap + 1; ++h)
            {
                for (int i(gapNxt[h + maxn]); i < maxn; i = gapNxt[i])
                {
                    height[i] = maxn;
                }
                gapNxt[h + maxn] = gapPrv[h + maxn] = h + maxn;
            }
            highestGap = oldH - 1;
        }
    }

    inline vec<int> scanMincut()
    {
        vec<int> D(height.size());
        vec<int> Q(height.size());
        Q[0] = s;
        D[s] = 1;
        int top = 0;
        int cap = 1;
        while (top < cap)
        {
            int v = Q[top++];
            for (int ie = first_edge[v]; ie < first_edge[v + 1]; ++ie)
            {
                auto &e = edges[ie];
                if (!D[e.to] && e.f > 0)
                {
                    D[e.to] = D[v] + 1;
                    // if (e.to == t)
                    // throw "destination accessed";
                    Q[cap++] = e.to;
                }
            }
        }
        return D;
    }

    auto cutset(T *d) const
    {

        vec<vec<T>> cset;
        for (T u = 0; u < maxn; ++u)
            for (int ie = first_edge[u]; ie < first_edge[u + 1]; ++ie)
            {
                auto &e = edges[ie];
                    if (u > e.to && (d[u] == 0) != (0 == d[e.to]))
                        cset.push_back({u, e.to});
            }
        return cset;
    }

    inline void recoverUndirected()
    {
        for (auto &e : edges)
            edges[e.rev].f = e.f = (edges[e.rev].f + e.f) / 2;
    }

    inline void returnExcess()
    {
        work = 0;
        globalRelabel(false);
        for (; ~highest; --highest)
        {
            while (~lst[highest])
            {
                int v = lst[highest];
                lst[highest] = nxt[v];
                if (height[v] == highest)
                {
                    discharge(v);
                    if (work > maxn << 2)
                        globalRelabel(false);
                }
            }
        }
    }
    /* prework: call exactly once */
    inline void sortEdges()
    {
        for (int v(0); v < maxn; ++v)
        {
            sort(edges.begin() + first_edge[v], edges.begin() + first_edge[v + 1],
                 [](edge &l, edge &r)
                 { return l.to < r.to; });
            for (int i(first_edge[v]); i < first_edge[v + 1]; ++i)
            {
                auto &e = edges[i];
                edges[e.rev].rev = i;
            }
        }
    }

    inline T calc(int s, int t)
    {
        this->s = s;
        this->t = t;

        copy(first_edge.begin(), first_edge.begin() + maxn, arc.begin());
        fill(excess.begin(), excess.begin() + maxn, 0);
        excess[s] = INF, excess[t] = -INF;
        globalRelabel();

        for (int ie(first_edge[s]); ie < first_edge[s + 1]; ++ie)
            push(s, edges[ie]);

        for (; ~highest; --highest)
        {
            while (~lst[highest])
            {
                int v = lst[highest];
                lst[highest] = nxt[v];
                if (height[v] == highest)
                {
                    discharge(v);
                    if (work > maxn << 2)
                        globalRelabel();
                }
            }
        }

        return excess[t] + INF;
    }
};
HLPP<int> *hp;
int n, m;

vec<int> I, J;

void with_clean_data(int n, int m)
{
    using T = int;
    mp<int, T> P;
    I.assign(n + 1, -1);
    unordered_set<int> S;
    while (m--)
    {
        int u, v;
        T w;
        cin >> u >> v >> w;
        if (u == v)
            continue;
        if (u > v)
            swap(u, v);
        S.emplace(u);
        S.emplace(v);
        auto k = u * (n + 1) + v;
        auto iter = P.find(k);
        if (iter != P.end())
            iter->second += w;
        else
            P[k] = w;
    }
    int ctr = 0;
    J.reserve(S.size());
    for (auto k : S)
    {
        I[k] = ctr++;
        J.emplace_back(k);
    }
    int nn = S.size();
    int mm = P.size();
    vec<int> degs(nn);
    for (auto [k, w] : P)
    {
        int u = k / (n + 1);
        int v = k % (n + 1);
        ++degs[I[u]];
        ++degs[I[v]];
    }
    hp = new HLPP<int>(degs, nn, mm);

    for (auto [k, w] : P)
    {
        int u = k / (n + 1);
        int v = k % (n + 1);
        hp->addEdge(I[u], I[v], w, false);
    }
    hp->sortEdges();
}

vec<int> node, tmp1, tmp2;
vec<vec<int>> ans;
vec<int> d;

void work(int l, int r)
{
    if (l == r)
        return;
    int S = node[l], T = node[l + 1];
    int t = 0;
    int s = node[l], tt = node[l + 1];
    if (I[S] == -1)
    {
        d.assign(n, 0);
        d[S] = 1;
    }
    else if (I[T] == -1)
    {
        d.assign(n, 1);
        d[T] = 0;
    }
    else
    {
        d.assign(n, 0);
        int c = 0;
        t = hp->calc(I[S], I[T]);
        hp->returnExcess();
        for (auto i : hp->scanMincut())
            d[J[c++]] = i;
        hp->recoverUndirected();
        hp->cutset(&d.front());
    }
    ans[T][S] = ans[S][T] = t;
    int cnt1 = 0, cnt2 = 0;
    for (int i = l; i <= r; ++i)
        if (d[node[i]])
            tmp1[++cnt1] = node[i];
        else
            tmp2[++cnt2] = node[i];
    for (int i = 1; i <= cnt1; ++i)
        node[i + l - 1] = tmp1[i];
    for (int i = 1; i <= cnt2; ++i)
        node[cnt1 + l + i - 1] = tmp2[i];
    work(l, l + cnt1 - 1);
    work(l + cnt1, r); //分治
    for (int i = 1; i <= cnt1; ++i)
        for (int j = 1; j <= cnt2; ++j)
        {
            int ii = node[i + l - 1], jj = node[j + cnt1 + l - 1];
            ans[jj][ii] = ans[ii][jj] = min(min(ans[ii][s], ans[s][tt]), ans[tt][jj]);
        } //每个点都要处理
    return;
}

int main()
{
    ios::sync_with_stdio(0);
    cin >> n >> m;
    ++n;
    ans.assign(n, vec<int>(n, ~0u >> 1));
    node.assign(n, 0);
    tmp1.assign(n, 0);
    tmp2.assign(n, 0);
    int x, y;
    with_clean_data(n, m);
    iota(node.begin(), node.end(), 0);
    work(0, n - 1);
    int que;
    cin >> que;
    while (que--)
    {
        cin >> x >> y;
        printf("%d\n", ans[x][y]);
    }
    // delete hp;
    return 0;
}