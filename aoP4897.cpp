// #include <bits/stdc++.h>
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
#define us unordered_set
#define mp unordered_map

char buf[1 << 23], *p1 = buf, *p2 = buf, obuf[1 << 23], *O = obuf;
#define getchar() (p1 == p2 && (p2 = (p1 = buf) + fread(buf, 1, 1 << 21, stdin), p1 == p2) ? EOF : *p1++)

template <class T>
inline void qr(T &n)
{
    n = 0;
    int c;

    while (!isdigit(c = getchar()))
        ;

    while (isdigit(c))
    {
        n = (n * 10) + (c ^ 0x30);
        c = getchar();
    }
}

template <typename T>
class que
{
    unique_ptr<T[]> _data{nullptr};
    size_t _back{0};
    size_t _front{0};

public:
    que() = default;
    explicit que(std::size_t size) : _data(std::make_unique<T[]>(size)) {}
    que(que<T> &&other) noexcept = default;
    que<T> &operator=(que<T> &&other) noexcept = default;
    que(const que<T> &other) = delete;
    que<T> &operator=(const que<T> &other) = delete;
    void push(T val) noexcept { _data[_back++] = val; }
    T pop() noexcept { return _data[_front++]; }
    bool empty() const noexcept { return _front == _back; }
    void reset() noexcept { _front = _back = 0; }
};

template <typename node>
struct li
{
    node _head = node{}, _tail = node{};
    size_t _size{0};
    li()
    {
        _head.next = &_tail;
        _tail.prev = &_head;
    }
    node *pop()
    {
        auto *ret = _head.next;
        _head.next = _head.next->next;
        _head.next->prev = &_head;
        --_size;
        return ret;
    }
    void push(node *n)
    {
        n->next = &_tail;
        n->prev = _tail.prev;
        _tail.prev->next = n;
        _tail.prev = n;
        ++_size;
    }
    void remove(node *n)
    {
        n->prev->next = n->next;
        n->next->prev = n->prev;
        --_size;
    }
    node *front() const { return _head.next; }
    bool empty() const { return _size == 0; }
    void clear()
    {
        _head.next = &_tail;
        _tail.prev = &_head;
        _size = 0;
    }
    size_t size() const { return _size; }
};

template <typename T, typename U>
struct cached_edge
{
    cached_edge(T to, U c, T rev) : to(to), rev(rev), cap(c), rcap(c) {}
    bool operator==(const cached_edge<T, U> &other) { return to == other.to && rev == other.rev && cap == other.cap && rcap == other.rcap; };
    T to, rev;
    U cap, rcap;
};

namespace ahuja_orlin
{
    template <typename T, typename U>
    class max_flow_instance
    {
        struct vertex
        {
            vertex *next, *prev;
            U excess{0};
            T label;
        };

        struct label_info
        {
            li<vertex> act, iact;
        };

    private:
        using pa = pair<T, T>;
        vector<vector<cached_edge<T, U>>> &rnet;
        unique_ptr<label_info[]> lab;
        unique_ptr<vertex[]> verts;
        que<pa> distq;
        T _source, _sink, lact{0}, hact{0}, hivert{0}, _relabel_progress{0}, _relabel_threshold;
        U _max_cap;

        // statistics
        uint64_t _push_cnt{0}, _relabel_cnt{0}, _gap_cnt{0}, _gap_nodes{0}, _global_relabel_cnt{0};

    public:
        /* call prepare before use */
        max_flow_instance(vector<vector<cached_edge<T, U>>> &graph)
            : rnet(graph),
              lab(make_unique<label_info[]>(rnet.size() + 1)),
              verts(make_unique<vertex[]>(rnet.size())),
              distq(que<pa>{rnet.size()})
        {
            T m = 0;
            for (size_t i = 0; i < rnet.size(); ++i)
                m += rnet[i].size();
            _relabel_threshold = rnet.size() * ALPHA + m / 2;
        }
        max_flow_instance() {}
        U find_max_flow()
        {
            find_max_flow_inner();
            return verts[_sink].excess;
        }
        void preflow_to_flow()
        {
            swap(_source, _sink);
            hivert = rnet.size();
            find_max_flow_inner();
            swap(_source, _sink);
            verts[_source].excess = verts[_sink].excess = 0;
        }

        auto stcut() const
        {
            vec<int> D(rnet.size());
            vec<int> Q(rnet.size());
            Q[0] = _source;
            D[_source] = 1;
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
        void prepare(T s, T t)
        {
            _source = s;
            _sink = t;
            for (auto &u : rnet)
                for (auto &edge : u)
                    edge.cap = edge.rcap = edge.cap + edge.rcap >> 1;
            init();
        }

    private:
        static constexpr T ALPHA = 6, BETA = 12;
        static constexpr double GLOBAL_RELABEL_FREQ = 0.5;

        void init()
        {
            for (int u = 0; u < rnet.size() + 1; ++u)
            {
                lab[u].act.clear();
                lab[u].iact.clear();
            }

            _max_cap = 0;
            for (auto &edge : rnet[_source])
            {
                _max_cap = max(_max_cap, edge.cap);
                verts[edge.to].excess = edge.cap;
                edge.rcap += edge.cap;
                rnet[edge.to][edge.rev].cap += edge.cap;
                rnet[edge.to][edge.rev].rcap -= edge.cap;
                edge.cap = 0;
            }
            _push_cnt = rnet[_source].size();
            hivert = 1;
        }

        void find_max_flow_inner()
        {
            auto K = static_cast<U>(ceil(log2(_max_cap)));
            global_relabel(_max_cap);
            for (U k = 0; k <= K; ++k)
            {
                auto delta = static_cast<U>(pow(2, K - k));

                lact = rnet.size();
                hact = 1;

                for (size_t i = 0; i <= hivert; ++i)
                {
                    lab[i].act.clear();
                    lab[i].iact.clear();
                }

                for (size_t i = 0; i < rnet.size(); ++i)
                {
                    if (verts[i].excess > delta / 2 && i != _source &&
                        i != _sink && verts[i].label < rnet.size())
                    {
                        lact = min(lact, verts[i].label);
                        hact = max(hact, verts[i].label);
                        lab[verts[i].label].act.push(&verts[i]);
                    }
                    else
                        lab[verts[i].label].iact.push(&verts[i]);
                }

                while (lact <= hact)
                {
                    if (lab[lact].act.empty())
                    {
                        ++lact;
                        continue;
                    }

                    auto vertex = get_vertex_idx(lab[lact].act.front());
                    process(vertex, delta);

                    if (_relabel_progress * GLOBAL_RELABEL_FREQ >= _relabel_threshold)
                    {
                        _relabel_progress = 0;
                        global_relabel(delta);
                    }
                }
            }
        }

        T get_vertex_idx(vertex *n) { return distance(verts.get(), n); }

        inline void process(const T vertex, const U delta)
        {
            const auto label = verts[vertex].label;
            if (push(vertex, label, delta))
                return;
            relabel(vertex, label);
        }

        inline bool push(const T vertex, const T label, const U delta)
        {
            for (auto &edge : rnet[vertex])
                if (edge.cap > 0 && label == verts[edge.to].label + 1)
                {
                    ++_push_cnt;
                    auto flow = min(verts[vertex].excess, edge.cap);
                    if (edge.to != _sink)
                        flow = min(flow, delta - verts[edge.to].excess);

                    verts[vertex].excess -= flow;
                    verts[edge.to].excess += flow;
                    edge.cap -= flow;
                    edge.rcap += flow;
                    rnet[edge.to][edge.rev].rcap -= flow;
                    rnet[edge.to][edge.rev].cap += flow;

                    bool ret = false;
                    if (verts[vertex].excess <= delta / 2)
                    {
                        lab[label].act.remove(&verts[vertex]);
                        lab[label].iact.push(&verts[vertex]);
                        ret = true;
                    }

                    if (verts[edge.to].excess > delta / 2 && edge.to != _source &&
                        edge.to != _sink)
                    {
                        lab[label - 1].iact.remove(&verts[edge.to]);
                        lab[label - 1].act.push(&verts[edge.to]);
                        --lact;
                        ret = true;
                    }
                    if (ret)
                        return true;
                }
            return false;
        }

        inline void relabel(const T vertex, const T current_label)
        {
            ++_relabel_cnt;
            _relabel_progress += BETA;
            auto new_label = calculate_new_label(vertex);
            lab[current_label].act.remove(&verts[vertex]);
            verts[vertex].label = new_label;

            if (new_label != rnet.size())
            {
                hivert = max(hivert, new_label);
                hact = max(hact, new_label);
                lab[new_label].act.push(&verts[vertex]);
            }

            if (lab[current_label].act.empty() &&
                lab[current_label].iact.empty())
            {
                gap_relabel(current_label);
                verts[vertex].label = rnet.size();
            }
        }

        inline T calculate_new_label(const T vertex)
        {
            T increase_to = rnet.size() - 1;
            for (auto &edge : rnet[vertex])
            {
                if (edge.cap == 0)
                    continue;
                increase_to = min(increase_to, verts[edge.to].label);
            }
            _relabel_progress += rnet[vertex].size();
            return increase_to + 1;
        }

        void global_relabel(const U delta)
        {
            ++_global_relabel_cnt;
            auto not_reached = rnet.size();
            for (size_t i = 0; i < rnet.size(); ++i)
                verts[i].label = not_reached;
            distq.reset();
            distq.push(make_pair(_sink, 0));
            verts[_sink].label = 0;

            for (size_t i = 0; i <= hivert; ++i)
            {
                lab[i].act.clear();
                lab[i].iact.clear();
            }

            lact = rnet.size();
            hact = hivert = 1;

            while (!distq.empty())
            {
                auto cur = distq.pop();
                auto current_vertex = cur.first;
                auto current_distance = cur.second;
                hivert = max(hivert, current_distance);
                for (auto &edge : rnet[current_vertex])
                    if (edge.rcap > 0 && verts[edge.to].label == not_reached)
                    {
                        verts[edge.to].label = current_distance + 1;
                        distq.push(make_pair(edge.to, current_distance + 1));
                        if (edge.to != _source)
                        {
                            auto *node = &verts[edge.to];
                            if (verts[edge.to].excess > delta / 2)
                            {
                                lact = min(lact, verts[edge.to].label);
                                hact = max(hact, verts[edge.to].label);
                                lab[verts[edge.to].label].act.push(node);
                            }
                            else
                                lab[verts[edge.to].label].iact.push(node);
                        }
                    }
            }
        }

        void gap_relabel(const T gap_height)
        {
            ++_gap_cnt;
            for (auto chei = gap_height + 1; chei <= hivert; ++chei)
            {
                while (!lab[chei].act.empty())
                {
                    ++_gap_nodes;
                    auto *ptr = lab[chei].act.pop();
                    auto vertex_idx = get_vertex_idx(ptr);
                    verts[vertex_idx].label = rnet.size();
                }
                while (!lab[chei].iact.empty())
                {
                    ++_gap_nodes;
                    auto *ptr = lab[chei].iact.pop();
                    auto vertex_idx = get_vertex_idx(ptr);
                    verts[vertex_idx].label = rnet.size();
                }
            }
            hivert = hact = gap_height - 1;
        }
    };
}

struct Solver
{
    vec<vec<cached_edge<int, int>>> G;
    vec<int> I, J;

    void with_clean_data(int n, int m)
    {
        using T = int;
        mp<int, T> P;
        I.assign(n + 1, -1);
        us<int> S;
        while (m--)
        {
            int u, v;
            T w;
            qr(u);
            qr(v);
            qr(w);
            // cin >> u >> v >> w;
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
        G.assign(nn, {});
        for (auto [k, w] : P)
        {
            int u = k / (n + 1);
            int v = k % (n + 1);
            int iu = I[u];
            int iv = I[v];
            G[iu].emplace_back(iv, w, G[iv].size());
            G[iv].emplace_back(iu, w, G[iu].size() - 1);
        }
    }
    unique_ptr<int[]> node, tmp1, d;
    vec<vec<int>> ans;
    int _n, _m;
    Solver(int n, int m) : node(make_unique<int[]>(n)), tmp1(make_unique<int[]>(n)), d(make_unique<int[]>(n)),
                           //    ans(n, vec<int>(n, ~0u >> 1)),
                           _n(n), _m(m)
    {
        init();
    }
    void init()
    {
        iota(node.get(), node.get() + _n, 0);
        ans.assign(_n, vec<int>(_n, ~0u >> 1));
    }
    void work(int l, int r, ahuja_orlin::max_flow_instance<int, int> &ao)
    {
        if (l == r)
            return;
        int S = node[l], T = node[l + 1];
        int t = 0;
        int s = node[l], tt = node[l + 1];
        if (I[S] == -1)
        {
            fill(d.get(), d.get() + _n, 0);
            d[S] = 1;
        }
        else if (I[T] == -1)
        {
            fill(d.get(), d.get() + _n, 1);
            d[T] = 0;
        }
        else
        {
            fill(d.get(), d.get() + _n, 0);
            int c = 0;
            ao.prepare(I[S], I[T]);
            t = ao.find_max_flow();
            ao.preflow_to_flow();
            for (auto i : ao.stcut())
                d[J[c++]] = i;
        }
        ans[T][S] = ans[S][T] = t;
        int cnt1 = 0, cnt2 = r - l;
        for (int i = l; i <= r; ++i)
            if (d[node[i]])
                tmp1[cnt1++] = node[i];
            else
                tmp1[cnt2--] = node[i];
        copy(tmp1.get(), tmp1.get() + r - l + 1, node.get() + l);
        cnt2 = r - l - cnt2;
        work(l, l + cnt1 - 1, ao);
        work(l + cnt1, r, ao);
        for (int i = 1; i <= cnt1; ++i)
            for (int j = 1; j <= cnt2; ++j)
            {
                int ii = node[i + l - 1], jj = node[j + cnt1 + l - 1];
                ans[jj][ii] = ans[ii][jj] = min(min(ans[ii][s], ans[s][tt]), ans[tt][jj]);
            }
        return;
    } //-fsanitize=address
    void solve()
    {
        with_clean_data(_n, _m);
        ahuja_orlin::max_flow_instance<int, int> ao(G);
        work(0, _n - 1, ao);
        int que, x, y, aq;
        qr(que);
        aq = que;
        while (que--)
        {
            qr(x);
            qr(y);
            printf("%d\n", ans[x][y]);
        }
    }
};

int main()
{
    int n, m;
    qr(n);
    qr(m);
    ++n;
    auto sol = Solver(n, m);
    sol.solve();
    return 0;
}