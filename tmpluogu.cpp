

// #include <bits/stdc++.h>
#include <memory>
#include <cmath>
#include <cstdint>
#include <utility>
#include <limits>
#include <queue>
#include <vector>
#include <list>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <numeric>
using namespace std;

#define vec vector
#define us unordered_set
#define mp unordered_map

int n, m, s, t;
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

// template <typename T = uint32_t, typename U = uint32_t>
template <typename T, typename U>
struct cached_edge
{
    cached_edge(T to, U c, T rev) : to(to), rev(rev), cap(c), rcap(c) {}
    cached_edge(T to, U c, T rev, U rc) : to(to), rev(rev), cap(c), rcap(rc) {}
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
            li<vertex> largeset, mediumset, smallset;
        };

    private:
        using pa = pair<T, T>;
        vector<vector<cached_edge<T, U>>> &rnet;
        unique_ptr<label_info[]> lab;
        unique_ptr<vertex[]> verts;
        que<pa> distq;
        T _source, _sink, hivert{0}, _relabel_progress{0}, _relabel_threshold;
        T lactl{0}, hactl{0}, lactm{0}, hactm{0};
        U _max_cap, K, karg, delta;

    public:
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

            // for (auto &node_edges : rnet)
            // {
            //     sort(node_edges.begin(), node_edges.end(), [&](const cached_edge<T, U> &a, const cached_edge<T, U> &b) -> bool
            //          { return a.cap + a.rcap > b.cap + b.rcap; });
            //     for (auto &e : node_edges)
            //     {
            //         auto &re = rnet[e.to][(long long)(e.rev)];
            //         long long dist = &e - &node_edges.front();
            //         re.rev = (cached_edge<T, U> *)dist;
            //     }
            // }
            _max_cap = 0;
            for (auto &node_edges : rnet)
                for (auto &e : node_edges)
                {
                    _max_cap = max(_max_cap, e.cap);
                    // auto &re = rnet[e.to][(long long)(e.rev)];
                    // e.rev = &re;
                }
            K = static_cast<U>(ceil(log2(ceil(_max_cap)))); // 流值相关，过大可能有溢出风险。logU的来源，不(能很好)支持浮点
            auto kpow = static_cast<U>(log2(2 + ceil(log2(_max_cap) / log2(log2(_max_cap)))));
            karg = static_cast<U>(ceil(pow(2, kpow)));
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

        auto cutset(T *d) const
        {
            // cerr << "[s,t]:" << _source << ", " << _sink << endl;
            // cerr << "d:";
            // for (int u = 0; u < rnet.size(); ++u)
            // cerr << d[u] << " ";
            // cerr << endl;
            vec<vec<T>> cset;
            for (int u = 0; u < rnet.size(); ++u)
                for (auto &edge : rnet[u])
                {
                    // cerr << "(" << u << ", " << edge.to << "), cap:" << edge.cap << ", rcap:" << edge.rcap << endl;
                    // if (edge.cap == 0)
                    if (u > edge.to && (d[u] == 0) != (0 == d[edge.to]))
                        cset.push_back({u, edge.to});
                }
            return cset;
        }

        void prepare(T s, T t)
        {
            _source = s;
            _sink = t;
            // for (auto &u : rnet)
            //     for (auto &edge : u)
            //         edge.cap = edge.rcap = (edge.cap + edge.rcap) / 2;
            init();
        }

    private:
        static constexpr T ALPHA = 6, BETA = 12;
        static constexpr double GLOBAL_RELABEL_FREQ = 0.5;

        void init()
        {
            for (int u = 0; u < rnet.size() + 1; ++u)
            {
                lab[u].largeset.clear();
                lab[u].mediumset.clear();
                lab[u].smallset.clear();
            }
            // _max_cap = 0;

            for (auto &edge : rnet[_source]) // 先把流从源点提出来
            {
                // _max_cap = max(_max_cap, edge.cap + 1);

                verts[edge.to].excess = edge.cap;
                edge.rcap += edge.cap;
                rnet[edge.to][edge.rev].cap += edge.cap;
                rnet[edge.to][edge.rev].rcap -= edge.cap;
                edge.cap = 0;
            }
            hivert = 1;
        }

        void find_max_flow_inner()
        {

            global_relabel();
            for (delta = static_cast<U>(pow(2, K)); delta >= 1; delta /= karg)
            {
                // auto delta = static_cast<U>(pow(2, K - k)); // 流值相关

                lactm = lactl = rnet.size();
                hactm = hactl = 1;
                // if (_source == 1 && _sink == 2)
                // cerr << "breakpoint";
                for (size_t i = 0; i <= hivert; ++i)
                {
                    lab[i].largeset.clear();
                    lab[i].mediumset.clear();
                    lab[i].smallset.clear();
                }

                for (size_t i = 0; i < rnet.size(); ++i)
                {
                    if (i != _source && i != _sink && verts[i].label < rnet.size() && verts[i].excess > 0)
                    {
                        if (verts[i].excess >= delta / 2)
                        {
                            lactl = min(lactl, verts[i].label);
                            hactl = max(hactl, verts[i].label);
                            lab[verts[i].label].largeset.push(&verts[i]);
                            continue;
                        }
                        else if (verts[i].excess >= delta / karg)
                        {
                            lactm = min(lactm, verts[i].label);
                            hactm = max(hactm, verts[i].label);
                            lab[verts[i].label].mediumset.push(&verts[i]);
                            continue;
                        }
                    }
                    // else
                    lab[verts[i].label].smallset.push(&verts[i]);
                }
                while (lactl <= hactl || lactm <= hactm)
                {
                    while (lactl <= hactl)
                    {
                        if (lab[lactl].largeset.empty())
                        {
                            ++lactl;
                            continue;
                        }

                        auto vertex = get_vertex_idx(lab[lactl].largeset.front());
                        // if (vertex == 5)
                        // cerr << "breakpoint2";
                        process(vertex, false);

                        if (_relabel_progress * GLOBAL_RELABEL_FREQ >= _relabel_threshold)
                        {
                            _relabel_progress = 0;
                            global_relabel();
                        }
                    }
                    if (lactm <= hactm)
                    {
                        if (lab[hactm].mediumset.empty())
                        {
                            --hactm;
                            continue;
                        }

                        auto vertex = get_vertex_idx(lab[hactm].mediumset.front());
                        // if (vertex == 8)
                        // cerr << "bp3";
                        process(vertex, true);

                        if (_relabel_progress * GLOBAL_RELABEL_FREQ >= _relabel_threshold)
                        {
                            _relabel_progress = 0;
                            global_relabel();
                        }
                    }
                }
            }
        }

        T get_vertex_idx(vertex *n) { return distance(verts.get(), n); }

        inline void process(const T vertex, const bool med)
        {
            const auto label = verts[vertex].label;
            if (push(vertex, label, med))
                return;
            relabel(vertex, label, med); // 没推到至少delta/2
        }

        inline char point_type(const T vertex) const
        {
            if (verts[vertex].excess == 0)
                return 0;
            else if (verts[vertex].excess >= delta / 2)
                return 2;
            else if (verts[vertex].excess >= delta / karg)
                return 1;
            return 0;
        }

        /* 非饱和推流，不是推完vertex的excess，而是推到还剩delta/2以下 */
        inline bool push(const T vertex, const T label, const bool med)
        {
            for (auto &edge : rnet[vertex])
                if (edge.cap > 0 && label == verts[edge.to].label + 1)
                {
                    auto flow = min(verts[vertex].excess, edge.cap);
                    if (edge.to != _sink)
                        flow = min(flow, delta - verts[edge.to].excess); // 推向的点最多持有delta的超额流
                    char target_point_type_old = point_type(edge.to);

                    verts[vertex].excess -= flow;
                    verts[edge.to].excess += flow;
                    edge.cap -= flow;
                    edge.rcap += flow;
                    rnet[edge.to][edge.rev].rcap -= flow;
                    rnet[edge.to][edge.rev].cap += flow;

                    bool ret = false;
                    if (med)
                    {
                        if (verts[vertex].excess < delta / karg || verts[vertex].excess == 0)
                        {
                            lab[label].mediumset.remove(&verts[vertex]);
                            lab[label].smallset.push(&verts[vertex]);
                            ret = true;
                        }
                    }
                    else
                    {
                        if (verts[vertex].excess == 0)
                        {
                            lab[label].largeset.remove(&verts[vertex]);
                            lab[label].smallset.push(&verts[vertex]);
                            ret = true;
                        }
                        else if (verts[vertex].excess < delta / 2)
                        {
                            lab[label].largeset.remove(&verts[vertex]);
                            if (verts[vertex].excess >= delta / karg)
                                lab[label].mediumset.push(&verts[vertex]);
                            else
                                lab[label].smallset.push(&verts[vertex]);

                            ret = true;
                        }
                    }

                    if (edge.to != _source && edge.to != _sink)
                    {
                        char target_point_type_new = point_type(edge.to);
                        if (target_point_type_new != target_point_type_old)
                        {
                            if (target_point_type_old == 0)
                                lab[label - 1].smallset.remove(&verts[edge.to]);
                            else if (target_point_type_old == 1)
                                lab[label - 1].mediumset.remove(&verts[edge.to]);

                            if (target_point_type_new == 1)
                            {
                                lab[label - 1].mediumset.push(&verts[edge.to]);
                                lactm = min(lactm, label - 1);
                            }
                            else if (target_point_type_new == 2)
                            {
                                lab[label - 1].largeset.push(&verts[edge.to]);
                                lactl = min(lactl, label - 1);
                            }
                            ret = true;
                        }
                        // if (verts[edge.to].excess >= delta / 2)
                        // {
                        //     lab[label - 1].smallset.remove(&verts[edge.to]);
                        //     lab[label - 1].largeset.push(&verts[edge.to]);
                        //     lactl = min(lactl, label - 1);
                        //     // --lactl;
                        // }
                        // else if (verts[edge.to].excess >= delta / karg)
                        // {
                        //     lab[label - 1].smallset.remove(&verts[edge.to]);
                        //     lab[label - 1].mediumset.push(&verts[edge.to]);
                        //     lactm = min(lactm, label - 1);
                        //     // --lactm;
                        //     ret = true;
                        // }
                    }
                    if (ret)
                        return true;
                }
            return false;
        }

        inline void relabel(const T vertex, const T current_label, bool med)
        {
            _relabel_progress += BETA;
            auto new_label = calculate_new_label(vertex);
            if (med)
                lab[current_label].mediumset.remove(&verts[vertex]);
            else
                lab[current_label].largeset.remove(&verts[vertex]);
            verts[vertex].label = new_label;

            if (new_label != rnet.size()) // 并非不可达
            {
                hivert = max(hivert, new_label);
                if (med)
                {
                    hactm = max(hactm, new_label);
                    lab[new_label].mediumset.push(&verts[vertex]);
                }
                else
                {
                    hactl = max(hactl, new_label);
                    lab[new_label].largeset.push(&verts[vertex]);
                }
            }

            if (lab[current_label].largeset.empty() &&
                lab[current_label].mediumset.empty() &&
                lab[current_label].smallset.empty())
            {
                gap_relabel(current_label);
                verts[vertex].label = rnet.size();
            }
        }
        /* 抬到正好能漏到一个可用的邻接点 */
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

        void global_relabel()
        {
            auto not_reached = rnet.size();
            for (size_t i = 0; i < rnet.size(); ++i)
                verts[i].label = not_reached;
            distq.reset();
            distq.push(make_pair(_sink, 0));
            verts[_sink].label = 0;

            for (size_t i = 0; i <= hivert; ++i)
            {
                lab[i].largeset.clear();
                lab[i].mediumset.clear();
                lab[i].smallset.clear();
            }

            lactm = lactl = rnet.size();
            hactm = hactl = hivert = 1;

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
                            if (verts[edge.to].excess == 0)
                                lab[verts[edge.to].label].smallset.push(node);
                            else if (verts[edge.to].excess >= delta / 2) // 流值相关，只贴大于delta/2的
                            {
                                lactl = min(lactl, verts[edge.to].label);
                                hactl = max(hactl, verts[edge.to].label);
                                lab[verts[edge.to].label].largeset.push(node);
                            }
                            else if (verts[edge.to].excess >= delta / karg) // 流值相关，只贴大于delta/2的
                            {
                                lactm = min(lactm, verts[edge.to].label);
                                hactm = max(hactm, verts[edge.to].label);
                                lab[verts[edge.to].label].mediumset.push(node);
                            }
                            else
                                lab[verts[edge.to].label].smallset.push(node);
                        }
                    }
            }
        }
        /* 由于gap_height已经没有点了，所以把大于gap_height的节点全都设置为不可达，rnet.size()即为不可达标号 */
        void gap_relabel(const T gap_height)
        {
            for (auto chei = gap_height + 1; chei <= hivert; ++chei)
            {
                while (!lab[chei].largeset.empty())
                {
                    auto *ptr = lab[chei].largeset.pop();
                    auto vertex_idx = get_vertex_idx(ptr);
                    verts[vertex_idx].label = rnet.size();
                }
                while (!lab[chei].mediumset.empty())
                {
                    auto *ptr = lab[chei].mediumset.pop();
                    auto vertex_idx = get_vertex_idx(ptr);
                    verts[vertex_idx].label = rnet.size();
                }
                while (!lab[chei].smallset.empty())
                {
                    auto *ptr = lab[chei].smallset.pop();
                    auto vertex_idx = get_vertex_idx(ptr);
                    verts[vertex_idx].label = rnet.size();
                }
            }
            hivert = hactm = hactl = gap_height - 1;
        }
    };
}

vec<vec<cached_edge<int, long long>>> G;
vec<int> I, J;

void with_clean_data(int n, int m)
{
    using T = long long;
    mp<int, T> P;
    I.assign(n + 1, -1);
    us<int> S;
    while (m--)
    {
        int u, v;
        T w;
        cin >> u >> v >> w;
        if (u == v)
            continue;
        // if (u > v)
        // swap(u, v);
        S.emplace(u);
        S.emplace(v);
        auto k = u * (n + 2) + v;
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
        int u = k / (n + 2);
        int v = k % (n + 2);
        if (u < v)
        {
            int iu = I[u];
            int iv = I[v];
            auto it = P.find(v * (n + 2) + u);
            T rc;
            if (it != P.end())
                rc = it->second;
            else
                rc = 0;
            G[iu].emplace_back(iv, w, G[iv].size(), rc);
            G[iv].emplace_back(iu, rc, G[iu].size() - 1, w);
        }
        else
        {
            int iu = I[u];
            int iv = I[v];
            auto it = P.find(v * (n + 2) + u);
            T rc;
            if (it != P.end())
                continue;
            else
                rc = 0;
            G[iu].emplace_back(iv, w, G[iv].size(), rc);
            G[iv].emplace_back(iu, rc, G[iu].size() - 1, w);
        }
    }
}

int main()
{
    ios::sync_with_stdio(0);
    cin >> n >> m >> s >> t;
    if (s == t)
    {
        cout << 0 << '\n';
        return 0;
    }
    int x, y;
    with_clean_data(n, m);
    auto ao = ahuja_orlin::max_flow_instance<int, long long>(G);
    ao.prepare(I[s], I[t]);
    cout << ao.find_max_flow() << '\n';
    return 0;
}