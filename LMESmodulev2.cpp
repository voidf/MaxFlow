﻿// Large-MediumExcessScalingAlgorithm
// from arxiv:1910.04848 section 4 & 5

#define DBG

#ifndef DBG

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#endif

#include <cassert>
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
            li<vertex> nset[3];
        };

    private:
        using pa = pair<T, T>;
        // A+ 出边 A- 入边
        vector<vector<cached_edge<T, U>>> &rnet;
        unique_ptr<label_info[]> lab;
        unique_ptr<vertex[]> verts;
        unique_ptr<T[]> NextL;
        que<pa> distq;
        T _source, _sink, hivert{0}, _relabel_progress{0}, _relabel_threshold;
        // T lactl{0}, hactl{0}, lactm{0}, hactm{0};
        T MinL{0}, MaxML{0};
        U _max_cap, K, karg, delta, Q, epi, M;

    public:
        max_flow_instance(vector<vector<cached_edge<T, U>>> &graph)
            : rnet(graph),
              lab(make_unique<label_info[]>(rnet.size() + 1)),
              verts(make_unique<vertex[]>(rnet.size())),
              NextL(make_unique<T[]>(rnet.size())),
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
            // 1910.04848 Section 4 Page 5:
            // In order to balance the terms in the running time and optimize the overall running time,
            // Ahuja et al. chose k to be the least power of 2 that exceeds 2 + logU / log(logU)

            karg = static_cast<U>(ceil(pow(2, kpow))); // 论文中的k
            // Section 6
            Q = ceil(log(4 * rnet.size()) / log(karg));
            epi = pow(karg, -Q);
            M = pow(epi, -2);
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

            vec<vec<T>> cset;
            for (int u = 0; u < rnet.size(); ++u)
                for (auto &edge : rnet[u])
                    if (u > edge.to && (d[u] == 0) != (0 == d[edge.to]))
                        cset.push_back({u, edge.to});
            return cset;
        }

        void prepare(T s, T t)
        {
            _source = s;
            _sink = t;
            for (auto &u : rnet)
                for (auto &edge : u)
                    edge.cap = edge.rcap = (edge.cap + edge.rcap) / 2;
            init();
        }

    private:
        static constexpr T ALPHA = 6, BETA = 12;
        static constexpr double GLOBAL_RELABEL_FREQ = 0.5;

        void init()
        {
            for (int u = 0; u < rnet.size() + 1; ++u)
            {
                for (auto &x : lab[u].nset)
                    x.clear();
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
                MinL = rnet.size();
                MaxML = 1;
                for (size_t i = 0; i <= hivert; ++i)
                    for (auto &x : lab[i].nset)
                        x.clear();

                for (size_t i = 0; i < rnet.size(); ++i)
                {
                    if (i != _source && i != _sink && verts[i].label < rnet.size() && verts[i].excess > 0)
                    {
                        if (verts[i].excess >= delta / 2)
                        {
                            MaxML = max(MaxML, verts[i].label);
                            MinL = min(MinL, verts[i].label);
                            lab[verts[i].label].nset[2].push(&verts[i]);
                            continue;
                        }
                        else if (verts[i].excess >= delta / karg)
                        {
                            MaxML = max(MaxML, verts[i].label);
                            lab[verts[i].label].nset[1].push(&verts[i]);
                            continue;
                        }
                    }
                    lab[verts[i].label].nset[0].push(&verts[i]);
                }
                upd_nextL();
                // todo
                // if (MinL == rnet.size())
                // MinL = 0;
                while (MaxML > 0)
                {
                    while (MinL < rnet.size())
                    {
                        // if (MinL > rnet.size())
                        // cerr << "bp0";
                        if (lab[MinL].nset[2].empty())
                        {
                            // auto beforeMinL = MinL;
                            MinL = NextL[MinL];
                            // if (MinL < 0)
                            // cerr << "bp1";
                            // NextL[beforeMinL] = 0;
                            continue;
                        }

                        auto vertex = get_vertex_idx(lab[MinL].nset[2].front());
                        process(vertex, false);

                        if (_relabel_progress * GLOBAL_RELABEL_FREQ >= _relabel_threshold)
                        {
                            _relabel_progress = 0;
                            global_relabel();
                        }
                    }
                    if (lab[MaxML].nset[1].empty())
                    {
                        --MaxML;
                        continue;
                    }

                    auto vertex = get_vertex_idx(lab[MaxML].nset[1].front());
                    process(vertex, true);

                    if (_relabel_progress * GLOBAL_RELABEL_FREQ >= _relabel_threshold)
                    {
                        _relabel_progress = 0;
                        global_relabel();
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
                            lab[label].nset[1].remove(&verts[vertex]);
                            lab[label].nset[0].push(&verts[vertex]);
                            ret = true;
                        }
                    }
                    else
                    {
                        char typ = point_type(vertex);
                        if (typ < 2)
                        {
                            lab[label].nset[2].remove(&verts[vertex]);
                            lab[label].nset[typ].push(&verts[vertex]);
                            ret = true;
                        }
                    }

                    if (edge.to != _source && edge.to != _sink)
                    {
                        char target_point_type_new = point_type(edge.to);
                        if (target_point_type_new != target_point_type_old)
                        {
                            lab[label - 1].nset[target_point_type_old].remove(&verts[edge.to]);
                            if (target_point_type_new == 2)
                            {
                                MinL = label - 1;
                                // if (lab[label].nset[2].empty())
                                // NextL[MinL] = NextL[label];
                                // else
                                NextL[MinL] = label;
                            }
                            lab[label - 1].nset[target_point_type_new].push(&verts[edge.to]);

                            ret = true;
                        }
                        if (target_point_type_new == 2 && med)
                            MinL = label - 1;
                    }
                    if (ret)
                        return true;
                }
            return false;
        }

        inline void relabel(const T vertex, const T current_label, bool med)
        {
            _relabel_progress += BETA;
            // auto new_label = calculate_new_label(vertex);
            T increase_to = rnet.size() - 1;
            for (auto &edge : rnet[vertex])
            {
                if (edge.cap == 0)
                    continue;
                increase_to = min(increase_to, verts[edge.to].label);
            }
            _relabel_progress += rnet[vertex].size();
            auto new_label = increase_to + 1;

            lab[current_label].nset[2 - med].remove(&verts[vertex]);
            verts[vertex].label = new_label;

            if (new_label != rnet.size()) // 并非不可达
            {
                // update NextL, optimize?
                if (!med)
                    for (T c = current_label;; c = NextL[c])
                    {
                        if (NextL[c] > new_label)
                        {
                            NextL[new_label] = NextL[c];
                            NextL[c] = new_label;
                            break;
                        }
                        else if (NextL[c] == new_label)
                            break;
                        // else if (NextL[c] == 0)
                        // {
                        // NextL[c] = new_label;
                        // break;
                        // }
                    }

                hivert = max(hivert, new_label);
                MaxML = max(MaxML, new_label);
                lab[new_label].nset[2 - med].push(&verts[vertex]);
            }

            if (lab[current_label].nset[0].empty() &&
                lab[current_label].nset[1].empty() &&
                lab[current_label].nset[2].empty())
            {
                NextL[current_label] = 0;
                gap_relabel(current_label);
                verts[vertex].label = rnet.size();
            }
        }

        inline void upd_nextL()
        {
            fill(NextL.get(), NextL.get() + rnet.size(), rnet.size());
            T prv = rnet.size();
            for (T i = hivert; i >= 0; --i)
            {
                if (!lab[i].nset[2].empty())
                {
                    NextL[i] = prv;
                    prv = i;
                }
            }
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
                for (auto &x : lab[i].nset)
                    x.clear();
            }

            MaxML = 0;
            MinL = rnet.size();

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
                            char typ = point_type(edge.to);
                            lab[verts[edge.to].label].nset[typ].push(node);
                            switch (typ)
                            {
                            case 0:
                                break;
                            case 2:
                                MinL = min(MinL, verts[edge.to].label);
                            case 1:
                                MaxML = max(MaxML, verts[edge.to].label);
                                break;
                            }
                        }
                    }
            }
            upd_nextL();
        }
        /* 由于gap_height已经没有点了，所以把大于gap_height的节点全都设置为不可达，rnet.size()即为不可达标号 */
        inline void gap_relabel(const T gap_height)
        {
            for (auto chei = gap_height + 1; chei <= hivert; ++chei)
            {
                for (auto &x : lab[chei].nset)
                {
                    while (!x.empty())
                    {
                        auto *ptr = x.pop();
                        auto vertex_idx = get_vertex_idx(ptr);
                        verts[vertex_idx].label = rnet.size();
                    }
                }
                NextL[chei] = rnet.size();
            }
            MaxML = hivert = gap_height - 1;
        }
    };
}

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

    void work(int l, int r, ahuja_orlin::max_flow_instance<int, U> &ao)
    {
        if (l == r)
            return;
        int S = node[l], T = node[l + 1];
        U t = 0;

        fill(d.get(), d.get() + _n, 0);
        int c = 0;
        ao.prepare(S, T);
        t = ao.find_max_flow();
        ao.preflow_to_flow();
        for (auto i : ao.stcut())
            d[c++] = i;
        ans[T][S] = ans[S][T] = t;
        int pos = ansedge[T][S] = ansedge[S][T] = cutsets.size();
        // cerr << "DEBUG: [before] cutsets size:" << cutsets.size() << endl;
        cutsets.emplace_back(ao.cutset(d.get()));
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
        work(l, l + cnt1 - 1, ao);
        work(l + cnt1, r, ao);
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
        ahuja_orlin::max_flow_instance<int, U> ao(G);
        work(0, _n - 1, ao);
    }
};

#ifdef DBG

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
    // auto sol = Solver<double>(3);
    // sol.add_edge(0, 1, .5);
    // sol.add_edge(0, 2, .5);
    // sol.add_edge(1, 2, .5);
    // sol.solve();
    // cerr << sol.ans[0][1] << endl;

    // auto sol2 = Solver<double>(3);
    // sol2.add_edge(0, 1, 1);
    // sol2.add_edge(0, 2, 1);
    // sol2.add_edge(1, 2, 1);
    // sol2.solve();
    // cerr << sol2.ans[0][1] << endl;
    return 0;
}

#else

template <class T>
PyObject *vec2int_tuple(vec<T> &src)
{
    int n = src.size();
    PyObject *uarr = PyTuple_New(n);
    if (!uarr)
        throw logic_error("Unable to allocate memory for 1d array");

    for (int i = 0; i < n; ++i)
    {
        PyObject *num = extractor<T>::pack(src[i]);
        if (!num)
        {
            Py_DECREF(uarr);
            throw logic_error("Unable to allocate memory for Python int");
        }
        PyTuple_SET_ITEM(uarr, i, num);
    }
    return uarr;
}
template <class T>
PyObject *vec2d2int_tuple(vec<vec<T>> &src)
{
    int n = src.size();
    PyObject *ans = PyTuple_New(n);
    if (!ans)
        throw logic_error("Unable to allocate memory for 2d array");
    for (int i = 0; i < n; ++i)
    {
        PyObject *uarr = NULL;
        try
        {
            uarr = vec2int_tuple(src[i]);
        }
        catch (logic_error &e)
        {
            Py_DECREF(ans);
            throw e;
        }
        PyTuple_SET_ITEM(ans, i, uarr);
    }
    return ans;
}

template <class T>
PyObject *vec3d2int_tuple(vec<vec<vec<T>>> &src)
{
    int n = src.size();
    PyObject *ans = PyTuple_New(n);
    if (!ans)
        throw logic_error("Unable to allocate memory for 3d array");
    for (int i = 0; i < n; ++i)
    {
        PyObject *uarr = NULL;
        try
        {
            uarr = vec2d2int_tuple(src[i]);
        }
        catch (logic_error &e)
        {
            Py_DECREF(ans);
            throw e;
        }
        PyTuple_SET_ITEM(ans, i, uarr);
    }
    return ans;
}

template <class T, int Dim = 1>
struct extractor
{
    static T extract(PyObject *obj);
    static PyObject *pack(T val);
};

template <>
static int extractor<int>::extract(PyObject *obj) { return PyLong_AsLong(obj); }
template <>
static double extractor<double>::extract(PyObject *obj) { return PyFloat_AsDouble(obj); }
template <>
static PyObject *extractor<int>::pack(int val) { return PyLong_FromLong(val); }
template <>
static PyObject *extractor<double>::pack(double val) { return PyFloat_FromDouble(val); }

template <class U>
static PyObject *ao_main(PyObject *self, PyObject *args)
{
    // PyObject *result = NULL; // exception
    PyObject *li;
    if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &li))
    {
        PyErr_SetString(PyExc_TypeError, "edge list must be a list of list[u, v, w].");
        return NULL;
    }
    auto m = PyList_Size(li);
    int maxn = 0;
    vec<tuple<int, int, U>> Ebuffer;
    for (Py_ssize_t i = 0; i < m; ++i)
    {
        PyObject *elem = PyList_GetItem(li, i);
        if (!PyList_Check(elem) || PyList_Size(elem) != 3)
        {
            PyErr_SetString(PyExc_TypeError, "edge list element must be a list following the form [u, v, w].");
            return NULL;
        }
        int u, v;
        U w;
        u = PyLong_AsLong(PyList_GetItem(elem, 0));
        v = PyLong_AsLong(PyList_GetItem(elem, 1));
        w = extractor<U>::extract(PyList_GetItem(elem, 2));
        // cerr << "DEBUG: edge(" << u << ", " << v << ", " << w << ")" << endl;
        maxn = max({u, v, maxn});
        Ebuffer.emplace_back(u, v, w);
    }
    auto n = maxn + 1;
    // cerr << "DEBUG: n:" << n << endl;
    auto sol = Solver<U>(n);
    // cerr << "DEBUG: Solver:" << &sol << endl;
    for (auto &[u, v, w] : Ebuffer)
        sol.add_edge(u, v, w);
    // cerr << "DEBUG: before solve" << endl;
    sol.solve();
    // cerr << "DEBUG: after solve" << endl;

    PyObject *ans = vec2d2int_tuple(sol.ans);
    PyObject *ans_cutset_index = vec2d2int_tuple(sol.ansedge);
    PyObject *cutsets = vec3d2int_tuple(sol.cutsets);

    PyObject *ret = PyTuple_New(3);
    if (!ret)
        throw logic_error("Unable to allocate memory for return tuple(ans, ans_cutset_index, cutsets)");
    PyTuple_SET_ITEM(ret, 0, ans);
    PyTuple_SET_ITEM(ret, 1, ans_cutset_index);
    PyTuple_SET_ITEM(ret, 2, cutsets);
    return ret;
}

static char ao_docs[] =
    "ahuja_orlin O(n * m) solve maxflow with residual network \n"
    "output n * n tuple of pair wise mincut answers";

static PyMethodDef cmds[] = {
    {"solve_int", (PyCFunction)ao_main<int>,
     METH_VARARGS /*（METH_VARARGS,MET_KEYWORDSMET_KEYWORDS关键字参数,MET_NOARGS无参数*/,
     ao_docs},
    {"solve_double", (PyCFunction)ao_main<double>, // 如果流值不全是整数，不保证正确性
     METH_VARARGS,
     ao_docs},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef ao_definition = {
    PyModuleDef_HEAD_INIT,
    "ao",
    "ahuja_orlin copied from https://github.com/touqir14/MaxFlow",
    -1,
    cmds};

PyMODINIT_FUNC PyInit_ao(void) // 必须是PyInit_{setup.py中提供的模块名}
{
    Py_Initialize();
    return PyModule_Create(&ao_definition);
}
// C:\\etc\\VisualStudio\\2019\\Community\\VC\\Tools\\MSVC\\14.29.30133\\bin\\HostX86\\x64\\cl.exe

#endif

/*
4 5
1 2 2
2 3 2
4 2 3
4 3 1
1 3 1
3
1 4
2 4
2 3


8 15
1 7 641
1 3 307
1 4 820
1 8 914
2 6 445
2 3 907
2 7 432
3 7 132
3 8 299
3 5 629
4 5 742
4 6 925
4 7 327
5 7 390
7 8 747
28
1 2
1 3
1 4
1 5
1 6
1 7
1 8
2 3
2 4
2 5
2 6
2 7
2 8
3 4
3 5
3 6
3 7
3 8
4 5
4 6
4 7
4 8
5 6
5 7
5 8
6 7
6 8
7 8



4 5
1 4 583
1 2 132
2 3 695
2 4 441
3 4 612
6
1 2
1 3
1 4
2 3
2 4
3 4

*/