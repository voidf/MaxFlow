﻿#define DBG

#ifndef DBG

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#endif

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
    int _back{0};
    int _front{0};

public:
    que() = default;
    explicit que(int size) : _data(std::make_unique<T[]>(size)) {}
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
    cached_edge() = default;
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
        unique_ptr<label_info[]> lab;
        unique_ptr<vertex[]> verts;
        que<pa> distq;
        T _source, _sink, lact{0}, hact{0}, hivert{0}, _relabel_progress{0}, _relabel_threshold, n;
        U _max_cap;

        vec<cached_edge<T, U>> &R;
        T *eidx;
        T *deg;

        // statistics
        uint64_t _push_cnt{0}, _relabel_cnt{0}, _gap_cnt{0}, _gap_nodes{0}, _global_relabel_cnt{0};

    public:
        max_flow_instance(T n, T *_deg, T *_eid, vector<cached_edge<T, U>> &graph)
            : R(graph),
              deg(_deg),
              n(n),
              eidx(_eid),
              lab(make_unique<label_info[]>(n + 1)),
              verts(make_unique<vertex[]>(n)),
              distq(que<pa>{n})
        {
            _relabel_threshold = n * ALPHA + R.size() / 2;
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
            hivert = n;
            find_max_flow_inner();
            swap(_source, _sink);
            verts[_source].excess = verts[_sink].excess = 0;
        }

        auto stcut() const
        {
            vec<int> D(n);
            vec<int> Q(n);
            Q[0] = _source;
            D[_source] = 1;
            int top = 0;
            int cap = 1;
            while (top < cap)
            {
                int v = Q[top++];
                for (T ei = eidx[v]; ei < eidx[v + 1]; ++ei)
                {
                    auto &e = R[ei];
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
            for (T ei = 0, u = 0; ei < R.size(); ++ei)
            {
                while (ei >= eidx[u + 1])
                    ++u;
                auto &edge = R[ei];
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
            for (auto &edge : R)
                edge.cap = edge.rcap = (edge.cap + edge.rcap) / 2;
            init();
        }

    private:
        static constexpr T ALPHA = 6, BETA = 12;
        static constexpr double GLOBAL_RELABEL_FREQ = 0.5;

        void init()
        {
            for (int u = 0; u < n + 1; ++u)
            {
                lab[u].act.clear();
                lab[u].iact.clear();
            }

            _max_cap = 0;
            for (T ei = eidx[_source]; ei < eidx[_source + 1]; ++ei) // 先把流从源点提出来
            {
                auto &edge = R[ei];
                _max_cap = max(_max_cap, edge.cap);
                verts[edge.to].excess = edge.cap;
                edge.rcap += edge.cap;
                R[edge.rev].cap += edge.cap;
                R[edge.rev].rcap -= edge.cap;
                edge.cap = 0;
            }
            _push_cnt = deg[_source];
            hivert = 1;
        }

        void find_max_flow_inner()
        {
            auto K = static_cast<U>(ceil(log2(ceil(_max_cap)))); // 流值相关，过大可能有溢出风险。logU的来源，不(能很好)支持浮点
            global_relabel(_max_cap);
            for (U k = 0; k <= K; ++k)
            {
                auto delta = static_cast<U>(pow(2, K - k)); // 流值相关

                lact = n;
                hact = 1;

                for (size_t i = 0; i <= hivert; ++i)
                {
                    lab[i].act.clear();
                    lab[i].iact.clear();
                }

                for (size_t i = 0; i < n; ++i)
                {
                    if (verts[i].excess > delta / 2 && i != _source &&
                        i != _sink && verts[i].label < n)
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
            relabel(vertex, label); // 没推到至少delta/2
        }
        /* 非饱和推流，不是推完vertex的excess，而是推到还剩delta/2以下 */
        inline bool push(const T vertex, const T label, const U delta)
        {
            for (T ei = eidx[vertex]; ei < eidx[vertex + 1]; ++ei)
            {
                auto &edge = R[ei];
                if (edge.cap > 0 && label == verts[edge.to].label + 1)
                {
                    ++_push_cnt;
                    auto flow = min(verts[vertex].excess, edge.cap);
                    if (edge.to != _sink)
                        flow = min(flow, delta - verts[edge.to].excess); // 推向的点最多持有delta的超额流

                    verts[vertex].excess -= flow;
                    verts[edge.to].excess += flow;
                    edge.cap -= flow;
                    edge.rcap += flow;
                    R[edge.rev].rcap -= flow;
                    R[edge.rev].cap += flow;

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

            if (new_label != n) // 并非不可达
            {
                hivert = max(hivert, new_label);
                hact = max(hact, new_label);
                lab[new_label].act.push(&verts[vertex]);
            }

            if (lab[current_label].act.empty() &&
                lab[current_label].iact.empty())
            {
                gap_relabel(current_label);
                verts[vertex].label = n;
            }
        }
        /* 抬到正好能漏到一个可用的邻接点 */
        inline T calculate_new_label(const T vertex)
        {
            T increase_to = n - 1;
            for (T ei = eidx[vertex]; ei < eidx[vertex + 1]; ++ei)
            {
                auto &edge = R[ei];
                if (edge.cap == 0)
                    continue;
                increase_to = min(increase_to, verts[edge.to].label);
            }
            _relabel_progress += deg[vertex];
            return increase_to + 1;
        }

        void global_relabel(const U delta)
        {
            ++_global_relabel_cnt;
            auto not_reached = n;
            for (size_t i = 0; i < n; ++i)
                verts[i].label = not_reached;
            distq.reset();
            distq.push(make_pair(_sink, 0));
            verts[_sink].label = 0;

            for (size_t i = 0; i <= hivert; ++i)
            {
                lab[i].act.clear();
                lab[i].iact.clear();
            }

            lact = n;
            hact = hivert = 1;

            while (!distq.empty())
            {
                auto cur = distq.pop();
                auto current_vertex = cur.first;
                auto current_distance = cur.second;
                hivert = max(hivert, current_distance);
                for (T ei = eidx[current_vertex]; ei < eidx[current_vertex + 1]; ++ei)
                {
                    auto &edge = R[ei];
                    if (edge.rcap > 0 && verts[edge.to].label == not_reached)
                    {
                        verts[edge.to].label = current_distance + 1;
                        distq.push(make_pair(edge.to, current_distance + 1));
                        if (edge.to != _source)
                        {
                            auto *node = &verts[edge.to];
                            if (verts[edge.to].excess > delta / 2) // 流值相关，只贴大于delta/2的
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
        }
        /* 由于gap_height已经没有点了，所以把大于gap_height的节点全都设置为不可达，rnet.size()即为不可达标号 */
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
                    verts[vertex_idx].label = n;
                }
                while (!lab[chei].iact.empty())
                {
                    ++_gap_nodes;
                    auto *ptr = lab[chei].iact.pop();
                    auto vertex_idx = get_vertex_idx(ptr);
                    verts[vertex_idx].label = n;
                }
            }
            hivert = hact = gap_height - 1;
        }
    };
}

template <class U>
struct Solver
{
    vec<cached_edge<int, U>> edg;
    vec<tuple<int, int, U>> tmpedge;

    inline void add_edge(int u, int v, U w)
    {
        tmpedge.emplace_back(u, v, w);
        ++deg[u];
        ++deg[v];
    }

    unique_ptr<int[]> node, tmp1, d, deg, eid, ceid;
    vec<vec<U>> ans;
    vec<vec<int>> ansedge;
    vec<vec<vec<int>>> cutsets;
    int _n, _m;
    Solver(int n, int m) : node(make_unique<int[]>(n)),
                           tmp1(make_unique<int[]>(n)),
                           d(make_unique<int[]>(n)),
                           deg(make_unique<int[]>(n)),
                           eid(make_unique<int[]>(n + 1)),
                           ceid(make_unique<int[]>(n)),
                           //    ans(n, vec<int>(n, ~0u >> 1)),
                           _n(n),
                           _m(m)
    {
        init();
    }
    void init()
    {
        iota(node.get(), node.get() + _n, 0);
        fill(deg.get(), deg.get() + _n, 0);
        ans.assign(_n, vec<U>(_n, ~0u >> 1));
        ansedge.assign(_n, vec<int>(_n, ~0u >> 1));
        edg.resize(_m * 2);
        tmpedge.reserve(_m);
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
        // cutsets.emplace_back(ao.cutset(d.get()));
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
        for (int i(0), cnt(0); i < _n; ++i)
        {
            eid[i] = cnt;
            cnt += deg[i];
        }
        copy(eid.get(), eid.get() + _n, ceid.get());
        for (auto &[u, v, w] : tmpedge)
        {
            edg[ceid[u]++] = cached_edge(v, w, ceid[v]);
            edg[ceid[v]++] = cached_edge(u, w, ceid[u] - 1);
        }
        eid[_n] = edg.size();
        tmpedge.clear();

        ahuja_orlin::max_flow_instance<int, U> ao(_n, deg.get(), eid.get(), edg);
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
    auto sol = Solver<U>(n, m);
    while (m--)
    {
        int u, v;
        U w;
        cin >> u >> v >> w;
        sol.add_edge(u - 1, v - 1, w);
    }
    sol.solve();
    cin >> n;
    while (n--)
    {
        int x, y;
        cin >> x >> y;
        cout << sol.ans[x - 1][y - 1] << '\n';
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
    auto sol = Solver<U>(n, m);
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