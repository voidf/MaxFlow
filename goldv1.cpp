// #include <bits/stdc++.h>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <cassert>
#include <queue>
#include <numeric>
#include <list>
#include <cmath>
#include <chrono>

using namespace std;
#define vec vector
#define mp unordered_map
#define us unordered_set


struct GoldbergRao
{
    using T = int;
    const T INF = ~0u >> 1;
    struct Edge
    {
        T capacity;
        T flow = 0;
        T len = -1;
        T _len = -1;
        Edge(T cap) : capacity(cap), flow(0) {}
        Edge() : capacity(0) {}
        void init() { flow = 0; }
    };

    struct Graph
    {
        vec<mp<int, Edge>> E; // 需要查给定uv边，这里暂时用哈希表存
        Graph(int n) : E(n) {}
        void add_edge(int u, int v, T cap) // undirected
        {
            E[u].emplace(v, cap);
            E[v].emplace(u, cap);
        }
    };

    struct ContractedGraph : Graph
    {
        vec<vec<int>> pred;    // 反图
        vec<vec<int>> members; // 有单点查in操作
        vec<int> distances;
        vec<char> blocked;       // bool数组
        vec<T> excess;
        // 边
        vec<mp<int, vec<pair<int, int>>>> emembers;
        // 图变量
        int start_mapping, end_mapping;
        ContractedGraph(int n) : Graph(n),
                                 members(n),
                                 distances(n, 0),
                                 blocked(n, false),
                                 emembers(n, mp<int, vec<pair<int, int>>>()),
                                 pred(n, vec<int>()),
                                 excess(n, 0) {}
        mp<int, Edge>::iterator add_edge_directed(int u, int v, T cap)
        {
            pred[v].emplace_back(u); // 有重边的话这里要改set
            auto [itr, ok] = E[u].emplace(v, cap);
            return itr;
        }
    };

    Graph R;
    int N, M;
    int s, t;
    vec<T> dist;
    vec<T> nflow;
    vec<vec<int>> out_children, in_children;
    // int ctr_bf = 0;
    // int ctr_translate = 0;
    // int ctr_li = 0;
    // long long ctr_scc = 0;
    // long long ctr_sccedge = 0;
    // double t_cond = 0;
    // double t_bf = 0;
    // double t_mc = 0;
    // double t_dm = 0;
    // double t_tr = 0;
    // double t_tarjan = 0;

    GoldbergRao(int n, int m) : N(n), M(m), R(n)
    {
        dist.assign(N, 0);
        out_children.assign(N, {});
        in_children.assign(N, {});
    }

    void rec()
    {
        for (auto &x : R.E)
        {
            for (auto &[to, ed] : x)
            {
                ed.flow = 0;
            }
        }
        // ctr_bf = 0;
        // ctr_translate = 0;
        // ctr_li = 0;
        // ctr_scc = 0;
        // ctr_sccedge = 0;

        // t_cond = 0;
        // t_bf = 0;
        // t_mc = 0;
        // t_dm = 0;
        // t_tr = 0;
        // t_tarjan = 0;
    }

    static inline T iceil(T a, T b) { return a / b + (a % b > 0); }

    // const
    bool is_at_capacity(int u, int v) const { return get_residual_cap(u, v) == 0; }

    T get_residual_cap(int u, int v, bool include_reverse_flow = true) const
    {
        T val = R.E[u].at(v).capacity - R.E[u].at(v).flow;
        if (include_reverse_flow)
            val += R.E[v].at(u).flow;
        return val;
    }
    // 瓶颈2
    void construct_distance_metric(int t) // 只找到了_length的引用，所以这里不判断了
    {
        // auto begintime = chrono::steady_clock::now();
        dist.assign(N, INF);
        dist[t] = 0;
        vec<us<int>> buckets(N);
        int bucket_idx = 0;
        buckets[0].emplace(t);
        while (1)
        {
            while (bucket_idx < N && buckets[bucket_idx].size() == 0)
                ++bucket_idx;
            if (bucket_idx == N)
                break;
            int vertex = *buckets[bucket_idx].begin();
            buckets[bucket_idx].erase(buckets[bucket_idx].begin());
            for (auto &[neighbor, _] : R.E[vertex]) // 这里是遍历入边，但是建边的时候已经建了反边，也能一样
            {
                if (is_at_capacity(neighbor, vertex))
                    continue;
                T length_neighbor = R.E[neighbor][vertex]._len;
                T dist_vertex = dist[vertex];
                T dist_neighbor = dist[neighbor];
                if (dist_neighbor == INF || dist_neighbor > dist_vertex + length_neighbor)
                {
                    if (dist_neighbor != INF)
                        buckets[dist_neighbor].erase(neighbor);
                    dist_neighbor = dist_vertex + length_neighbor;
                    dist[neighbor] = dist_neighbor;
                    buckets[dist_neighbor].emplace(neighbor);
                }
            }
        }
        // t_dm += chrono::duration_cast<chrono::nanoseconds>((chrono::steady_clock::now() - begintime)).count();
    }
    // O(m)
    T min_canonical_cut(int start_node) const
    {
        // auto begintime = chrono::steady_clock::now();

        T max_distance = dist[start_node];
        if (max_distance == 0)
            return INF;
        vec<T> dist_arr(max_distance, 0);
        for (int u = 0; u < N; ++u)
            for (const auto &[v, a] : R.E[u])
                if (is_at_capacity(u, v) == 0 && dist[u] == dist[v] + 1 && dist[v] < max_distance)
                    dist_arr[dist[v]] += get_residual_cap(u, v);
        // t_mc += chrono::duration_cast<chrono::nanoseconds>((chrono::steady_clock::now() - begintime)).count();
        return *min_element(dist_arr.begin(), dist_arr.end());
    }

    bool is_admissible_edge(int u, int v) const { return dist[u] == dist[v] + R.E[u].at(v).len; }

    ContractedGraph condensation()
    {
        // auto begintime = chrono::steady_clock::now();

        vec<vec<int>> scc;
        // ++ctr_scc;
        int st_mp, ed_mp;
        // int si;
        // scc
        vec<int> dfn(N, 0);
        vec<int> dfrom(N, INF); // 用于处理上一个dfs过的儿子
        vec<int> low(N, INF);
        vec<int> belongs(N, 0);           // 已经被缩点的才打1
        vec<mp<int, Edge>::iterator> cur; // 非递归scc用当前弧，可以干掉看看是不是更快
        vec<int> sta;
        vec<int> dmp; // 成员栈，scc还是得开两个栈
        dmp.reserve(N);
        sta.reserve(N);
        cur.reserve(N);
        int cbelongs = 0; // 这个不一定是最终个数，得到的belongs数组相当于一个并查集，还得拉直才是强连通分量

        for (int i = 0; i < N; ++i)
            cur.emplace_back(R.E[i].begin());

        int dtime = 0;
        for (int i = 0; i < N; ++i)
        {
            if (belongs[i])
                continue;
            // assert(sta.size() == 0);
            sta.emplace_back(i);
            // si = 0;
            dfn[i] = low[i] = ++dtime;
            while (sta.size())
            {
                // int v = sta[si];
                int v = sta.back();
                bool done = 1;
                while (cur[v] != R.E[v].end())
                {
                    // ++ctr_sccedge;
                    auto [w, a] = *cur[v];
                    if (dfrom[v] != INF)
                    {
                        low[v] = min(low[dfrom[v]], low[v]);
                        dfrom[v] = INF; // 确保只消费一次
                    }
                    ++cur[v];
                    // if (dfrom[v] == w)
                    // continue;
                    if (!is_at_capacity(v, w) && a.len == 0 && is_admissible_edge(v, w))
                    {
                        if (dfn[w])
                        {
                            if (!belongs[w])
                                low[v] = min(low[v], dfn[w]);
                        }
                        else
                        {
                            sta.emplace_back(w);
                            // ++si;
                            done = 0;
                            low[w] = dfn[w] = ++dtime;
                            dfrom[v] = w;
                            break;
                        }
                    }
                }
                if (done)
                {
                    if (dfrom[v] != INF)
                    {
                        low[v] = min(low[dfrom[v]], low[v]);
                        dfrom[v] = INF;
                    }
                    dmp.emplace_back(sta.back());
                    sta.pop_back();

                    // --si;
                    if (low[v] == dfn[v])
                    {
                        scc.emplace_back();
                        // ++cbelongs;
                        while (dmp.size() && dfn[dmp.back()] >= dfn[v])
                        {
                            int c = dmp.back();
                            dmp.pop_back();
                            scc.back().emplace_back(c);
                            belongs[c] = scc.size();
                            if (c == this->s)
                                st_mp = scc.size() - 1;
                            else if (c == this->t)
                                ed_mp = scc.size() - 1;
                            // if (c == v)
                            // break;
                        }
                    }
                }
            }
        }
        // t_tarjan += chrono::duration_cast<chrono::nanoseconds>((chrono::steady_clock::now() - begintime)).count();

        // auto scc = length_strongly_connected_components();
        if (R.E.size() == 0)
            return ContractedGraph(0);
        int i = 0;
        vec<int> mapping(N, -1); /**/
        int number_of_components = scc.size();

        ContractedGraph C(number_of_components); // 这里是有向图
        C.start_mapping = st_mp;
        C.end_mapping = ed_mp;
        for (auto &component : scc)
        {
            // C.members[i].insert(component.begin(), component.end());
            C.members[i] = component;
            C.distances[i] = dist[*component.begin()];
            for (auto n : component)
                mapping[n] = i;
            ++i;
        }
        for (int u = 0; u < N; ++u)
            for (auto &[v, a] : R.E[u]) // O(m)
            {
                int mu = mapping[u], mv = mapping[v];
                if (mu != mv && !is_at_capacity(u, v) && is_admissible_edge(u, v))
                {
                    auto it = C.emembers[mu].find(mv);
                    if (it != C.emembers[mu].end())
                    {
                        it->second.emplace_back(u, v);
                        C.E[mu][mv].capacity += get_residual_cap(u, v);
                        // assert(R.E[u][v].len == C.E[mu][mv].len);
                    }
                    else
                    {
                        auto em = C.add_edge_directed(mu, mv, get_residual_cap(u, v));
                        // em->second.flow = 0; // 其实没必要写
                        em->second.len = R.E[u][v].len;
                        C.emembers[mu][mv].emplace_back(u, v);
                    }
                }
            }

        return C;
    }

    ContractedGraph construct_graph_contraction(int start_node, int end_node)
    {
        // auto begintime = chrono::steady_clock::now();

        auto condensed_graph = condensation();
        for (int scc = 0; scc < condensed_graph.E.size(); ++scc)
        {
            auto &scc_members = condensed_graph.members[scc];
            int rep_vertex = *scc_members.begin();
            // if (scc_members.count(start_node))
            //     condensed_graph.start_mapping = scc;
            // if (scc_members.count(end_node))
            //     condensed_graph.end_mapping = scc;
            if (scc_members.size() < 2)
                continue;

            // out tree
            vec<int> children_queue({rep_vertex});
            list<int> not_visited; // 魔改set的写法
            for (auto i : scc_members)
                if (i != rep_vertex)
                    not_visited.emplace_back(i);
            while (children_queue.size())
            {
                int curr_vertex = children_queue.back();
                children_queue.pop_back();
                auto ni = not_visited.begin();
                while (ni != not_visited.end())
                {
                    auto neighbor = *ni;
                    if (R.E[curr_vertex].count(neighbor) == 0 || is_at_capacity(curr_vertex, neighbor))
                        ;
                    else if (R.E[curr_vertex][neighbor].len == 0)
                    {
                        children_queue.emplace_back(neighbor);
                        out_children[curr_vertex].emplace_back(neighbor);
                        auto bak = ni;
                        ++ni;
                        not_visited.erase(bak);
                        continue; // 避免重复移动ni
                    }
                    ++ni;
                }
            }
            // in tree
            children_queue.clear();
            children_queue.emplace_back(rep_vertex);
            not_visited.clear();
            for (auto i : scc_members)
                if (i != rep_vertex)
                    not_visited.emplace_back(i);
            while (children_queue.size())
            {
                int curr_vertex = children_queue.back();
                children_queue.pop_back();
                auto ni = not_visited.begin();
                while (ni != not_visited.end())
                {
                    auto neighbor = *ni;
                    if (R.E[neighbor].count(curr_vertex) == 0 || is_at_capacity(neighbor, curr_vertex))
                        ;
                    else if (R.E[neighbor][curr_vertex].len == 0)
                    {
                        children_queue.emplace_back(neighbor);
                        in_children[curr_vertex].emplace_back(neighbor);
                        auto bak = ni;
                        ++ni;
                        not_visited.erase(bak);
                        continue;
                    }
                    ++ni;
                }
            }
        }
        // t_cond += chrono::duration_cast<chrono::nanoseconds>((chrono::steady_clock::now() - begintime)).count();

        return condensed_graph;
    }

    T flow_value(ContractedGraph &C, int start_node, int end_node)
    {
        T X = 0;
        for (auto u : C.pred[end_node])
            X += C.E[u][end_node].flow;
        for (auto &[w, a] : C.E[end_node])
            X -= a.flow;
        return X;
    }

    list<int> topological_sort(const vec<mp<int, Edge>> &E)
    {
        vec<int> indegree(E.size(), 0);
        list<int> L;   // 拓扑序点的标号
        queue<int> _q; // 暂存标号
        for (int u = 0; u < E.size(); ++u)
            for (auto &[v, a] : E[u])
                ++indegree[v];
        for (int u = 0; u < E.size(); ++u)
            if (indegree[u] == 0)
                _q.emplace(u);
        while (_q.size())
        {
            int u = _q.front();
            L.emplace_back(u);
            _q.pop();
            for (auto &[v, a] : E[u])
                if (--indegree[v] == 0)
                    _q.emplace(v);
        }
        // assert(L.size() == E.size());
        return L;
    }

    T limit_flow(ContractedGraph &C, int start_node, int end_node, T maximum_flow_to_route)
    {
        auto X = flow_value(C, start_node, end_node);
        if (X <= maximum_flow_to_route)
            return X;
        C.excess[end_node] = X - maximum_flow_to_route;
        auto topo = topological_sort(C.E);
        for (auto vp = topo.rbegin(); vp != topo.rend(); ++vp)
        {
            int v = *vp;
            for (auto u : C.pred[v])
            {
                if (C.excess[v] == 0)
                    break;
                T delta = min(C.E[u][v].flow, C.excess[v]);
                C.excess[v] -= delta;
                C.excess[u] += delta;
                C.E[u][v].flow -= delta;
            }
        }
        return maximum_flow_to_route;
    }

    T compute_blocking_flow(ContractedGraph &C, int start_node, int end_node, T maximum_flow_to_route)
    {
        // auto begintime = chrono::steady_clock::now();
        // ++ctr_bf;
        // list<int> actives;

        C.blocked.assign(C.E.size(), false);
        C.excess.assign(C.E.size(), 0);
        for (int u = 0; u < C.E.size(); ++u)
            for (auto &[v, a] : C.E[u])
            {
                if (u == start_node)
                    C.excess[v] = a.flow = a.capacity;
                else
                    a.flow = 0;
            }
        auto push = [&](int u, int v)
        {
            T delta = min(C.excess[u], C.E[u][v].capacity - C.E[u][v].flow);
            C.E[u][v].flow += delta; // 这里没给反边加流量，有点怪
            C.excess[v] += delta;
            // if (C.excess[v] > 0 && v != start_node && v != end_node)
                // actives.emplace_back(v);
            C.excess[u] -= delta;
        };
        auto pull = [&](int u, int v)
        {
            T delta = min(C.excess[v], C.E[u][v].flow);
            C.E[u][v].flow -= delta;
            C.excess[v] -= delta;
            C.excess[u] += delta;
            // if (C.excess[u] > 0 && u != start_node && u != end_node)
                // actives.emplace_back(u);
        };
        auto discharge = [&](int v)
        {
            if (!C.blocked[v])
                for (auto &[w, a] : C.E[v])
                    if (!C.blocked[w])
                    {
                        push(v, w); // v -> w
                        if (C.excess[v] == 0)
                            return;
                    }
            C.blocked[v] = true;
            for (auto u : C.pred[v])
            {
                pull(u, v); // u <- v
                if (C.excess[v] == 0)
                    return;
            }
            throw "Unexpected discharge behaviour";
        };
        // toposort
        auto L = topological_sort(C.E);
        // end toposort
        // actives.reserve(C.E.size());
        
        // for (auto i : L)
        // for (int i=0;i<C.E.size();++i)
        //     if (i != start_node && i != end_node && C.excess[i] > 0)
        //         actives.emplace_back(i);

        auto first_active = [&]() -> list<int>::iterator
        {
            auto vp = L.begin();
            while (vp != L.end())
            {
                if (*vp != start_node && *vp != end_node && C.excess[*vp] > 0)
                    return vp;
                ++vp;
                // ++ctr_li;
            }
            return vp;
        };
        auto vp = first_active();
        while (vp != L.end())
        // while (actives.size())
        {
            // ++ctr_li;
            // auto vv = actives.back();
            // actives.pop_back();
            // discharge(vv);
            discharge(*vp);
            if (C.blocked[*vp])
            {
                L.emplace_front(*vp);
                L.erase(vp);
            }
            vp = first_active();
        }
        // t_bf += chrono::duration_cast<chrono::nanoseconds>((chrono::steady_clock::now() - begintime)).count();

        return limit_flow(C, start_node, end_node, maximum_flow_to_route);
    }

    void update_flow(int u, int v, T flow_val)
    {
        Edge &a = R.E[u][v];
        Edge &r = R.E[v][u];
        T new_flow = flow_val + a.flow - r.flow;
        if (new_flow < 0)
        {
            // assert(-new_flow <= r.capacity);
            a.flow = 0;
            r.flow = -new_flow;
        }
        else
        {
            // assert(new_flow <= a.capacity);
            a.flow = new_flow;
            r.flow = 0;
        }
    }

    T route_in_flow_tree(int curr_vertex)
    {
        T flow = max(nflow[curr_vertex], 0);
        for (auto child : in_children[curr_vertex])
        {
            T child_flow = route_in_flow_tree(child);
            if (child_flow != 0)
                update_flow(child, curr_vertex, child_flow);
            flow += child_flow;
        }
        in_children[curr_vertex].clear();
        return flow;
    }

    T route_out_flow_tree(int curr_vertex)
    {
        T flow = max(-nflow[curr_vertex], 0);
        for (auto child : out_children[curr_vertex])
        {
            T child_flow = route_out_flow_tree(child);
            if (child_flow != 0)
                update_flow(curr_vertex, child, child_flow);
            flow += child_flow;
        }
        out_children[curr_vertex].clear();
        nflow[curr_vertex] = 0;
        return flow;
    }

    void translate_flow_from_contraction_to_original(ContractedGraph &C, const int start_node, const int end_node, const T flow_routed)
    {
        // auto begintime = chrono::steady_clock::now();
        // ++ctr_translate;
        nflow.assign(N, 0); // original
        nflow[start_node] = flow_routed;
        nflow[end_node] = -flow_routed;
        for (int cu = 0; cu < C.E.size(); ++cu)
            for (auto &[cv, ca] : C.E[cu])
            {
                T remaining_edge_flow = ca.flow;
                if (remaining_edge_flow == 0)
                    continue;
                for (auto [start_vert, end_vert] : C.emembers[cu][cv])
                {
                    T flow_to_route = min(remaining_edge_flow, get_residual_cap(start_vert, end_vert));
                    update_flow(start_vert, end_vert, flow_to_route);
                    nflow[start_vert] -= flow_to_route; /**/
                    nflow[end_vert] += flow_to_route;
                    if (flow_to_route == remaining_edge_flow)
                        break;
                    remaining_edge_flow -= flow_to_route;
                }
            }
        for (int v = 0; v < C.E.size(); ++v)
            if (C.members[v].size() >= 2)
            {
                int rep = *C.members[v].begin();
                T flow_in = route_in_flow_tree(rep);
                T flow_out = route_out_flow_tree(rep);
                // assert(flow_in == flow_out);
            }
        // t_tr += chrono::duration_cast<chrono::nanoseconds>((chrono::steady_clock::now() - begintime)).count();
        
    }

    vec<int> stcut()
    {
        vec<int> D(N, 0);
        vec<int> Q(N, 0);
        Q[0] = s;
        D[s] = 1;
        int top = 0;
        int cap = 1;
        while (top < cap)
        {
            int v = Q[top++];
            for (auto &[to, e] : R.E[v])
            {
                if (!D[to] && e.flow < e.capacity)
                {
                    D[to] = D[v] + 1;
                    // assert(to != t);
                    Q[cap++] = to;
                }
            }
        }
        return D;
    }

    T solve(const int s, const int t)
    {
        this->s = s;
        this->t = t;

        T sum_capacity = 0;
        for (auto &ee : R.E)
            for (auto &[v, a] : ee)
                sum_capacity += a.capacity;
        T error_bound = sum_capacity;
        /**/ int num_iterations_in_phase = ceil(min(sqrt(M), pow(N, 2. / 3.)));
        T total_routed_flow = 0;
        T prev_error_bound = INF;
        while (error_bound >= 1) // log2(sum(U))
        {
            // assert(prev_error_bound == INF || error_bound <= prev_error_bound / 2);
            /**/ T flow_to_route = iceil(error_bound, num_iterations_in_phase);
            for (int _ = 0; _ < 8 * num_iterations_in_phase; ++_) // min(sqrt(M), pow(N, 2. / 3.)
            {
                for (int u = 0; u < N; ++u)
                    for (auto &[v, a] : R.E[u]) // O(m)
                    {
                        if (is_at_capacity(u, v))
                        {
                            a.len = a._len = -1;
                            continue;
                        }
                        a._len = get_residual_cap(u, v) >= flow_to_route * 3 ? 0 : 1;
                    }
                construct_distance_metric(t);
                if (dist[s] == INF)
                    return total_routed_flow;
                T max_flow_upper_bound = min_canonical_cut(s);
                if (max_flow_upper_bound <= error_bound / 2)
                {
                    prev_error_bound = error_bound;
                    while (max_flow_upper_bound <= error_bound / 2 && error_bound >= 1)
                        error_bound /= 2;
                    break;
                }
                for (int u = 0; u < N; ++u)
                    for (auto &[v, a] : R.E[u])
                    {
                        if (is_at_capacity(u, v))
                            continue;
                        T resid_cap = get_residual_cap(u, v);
                        T resid_cap_reverse = get_residual_cap(v, u);
                        if (2 * flow_to_route <= resid_cap &&
                            resid_cap < 3 * flow_to_route &&
                            flow_to_route <= resid_cap_reverse &&
                            dist[u] == dist[v])
                            a.len = 0;
                        else
                            a.len = a._len;
                    }
                // 这个算法实际耗时瓶颈是scc
                auto contracted_graph = construct_graph_contraction(s, t);
                T flow_routed = flow_to_route;
                if (contracted_graph.start_mapping != contracted_graph.end_mapping)
                    flow_routed = compute_blocking_flow(contracted_graph, contracted_graph.start_mapping, contracted_graph.end_mapping, flow_to_route);
                total_routed_flow += flow_routed;
                if (flow_routed == 0)
                    return total_routed_flow;
                translate_flow_from_contraction_to_original(contracted_graph, s, t, flow_routed);
            }
        }
        throw "23333";
        // return total_routed_flow;
    }
} * gr;

// vec<int> I{-1, 0, 1, 2, 3, 4, 5, 6, 7},
// J{1, 2, 3, 4, 5, 6, 7, 8};
vec<int> I, J;

void with_clean_data(int n, int m)
{
    using T = int;
    mp<int, T> P;
    us<int> S;
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
    I.assign(n + 1, -1);
    J.reserve(S.size());
    for (auto k : S)
    {
        I[k] = ctr++;
        J.emplace_back(k);
    }
    int nn = S.size();
    int mm = P.size();

    gr = new GoldbergRao(nn, mm);

    for (auto [k, w] : P)
    {
        int u = k / (n + 1);
        int v = k % (n + 1);
        gr->R.add_edge(I[u], I[v], w);
    }
}

vec<int> node, tmp1, tmp2;
vec<vec<int>> ans;
vec<int> d;
int n, m;

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
        gr->rec();
        t = gr->solve(I[S], I[T]);
        // cerr << "ctr_bf:" << gr->ctr_bf << " ";
        // cerr << "ctr_translate:" << gr->ctr_translate << " ";
        // cerr << "ctr_li:" << gr->ctr_li << " ";
        // cerr << "ctr_scc:" << gr->ctr_scc << " ";
        // cerr << "ctr_sccedge:" << gr->ctr_sccedge << " ";
        // cerr << endl;

        // cerr << "blockingflow:" << gr->t_bf << " ";
        // cerr << "translate:" << gr->t_tr << " ";
        // cerr << "scc:" << gr->t_cond << " ";
        // cerr << "distance_metric:" << gr->t_dm << " ";
        // cerr << "mini cut:" << gr->t_mc << " ";
        // cerr << "tarjan:" << gr->t_tarjan << " ";
        // cerr << endl;

        // hp->returnExcess();
        for (auto i : gr->stcut())
            d[J[c++]] = i;
        // hp->recoverUndirected();
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