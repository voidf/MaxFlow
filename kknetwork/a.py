import copy
from itertools import chain
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

from networkx.drawing.layout import _process_params, random_layout, rescale_layout

def _kamada_kawai_costfn(pos_vec, np, invdist, meanweight, dim):
    # Cost-function and gradient for Kamada-Kawai layout algorithm
    # invdist: 1 / 距离矩阵
    nNodes = invdist.shape[0]
    pos_arr = pos_vec.reshape((nNodes, dim))

    delta = pos_arr[:, np.newaxis, :] - pos_arr[np.newaxis, :, :] # 本节点和其它节点的位置坐标之差 shape=(n, n, 2)
    nodesep = np.linalg.norm(delta, axis=-1) # 本节点和其它节点的欧氏距离 shape=(n, n)
    direction = np.einsum("ijk,ij->ijk", delta, 1 / (nodesep + np.eye(nNodes) * 1e-3)) # 变化方向，除掉之后每个delta模长为1 shape=(n, n, 2)

    offset = nodesep * invdist - 1.0
    offset[np.diag_indices(nNodes)] = 0

    cost = 0.5 * np.sum(offset**2)
    grad = np.einsum("ij,ij,ijk->ik", invdist, offset, direction) - np.einsum(
        "ij,ij,ijk->jk", invdist, offset, direction
    )

    # Additional parabolic term to encourage mean position to be near origin:
    sumpos = np.sum(pos_arr, axis=0) # shape=(2)， 分别为x，y方向的位置坐标的总和
    cost += 0.5 * meanweight * np.sum(sumpos**2)
    grad += meanweight * sumpos

    return (cost, grad.ravel())


def cost_fn(pos_vec, np, invdist, meanweight, dim):
    cost, grad = _kamada_kawai_costfn(pos_vec, np, invdist, meanweight, dim)
    

# 魔改一下kk
def kamada_kawai_layout(
    G, dist=None, weight="weight", scale=1, center=None, dim=2
):
    import numpy as np

    G, center = _process_params(G, center, dim)
    nNodes = len(G)
    if nNodes == 0:
        return {}

    if dist is None:
        dist = dict(nx.shortest_path_length(G, weight=weight))
    dist_mtx = 1e6 * np.ones((nNodes, nNodes))
    for row, nr in enumerate(G):
        if nr not in dist:
            continue
        rdist = dist[nr]
        for col, nc in enumerate(G):
            if nc not in rdist:
                continue
            dist_mtx[row][col] = rdist[nc]

    assert pos is not None
    pos_arr = np.array([pos[n] for n in G])

    pos = _kamada_kawai_solve(dist_mtx, pos_arr, dim)

    pos = rescale_layout(pos, scale=scale) + center
    return dict(zip(G, pos))


def _kamada_kawai_solve(dist_mtx, pos_arr, dim):
    # pos_arr: {"A": (114514, 1), "B": (1919810, 0)} 按每个节点给的位置信息
    # dist_mtx: [[0, 114514], [1919810, 1]] n^2的距离矩阵
    # dim: 我们的case中始终等于2，不用管
    # Anneal node locations based on the Kamada-Kawai cost-function,
    # using the supplied matrix of preferred inter-node distances,
    # and starting locations.

    import numpy as np
    import scipy as sp
    import scipy.optimize  # call as sp.optimize

    meanwt = 1e-3
    costargs = (np, 1 / (dist_mtx + np.eye(dist_mtx.shape[0]) * 1e-3), meanwt, dim)

    optresult = sp.optimize.minimize(
        _kamada_kawai_costfn,
        pos_arr.ravel(),
        method="L-BFGS-B",
        args=costargs,
        jac=True,
    )

    return optresult.x.reshape((-1, dim))


if __name__ == "__main__":
    dirty = False
    nodelist = ['A', 'B', 'C', 'D', 'E', 'F']

    def refresh_nodelist(G: nx.Graph):
        global nodelist, mapping, dirty
        if dirty:
            nodelist = list(sorted(G.nodes()))
            mapping = {k: v for v, k in enumerate(nodelist)}
            dirty = False
        return nodelist, mapping

    G = nx.Graph([
        ['A', 'B'],
        ['A', 'C'],
        ['A', 'D'],
        ['C', 'G'],
        ['D', 'E'],
        ['D', 'F'],
        ['E', 'F'],
    ])

    pos = nx.kamada_kawai_layout(G)
    print(pos)
    # nx.draw_networkx(G, pos)
    # plt.show()

    # get I matrix
    # I = np.array([[0] * len(G)] * len(G))
    # 引入更改
    ## 先删后加，中间得交集

    # deleted_nodes = ['E']
    # deleted_edges = [('F', 'D')]
    # inserted_nodes = ['H']
    # inserted_edges = [('A', 'G'), ('A', 'H')]
    
    deleted_nodes = []
    deleted_edges = [('F', 'D')]
    inserted_nodes = []
    inserted_edges = [('F', 'G')]
    # def modify_graph(
    #     deleted_nodes,
    #     deleted_edges,
    #     inserted_nodes,
    #     inserted_edges,
    # ):
    # global dirty
    for i in deleted_nodes:
        G.remove_node(i)
    for i, j in deleted_edges:
        G.remove_edge(i, j)
    G_union = copy.deepcopy(G)
    V_union = list(sorted(G.nodes())) 
    for i in inserted_nodes:
        G.add_node(i)
    for i, j in inserted_edges:
        G.add_edge(i, j)
    dirty = True
    # pos = nx.kamada_kawai_layout(G)
    # nx.draw_networkx(G, pos)
    # plt.show()
    # return V_union

    # V_union = modify_graph(deleted_nodes, deleted_edges, inserted_nodes, inserted_edges)

    def getI(G):
        nl, mp = refresh_nodelist(G)
        I = nx.adjacency_matrix(G, nl).todense()
        I = np.array(I, dtype=np.float64)
        # print(I)
        for nodeid, nodename in enumerate(nodelist):
            # for pos in range(len(I)):
            if len(G[nodename]): # 避免除零
                I[:, nodeid] /= len(G[nodename]) # I阵完全就是为了幂乘而得出的，平时用直接1 / len(G[i])即可
        # print(I)
        # II = np.copy(I)
        # for i in range(1000):
        #     II = np.matmul(II, I)
        #     print(II)
        return I # 节点影响力矩阵

    def get_P_move(G, I, V_union, pos, V_new, E_new, E_remv):
        P_pos = np.zeros((len(G), ), dtype=np.float64)
        for i in V_union:
            mi = mapping[i]
            for j in G[i]:
                if j in V_union:
                    mj = mapping[j]
                    distance_l2 = np.linalg.norm(pos[i] - pos[j]) # l2范数距离（欧氏距离
                    l1 = 1 # 理想距离，即最短路，必定为1，因为是邻接的
                    P_pos[mi] += 1 / len(G[j]) * abs(distance_l2 - l1) / l1
        P_move = np.zeros((len(G), ), dtype=np.float64)
        # P0 = np.copy(P_pos) # 式3-13，但这里并不只是删除边的点
        # 式3-11，新点+1，同时由于新0度点的P0是0，这步中相当于置1
        for i in V_new:
            mi = mapping[i]
            # P_move[mi] += 1
            P_move[mi] = 1 + P_pos[mi]
        
        for i, j in E_new:
            mi, mj = mapping[i], mapping[j]
            if i in V_union: # 公式3-10，实际上当i在V_union中时，i必不属于新点，所以不用考虑重复加的问题
                P_move[mi] = 1 / len(G[i]) + P_pos[mi]
                # P0[mi] += I[mi][mj]
            if j in V_union:
                P_move[mj] = 1 / len(G[j]) + P_pos[mj] # 无向不带权图，I[mj]整列非零格子都是度数的倒数，不用考虑u点遵从哪个的问题
                # P0[mj] += I[mj][mi]

        for i in set(chain(*E_remv)):
            mi = mapping[i]
            if i in V_union: # 公式3-13
                if P_move[mi] != 0: # 增加的边优先
                    continue
                # P0[mi] = I[mi][mj] + P_pos[mi]
                P_move[mi] = P_pos[mi]

        return P_move

    def get_P_imp_map(G, I=None, alpha=0.8):
        """我感觉多此一举的模拟退火mechanism，但是没它又不能算imp"""
        if I is None:
            I = getI(G)
        P_imp_map = np.zeros((len(G), len(G))) # 返回应该是个二维阵
        for node_index, node_name in enumerate(nodelist):
            p0 = np.zeros((len(G), ), dtype=np.float64)
            p0[node_index] = 1
            # P_imp_map[node_index] = np.matmul(p0, np.linalg.matrix_power(I, 200)) # 或者我们直接用矩阵快速幂
            p0 = np.matmul(p0, I)
            for n in range(20): # 我们先暂定迭代个n=20步
                P_imp_map[node_index] += p0 # 式3-9
                p0 = np.matmul(p0, alpha * I)

        return P_imp_map # Imp(v) = sum(P_imp_map[v])


    I = getI(G)
    print(I)
    
    P_move_0 = get_P_move(G, I, V_union, pos, inserted_nodes, inserted_edges, deleted_edges)
    # P0 = np.array([0,0,0,0,0,1,1], dtype=np.float64)
    # PP = np.copy(P_move_0)
    print(P_move_0)

    P_imp_map = get_P_imp_map(G, I)
    print('P_imp_map', P_imp_map)

    P_union_imp_map = get_P_imp_map(G_union)

    # 算S_ij
    S_ij_hat = np.zeros((len(V_union), len(V_union)), dtype=np.float64)
    for nodeidx_i, nodename_i in enumerate(V_union):
        # if nodename in G.nodes(): # union 里的一定在new里
        for nodeidx_j, nodename_j in enumerate(V_union):
            R_ij_union = P_union_imp_map[nodeidx_i][nodeidx_j] + P_union_imp_map[nodeidx_j][nodeidx_i]
            mi, mj = mapping[nodename_i], mapping[nodename_j]
            R_ij_new = P_imp_map[mi][mj] + P_imp_map[mj][mi]
            S_ij_hat[nodeidx_i][nodeidx_j] = min(R_ij_union, R_ij_new)

    S_ij = S_ij_hat / np.max(S_ij_hat)
    print('S_ij', S_ij)

    # for _ in range(100):
    #     PP = np.matmul(PP, I)
    #     print(PP)

    # print(PP)
    print()
    P_move_final = np.matmul(P_move_0, np.linalg.matrix_power(I, 100)) # 没有收敛就调大这个100，矩阵快速幂，非常快，但是I的幂本身是不收敛的，要小心溢出变成nan
    print(P_move_final)
    gamma = 1 # 我觉得这个γ是某篇引文里面的超参，全文只有一次出现，非常神秘
    Si_hat = gamma / (P_move_final + 0.1) + np.sum(P_imp_map, axis=1) # 后面这个项是大论文里的Imp
    beta_2 = 1 / max(Si_hat) # β2是一个用来归一化的常数
    Si = Si_hat * beta_2
    print(Si)
