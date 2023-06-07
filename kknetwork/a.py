import networkx as nx
import numpy as np
import matplotlib.pyplot as plt


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
    plt.show()

    # get I matrix
    # I = np.array([[0] * len(G)] * len(G))
    # 引入更改
    ## 先删后加，中间得交集
    G.remove_node('F')
    V_union = list(G.nodes()) 
    G.add_edge('A', 'G')
    V_new = ['G']
    E_new = [('A', 'G')]

    dirty = True

    def getI(G):
        nl, mp = refresh_nodelist(G)
        I = nx.adjacency_matrix(G, nl).todense()
        I = np.array(I, dtype=np.float64)
        # print(I)
        for nodeid, nodename in enumerate(nodelist):
            # for pos in range(len(I)):
            I[:][nodeid] /= len(G[nodename])
        # print(I)
        # II = np.copy(I)
        # for i in range(1000):
        #     II = np.matmul(II, I)
        #     print(II)
        return I

    def getP0(G, V_union, pos, V_new, E_new):
        I = getI(G)
        P_pos = np.zeros((len(G), ), dtype=np.float64)
        P0 = np.zeros((len(G), ), dtype=np.float64)
        for i in V_union:
            mi = mapping[i]
            for j in G[i]:
                mj = mapping[j]
                distance_l2 = np.linalg.norm(pos[i] - pos[j]) # l2范数距离（欧氏距离
                l1 = 1 # 理想距离，即最短路，必定为1，因为是邻接的
                P_pos[mi] += I[mi][mj] * abs(distance_l2 - l1) / l1
        # 新点置1
        for i in V_new:
            mi = mapping[i]
            P0[mi] = 1
        
        for i, j in E_new:
            mi, mj = mapping[i], mapping[j]
            if i in V_union: # 暴力
                P0[mi] = I[mi][mj] + P_pos[mi]

        return P0

    I = getI(G)
    P0 = getP0(G, V_union, pos, V_new, E_new)
    PP = np.copy(P0)

    print(P0)
    for _ in range(100):
        PP = np.matmul(PP, I)
        print(PP)


