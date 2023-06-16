from itertools import chain
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
    def modify_graph(
        deleted_nodes,
        deleted_edges,
        inserted_nodes,
        inserted_edges,
    ):
        global dirty
        for i in deleted_nodes:
            G.remove_node(i)
        for i, j in deleted_edges:
            G.remove_edge(i, j)
        V_union = list(G.nodes()) 
        for i in inserted_nodes:
            G.add_node(i)
        for i, j in inserted_edges:
            G.add_edge(i, j)
        dirty = True
        # pos = nx.kamada_kawai_layout(G)
        # nx.draw_networkx(G, pos)
        # plt.show()
        return V_union

    V_union = modify_graph(deleted_nodes, deleted_edges, inserted_nodes, inserted_edges)

    def getI(G):
        nl, mp = refresh_nodelist(G)
        I = nx.adjacency_matrix(G, nl).todense()
        I = np.array(I, dtype=np.float64)
        # print(I)
        for nodeid, nodename in enumerate(nodelist):
            # for pos in range(len(I)):
            if len(G[nodename]): # 避免除零
                I[:][nodeid] /= len(G[nodename])
        # print(I)
        # II = np.copy(I)
        # for i in range(1000):
        #     II = np.matmul(II, I)
        #     print(II)
        return I # 节点影响力矩阵

    def getP0(G, V_union, pos, V_new, E_new, E_remv):
        I = getI(G)
        P_pos = np.zeros((len(G), ), dtype=np.float64)
        for i in V_union:
            mi = mapping[i]
            for j in G[i]:
                if j in V_union:
                    mj = mapping[j]
                    distance_l2 = np.linalg.norm(pos[i] - pos[j]) # l2范数距离（欧氏距离
                    l1 = 1 # 理想距离，即最短路，必定为1，因为是邻接的
                    P_pos[mi] += I[mi][mj] * abs(distance_l2 - l1) / l1
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
                P_move[mi] = I[mi][mj] + P_pos[mi]
                # P0[mi] += I[mi][mj]
            if j in V_union:
                P_move[mj] = I[mj][mi] + P_pos[mj] # 无向不带权图，I[mj]整列非零格子都是度数的倒数，不用考虑u点遵从哪个的问题
                # P0[mj] += I[mj][mi]

        for i in set(chain(*E_remv)):
            mi = mapping[i]
            if i in V_union: # 公式3-13
                if P_move[mi] != 0: # 增加的边优先
                    continue
                # P0[mi] = I[mi][mj] + P_pos[mi]
                P_move[mi] = P_pos[mi]

        return P_move

    I = getI(G)
    print(I)
    
    P0 = getP0(G, V_union, pos, inserted_nodes, inserted_edges, deleted_edges)
    # P0 = np.array([0,0,0,0,0,1,1], dtype=np.float64)
    PP = np.copy(P0)
    print(P0)


    # for _ in range(100):
    #     PP = np.matmul(PP, I)
    #     print(PP)

    # print(PP)
    print()
    P_move_final = np.matmul(P0, np.linalg.matrix_power(I, 100)) # 没有收敛就调大这个100，矩阵快速幂，非常快，但是I的幂本身是不收敛的，要小心溢出变成nan
    print(P_move_final)
    gamma = 1 # 我觉得这个γ是某篇引文里面的超参，全文只有一次出现，非常神秘
    Si_hat = gamma / (P_move_final + 0.1)
    beta_2 = 1 / max(Si_hat) # β2是一个用来归一化的常数
    Si = Si_hat * beta_2
    print(Si)
