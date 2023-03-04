# 先说怎么用

DBG宏：定义了就是直接运行g++ aomodule.cpp，标准输入洛谷格式的数据，然后标准输出结果。

如果没有定义，则是当做Python模块，底下那#else一堆都是用来对接Python的，也就是把py的int对象和list对象转成c++格式的，不重要。

为了简便，我从main函数说起。

1. 初始化一个Solver类，把n给它，这一步是为了分配内存。

2. 用add_edge方法给它加边

3. solve方法是入口函数，它会完成全源最小割方案的计算，ans成员存全源最小割值，ansedge存全源最小割方案下标，cutsets存全源最小割具体方案。空间需要O(n*m)，实际上据随机图来说应该远小于它，毕竟你不太可能每次都选中所有的边（

我就是通过直接存下所有的方案从而避免了最小割树的构建，但是最小割树分治的思想和本质上只有O(n)的最小割的结论我是用了。

# Solver

solve函数会构建一个max_flow_instance，这是算法主体，传入邻接表G以便分配内存。

work函数是分治，即最小割树中O(n)的分治步骤，大意是求l与r的最小割，然后把他们根据在割集的哪边拆成两个集合，对应于最小割树上的加一条l与r的边。我们记cnt1为l集合的点数，cnt2为r集合的点数，因为我们可以先一次分治到底，然后把l集合塞到node数组的[0, cnt1]，r集合塞到[cnt1+1, cnt1+cnt2]中，然后注意到最小割树上这条边的影响，我们可以遍历两个集合，用这条边来更新答案。

我们枚举l集合点ii，r集合点jj，则有转移式：

$$
minans_{ii,jj} = \min\{minans_{ii,l}, minans_{l,r}, minans_{r, jj} \}
$$

对这个式子的直观理解，就是你考虑ii和l我们已经算出来它的最小割边（因为递归），jj和r同理，那么ii到jj的最小割边只有可能是ii->l、l->r、r->jj这三种答案其中之一，因为他们非通过l->r边不可，最小割树上有这条边，也是ii到jj的必经之路。

至于时间复杂度。我们考虑一个点，它不论怎么分治和递归，在以上的比较流程中都只会用上式和其它点比较一次，所以比较次数是严格n*(n-1)的，总复杂度n^2。

所以，有了这个框架，只要底层算法能够提供s-t割方案（或者更充分的残量网络），我们就能在O(n)*O(s-t割)+O(n^2)的时间内算出全源最小割，包括选边方案。

# ahuja orlin提出的excess scaling algorithm

论文中翻：https://www.bilibili.com/read/cv21334707

本质上是魔改push relabel，借用了超额流的概念。

注意每个节点带一个next和prev，这是为了结合后面的手写链表用，为了表示节点是active的还是inactive的，论文中有提到。

算法是双阶段设计，第一阶段调用find_max_flow算最大流，超额流留在点上，此时给出的残量网络是错的，不能用来最小割树。调用preflow_to_flow后，算法会交换源汇点，把多出来的超额流还给源点，此时的残量网络才是对的，所以在跑最小割树之前要调一下这个。

算法的核心在于find_max_flow_inner函数的那个K，它是一个logU的值，我们`for(int delta=2^K;delta>0;delta/=2)`每次推流的时候，只推不小于delta/2的流。论文中说这么操作有助于让超额流流均等地分布在网络中，而不是集中在某些点上，从而导致那些点反复和周围的点打太极，增加推流次数，影响计算时间。所以这个K就是logU*n^2的加项的复杂度来源。具体证明可以看我的论文翻译。

除了上面说这个操作，其它和push relabel大同小异。

# excess scaling algorithm

## 核心思想

鼓励超额流能够在网络中分布得更均匀，有利于促使中间点向汇点推流更容易。

push relabel推流：一次推完一条边的容量，使残量网络结构发生改变的，叫饱和推流。没推完边容量的叫做非饱和推流。

目前的push relabel推流，复杂度主要受$O(n^2m)$次的非饱和推流制约，即使是最高标号也只是把它优化为$O(n^2\sqrt{m})$次。

每次推流只推不小于delta/2的流，可以降低非饱和推流次数，将其优化为$O(n^2logU)$。从而在流量不大的情况下，复杂度转而变为O(nm)次的饱和推流所制约。

## 输入

邻接表G，源点s，汇点t

我们令G[i][j]表示边i->j的容量。

## 输出

最大流值flow，残量网络R

我们令R[i][j]表示边i->j的剩余容量，也叫做残量。

## 具体步骤

定义n为点数，m为边数，U为最大边容量值。

定义d为距离标号，d[i]衡量i到t的距离，我们的i->j推流只在d[i]==d[j]+1的时候做。

定义excess为超额流数组，excess[i]为i号点上的超额流流量值。

定义n+1个单链表label_bucket，他们的作用是桶，用来装标号分别为0~n的节点有哪些，比如li[0]用来装距离标号为0的节点。它们的加点和删点操作必须是O(1)的。删点操作需要给出一个点马上定位它在哪里。

对于邻接表，我们记每个节点的当前弧current_arc，因为我们希望一次一个点的重标号中每个点的每条邻边只扫一次，这个操作与复杂度相关。（实际实现中因为加了比没加慢，所以没加）

全局重标号
```python
# Goldberg (1987) 建议定期执行全局重标号有利于更快更顺利地推流，实际上确实优化很大，但与复杂度不相关。所以本示例代码中只在预处理中做一次
d = [n for _ in range(n)] # 先令所有点都不可达
d[t] = 0
# 除了源汇点之外，所有点的距离用汇点为起点的反向bfs决定
q = Queue()
q.push((t, 0))
while not q.empty():
    cur, dist = q.pop()
    for to in R[cur] if R[to][cur] == 0 and d[to] == n: # 这里是反向遍历，判断反边上还有没有容量
        d[to] = dist + 1
        q.push((to, d[to]))
```


预处理
```python
# 初始时，令残量网络R等于G
R = deepcopy(G)
# 先把源点邻接的点一次性推流完，把边容量用尽
for to in G[s]:
    flow = G[s][to]
    excess[to] = flow
    R[s][to] = 0 # 把残量网络上这条边的容量用尽
    R[to][s] += to
    
全局重标号()

```

最大流算法
```python
K = 1 + ceil(log2(U))

for x in label_bucket: # 初始化标签桶
    x.clear()

for k in range(K, 0, -1):
    delta = 2 ** k
    for i in range(n) if excess[i] > delta / 2:
        label_bucket[d[i]].add(i)
    level = 1
    while level < 2 * n:
        if len(label_bucket[level]) == 0:
            level += 1
        else:
            i = label_bucket[level].pop()
            推流_重标号操作(i)
```

推流_重标号操作(i)
```python
found = False
j = current_arc[i] # j为i的当前弧所指向的目的点
while not found and j is not None:
    if d[i] == d[j] + 1 and R[i][j] > 0:
        found = True
    else:
        j = current_arc[i] = next_arc(i, j) # 在邻接表中寻找以i为起点的i->j边的下一条边，如果没这条边则返回一个None，这个操作必须是O(1)的，与邻接表的具体实现有关，这里不给出具体代码
    if found:
        对弧(i,j)推min(excess[i], R[i][j], delta - excess[i])的流，同时立刻更新残量网络
        if excess[i] <= delta / 2: # 注意此处是更新后的excess
            label_bucket[d[i]].erase(i)
        if j != s and j != t and excess[j] > delta / 2:
            label_bucket[d[j]].add(j)
            level -= 1
    else:
        label_bucket[d[i]].erase(i)
        d[i] = min(d[j] + 1 for j in R[i] if R[i][j] > 0)
        label_bucket[d[i]].add(i)
        current_arc[i] = first_arc(i) # 找i点邻接表里的第一条边
```

------

入口
```python
预处理()
最大流算法()
maxflow = excess[t] # 最大流值
s, t = t, s # 交换源汇点以便退流
最大流算法() # 这一步将多出来的超额流还给源点，此时的残量网络是对的
s, t = t, s # 换回来
```




