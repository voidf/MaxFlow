import datetime
import os
from loguru import logger
logger.add(open('log4.txt', 'a'))

from cyaron import * # 引入CYaRon的库

# exe1 = 'copied'
# exe2 = 'lgmy'
# exe1 = 'gomory_hu_tree.py'
# exe1 = 'zxg.exe'
# exe2 = 'ao.exe'
# exe3 = 'diniczxgs.exe'
# exe1 = './diniczxgs.elf'
# exe2 = 'zxg.exe'

exes = [
    # 'zxg.exe',
    # 'ao.exe',
    # 'diniczxgs.exe',
    # 'aomodule.exe',
    'aomodulell.exe',
    'gold.exe',
    # 'aomodulev1.exe',
    # 'LMESmodule.exe',
    # 'LMESmodulev2.exe',
    # 'aomodulev2.exe',
    # 'aomodule_slower_arc.exe',
    # 'aomodule_slower_arcll.exe',
]

# 这是一个图论题的数据生成器，该题目在洛谷的题号为P1339
for T in range(400, 9999): # 即在[1, 4)范围内循环，也就是从1到3
    # fn = f'T{T:03d}'
    fn = f'T'
    test_data = IO(file_prefix=fn) # 生成 heat[1|2|3].in/out 三组测试数据
    tim = datetime.datetime.now()

    n = randint(100, 100)
    m = randint(n - 1, n * (n - 1)//2) # 边数
    # m = randint(n - 1, 2 * n) # 稀疏图
    # U = 1000_000_000
    U = 1000
    # s = randint(1, n) # 源点，随机选取一个
    # s = randint(1, n//2) # 源点，随机选取一个
    # t = randint(n//2+1, n) # 汇点，随机选取一个
    # t = randint(1, n) # 汇点，随机选取一个
    test_data.input_writeln(n, m) 
    # test_data.input_writeln(n, m, s, t) # 写入到输入文件里，自动以空格分割并换行

    graph = Graph.graph(n, m, weight_limit=(1, U),
        directed=False,
        self_loop=False, 
        repeated_edges=False) # 生成一个n点，m边的随机图，边权限制为5
    U = 1
    for ed in graph.iterate_edges():
        U = max(U, ed.weight)
    # while len(graph.edges[s]) == 0 or len(graph.edges[t]) == 0:
    #     graph = Graph.graph(n, m, weight_limit=(5, 300),
    #         directed=True,
    #         self_loop=False, 
    #         repeated_edges=False) # 生成一个n点，m边的随机图，边权限制为5
    
    test_data.input_writeln(graph) # 自动写入到输入文件里，默认以一行一组u v w的形式输出
    test_data.input_writeln(n*(n-1)>>1)
    for i in range(1, n+1):
        for j in range(i+1, n+1):
            test_data.input_writeln(i, j)


    test_data.flush_buffer()
    tim = datetime.datetime.now() - tim
    print(f'gen done. n={n}, m={m}, time:{tim}')

    answers = []
    log_result = []


    for exe in exes:
        t1 = datetime.datetime.now()
        res1 = os.popen(f'{exe} < {fn}.in').read().strip()
        answers.append(res1)
        t1 = datetime.datetime.now() - t1
        print(t1)
        log_result.append(f'{exe}:{t1.total_seconds()}')
    
    is_break = False
    for res in answers:
        if res!=answers[0]:
            print(f"{exe} ====>")
            print(res)
            print('DIFFERS =====')
            print(exes[0])
            print(f"{answers[0]} <====")
            is_break = True
    if is_break:
        break
    # print(i, n, m, s, t)
    logger.info(f"{T} n:{n}, m:{m}, U:{U} {' '.join(log_result)}")
    # if(i % 10 == 0): print(i)