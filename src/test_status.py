from myopic_mces import MCES_ILP, construct_graph

if __name__ == '__main__':
    longs1 = 'CCCCCCCCCCCCCCCCCCCC(=O)OCC(COC(=O)CCCCCCCCCCCC=CCCCCCCCC)OC(=O)CCCCCCCCCCCCCC'
    longs2 = 'CCCCCCCCCCCCCCCCCCCCCC(=O)OC(COC(=O)CCCCCCCCCCCCCCCCC)COC(=O)CCCCCCCC=CCCCCCCCC'
    g1, g2 = map(construct_graph, [longs1, longs2])

    for tl in 10, 15, 30, 300:
        out = MCES_ILP(g1, g2, 10., solver='CPLEX_PY', solver_options=dict(timeLimit=tl, threads=1))
        print(f'{tl=}, {out=}')
