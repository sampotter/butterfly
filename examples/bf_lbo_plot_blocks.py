import matplotlib.pyplot as plt; plt.ion()

from matplotlib.patches import Rectangle

from glob import glob

paths = ['Psi.txt'] + glob('W*.txt')

type2str = {
    6: 'MatIdentity',
    11: 'MatBlockCoo',
    13: 'MatBlockDiag',
    16: 'MatDenseReal'
}

type2color = {
    6: 'white',
    11: 'cyan',
    13: 'magenta',
    16: 'red'
}

def plot_blocks(path):
    f = open(path)
    fig = plt.gcf()
    ax = plt.gca()
    I, J = set(), set()
    for line in f:
        type_, i0, i1, j0, j1, depth = map(int, line.split(' '))
        I.add(i0)
        I.add(i1)
        J.add(j0)
        J.add(j1)
        m, n = i1 - i0, j1 - j0
        facecolor = type2color[type_]
        rect = Rectangle((j0, i0), n, m, facecolor=facecolor, edgecolor='k',
                         linewidth=1, zorder=depth)
        ax.add_patch(rect)
    I, J = sorted(list(I)), sorted(list(J))
    M, N = I[-1], J[-1]
    ax.set_xlim(0, N)
    ax.set_ylim(0, M)
    ax.set_xticks(J)
    ax.set_yticks(I)
    fac_name = path.split('.')[0]
    plt.get_current_fig_manager().set_window_title(fac_name)
    ax.invert_yaxis()
    ax.set_aspect('equal')
    f.close()

for path in paths:
    plt.figure()
    plot_blocks(path)
    plt.show()
