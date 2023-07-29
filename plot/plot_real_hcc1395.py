from matplotlib import pyplot as plt



def plot_hcc1395_recall():

    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    size = 0.2

    ax.pie([87, 13], counterclock=True, startangle=90, radius=1, colors=['#d2a993', '#ffffff'], wedgeprops=dict(width=size, edgecolor='w'))
    # [np-sigs, csv, inv, fewer-sig, others, novel]
    ax.pie([72, 28], counterclock=True, startangle=90, radius=1 - size, colors=['#95a6cf', '#ffffff'],  wedgeprops=dict(width=size, edgecolor='w'))

    ax.pie([6, 94], counterclock=True, startangle=90, radius=1 - size * 2, colors=['#b6c1d0', '#ffffff'],  wedgeprops=dict(width=size, edgecolor='w'))

    # ax.pie([30, 70], radius=1 - size,
    #        colors=['#cccccc', '#fdae6b'],
    #        labels=[0, 30], labeldistance=0.7, wedgeprops=dict(width=size, edgecolor='w'))
    plt.show()

    # plt.savefig
    plt.savefig("out_hcc1395_somatic_accuracy.png", dpi=1500, bbox_inches='tight',)


if __name__ == '__main__':
    plot_hcc1395_recall()