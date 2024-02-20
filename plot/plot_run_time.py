from matplotlib import pyplot as plt
from matplotlib import pyplot as plt
import math
import os
import numpy as np
import scipy as sc
import pylab as pl
import seaborn as sns
from scipy.stats import ranksums
import random


def plot_memory():

    x_axis = ["5x", "10x", "20x", "40x"]

    memory_hifi = [0.89, 1.28, 1.78, 2.29]
    memory_ont = [0.98, 1.84, 2.3, 3.68]

    # # STEP: plot
    r1 = np.arange(len(x_axis))

    # fig = pl.figure(figsize=(2.5, 3))  # llc_x, llc_y, width, height
    fig = pl.figure(figsize=(3, 1.2))  # llc_x, llc_y, width, height

    ax = fig.add_subplot(1, 1, 1)

    # # for HiFi
    marker = "."
    line = "-"
    ax.plot(r1, memory_hifi, label="PC_thread_8", color="#bf9000", marker=marker, ls=line)

    marker = "."
    line = "--"
    ax.plot(r1, memory_ont, label="PC_thread_8", color="#bf9000", marker=marker, ls=line)

    # ax.grid(b=True, color='grey', linestyle='-.', linewidth=0.5, alpha=0.2)
    ax.grid(color='grey', linestyle='-.', linewidth=0.5, alpha=0.2)

    # ax.legend(loc='best', fontsize=7,)

    # Remove axes splines
    for s in ['left', 'right']:
        ax.spines[s].set_visible(False)

    plt.xticks([r for r in range(len(x_axis))], x_axis, fontsize=8, rotation=340)

    # # for HiFi
    plt.savefig("out_memory.png", dpi=1500, bbox_inches='tight',)

    # # for ONT
    # plt.savefig("out_time.png", dpi=1500, bbox_inches='tight',)

def plot_time():

    x_axis = ["5x", "10x", "20x", "40x"]

    time_hifi = {"PC_thread_8": [9, 21, 33, 47],
               "Cluster_thread_24": [3, 8, 13, 18],}

    time_ont = {"PC_thread_8": [24, 38, 65, 135],
               "Cluster_thread_24": [15, 26, 36, 77],
}



    # time_ont = {"PC_thread_8": [24, 38, 65, 135],
    #            "Cluster_thread_24": [25, 26, 36, 77],
    #            "Cluster_thread_48": [23, 21, 26, 57]}

    # # STEP: plot
    r1 = np.arange(len(x_axis))

    # fig = pl.figure(figsize=(2.5, 3))  # llc_x, llc_y, width, height
    fig = pl.figure(figsize=(3, 1.2))  # llc_x, llc_y, width, height

    ax = fig.add_subplot(1, 1, 1)

    # # for HiFi
    marker = "."
    line = "-"
    ax.plot(r1, np.array(time_hifi["PC_thread_8"]), label="PC_thread_8", color="#78a2e7", marker=marker, ls=line)
    ax.plot(r1, np.array(time_hifi["Cluster_thread_24"]), label="Cluster_thread_24", color="#9e89d6", marker=marker, ls=line)
    # ax.plot(r1, np.array(time_hifi["Cluster_thread_48"]), label="Cluster_thread_48", color="#e98d76", marker=marker, ls=line)

    # # # for ONT
    # marker = "."
    # line = "--"
    # ax.plot(r1, np.array(time_ont["PC_thread_8"]), label="PC_thread_8", color="#78a2e7", marker=marker, ls=line)
    # ax.plot(r1, np.array(time_ont["Cluster_thread_24"]), label="Cluster_thread_24", color="#9e89d6", marker=marker, ls=line)
    # ax.plot(r1, np.array(time_ont["Cluster_thread_48"]), label="Cluster_thread_48", color="#e98d76", marker=marker, ls=line)

    # ax.grid(b=True, color='grey', linestyle='-.', linewidth=0.5, alpha=0.2)
    ax.grid(color='grey', linestyle='-.', linewidth=0.5, alpha=0.2)

    # ax.legend(loc='best', fontsize=7,)

    # Remove axes splines
    for s in ['left', 'right']:
        ax.spines[s].set_visible(False)

    plt.xticks([r for r in range(len(x_axis))], x_axis, fontsize=8, rotation=340)

    # # for HiFi
    ax.set_xticklabels([], fontsize=10)
    plt.savefig("out_time.png", dpi=1500, bbox_inches='tight',)

    # for ONT
    # plt.savefig("out_time.png", dpi=1500, bbox_inches='tight',)


if __name__ == '__main__':
    plot_time()
    # plot_memory()