import pandas as pd
from matplotlib import pyplot as plt

import os
import numpy as np
import scipy as sc
import pylab as pl
import seaborn as sns
from matplotlib.lines import Line2D



def plot_quartet_identical(data_file):
    data_frame = pd.read_table(data_file)

    fig = pl.figure(figsize=(5, 2))
    ax = fig.add_subplot(1, 1, 1)
    # ax.grid(b=True, color='grey', linestyle='-.', linewidth=0.5, alpha=0.2)
    ax.grid( color='grey', linestyle='-.', linewidth=0.5, alpha=0.2)

    # sns.boxplot(data=data_frame, x="tool", y="ratio", hue="tech",   whis=0, capprops={'color': 'none'}, medianprops={'color': 'none'}, boxprops={'facecolor': '#FFFFFF', 'edgecolor': 'none'}, whiskerprops={'color': 'none'}, showfliers=False)
    # sns.boxplot(data=data_frame, x="tool", y="ratio", hue="tech",whiskerprops={'color': 'none'}  )

    sns.swarmplot(data=data_frame, x="tool", y="ratio", hue="tech", hue_order=["HiFi"], dodge=True, color="black", size=10, palette=["#d2a993", "#b6c1d0"])

    for s in ['left', 'right']:
        ax.spines[s].set_visible(False)

    # ax.set_ylim([0.0, 0.15])  # xmin, xmax, ymin, ymax
    # ax.set_yticks(np.linspace(0.0, 0.15, 4))
    ax.set_ylim([0.0, 0.12])  # xmin, xmax, ymin, ymax
    ax.set_yticks(np.linspace(0.0, 0.12, 3))
    ax.set_xticklabels(tool_legend, fontsize=16, rotation=340)

    plt.xlabel("")
    plt.ylabel("Twin discordancy", fontsize=16,)

    plt.legend([],[], frameon=False)

    plt.savefig("out_trio_twin.png", dpi=1500,bbox_inches='tight')

    plt.savefig("out_trio_twin.svg", dpi=1500,bbox_inches='tight')

def plot_trio_mendelian(data_file):

    data_frame = pd.read_table(data_file)
    #
    # data_frame = data_frame[data_frame["tech"] == "ONT"]

    fig = pl.figure(figsize=(5, 3))
    ax = fig.add_subplot(1, 1, 1)

    # ax.grid(b=True, color='grey', linestyle='-.', linewidth=0.5, alpha=0.2)
    ax.grid(color='grey', linestyle='-.', linewidth=0.5, alpha=0.2)

    sns.boxplot(data=data_frame, x="tool", y="ratio", hue="tech", palette=[ "#d2a993", "#b6c1d0"])

    for s in ['left', 'right']:
        ax.spines[s].set_visible(False)

    ax.set_ylim([0.8, 1.0])  # xmin, xmax, ymin, ymax
    ax.set_yticks(np.linspace(0.8, 1.0, 3))
    # ax.set_xticklabels([], fontsize=6,)
    ax.set_xticklabels(tool_legend, fontsize=16, rotation=340)

    # ax.set_ylim([0.85, 1.0])  # xmin, xmax, ymin, ymax
    # ax.set_yticks(np.linspace(0.85, 1.0, 3))
    # ax.set_xticklabels(tool_legend, fontsize=16, rotation=340)


    plt.xlabel("")
    plt.ylabel("Mendelian Consistency", fontsize=16,)


    plt.legend([],[], frameon=False)
    # plt.legend( frameon=False, loc=1, ncol=2,  fontsize=10,)

    plt.savefig("out_trio_mendelian.png", dpi=1500, bbox_inches='tight')
    plt.savefig("out_trio_mendelian.svg", dpi=1500, bbox_inches='tight')
    # plt.savefig("out_trio_mendelian.pdf", dpi=1500, bbox_inches='tight')
    # plt.savefig("out_trio_mendelian.eps", dpi=1500, bbox_inches='tight')

if __name__ == '__main__':
    tool_legend = ["SVision-pro", "sniffles2", "SVision",  "cuteSV", "debreak"]


    # plot_trio_mendelian("data_trio_mendelian_hifi_s5_10_ont_s10.6trio.txt")

    plot_quartet_identical("data_trio_twin_discordant_hifi_s10_ont_s10.quartet.txt")
