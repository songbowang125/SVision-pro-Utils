import pandas as pd
from matplotlib import pyplot as plt

import os
import numpy as np
import scipy as sc
import pylab as pl
import seaborn as sns
from matplotlib.lines import Line2D



# def plot_quartet_identical(quartet_hifi):
#
#     survivor_res = [quartet_hifi[2], quartet_hifi[4], quartet_hifi[6]]
#     jasmine_res = [quartet_hifi[3], quartet_hifi[5], quartet_hifi[6]]
#
#     fig = pl.figure(figsize=(5, 2))
#     ax = fig.add_subplot(1, 1, 1)
#
#     # plt.plot(quartet_hifi, marker=".", color="#d2a993")
#
#     plt.scatter([0, 1], quartet_hifi[0: 2], marker="v", ms=10, color="#d8d8d8")
#
#     plt.scatter([2, 3, 4], survivor_res, marker="v", ms=10, color="#b6c1d0")
#     plt.scatter([2, 3, 4], jasmine_res, marker="v", ms=10, color="#d2a993")
#
#     for s in ['left', 'right']:
#         ax.spines[s].set_visible(False)
#
#     ax.set_ylim([0., 0.1])  # xmin, xmax, ymin, ymax
#
#     ax.set_yticks(np.linspace(0., 0.1, 3))
#     ax.grid(b=True, color='grey', linestyle='-.', linewidth=0.5, alpha=0.2)
#
#     # ax.set_yticklabels([round(val, 2) for val in np.linspace(0.6, 1.0, 5)], fontsize=8)
#
#     # plt.ylabel("Twin discordance")
#     plt.xticks(range(0, len(tool_legend)), tool_legend, fontsize=11)
#
#     pl.savefig("out_trio_twin.png", dpi=400, bbox_inches='tight',)


def plot_quartet_identical(data_file):
    data_frame = pd.read_table(data_file)

    fig = pl.figure(figsize=(5, 2))
    ax = fig.add_subplot(1, 1, 1)
    ax.grid(b=True, color='grey', linestyle='-.', linewidth=0.5, alpha=0.2)

    # sns.boxplot(data=data_frame, x="tool", y="ratio", hue="tech",   whis=0, capprops={'color': 'none'}, medianprops={'color': 'none'}, boxprops={'facecolor': '#FFFFFF', 'edgecolor': 'none'}, whiskerprops={'color': 'none'}, showfliers=False)
    # sns.boxplot(data=data_frame, x="tool", y="ratio", hue="tech",whiskerprops={'color': 'none'}  )

    sns.swarmplot(data=data_frame, x="tool", y="ratio", hue="tech", hue_order=["HiFi"], dodge=True, color="black", size=10, palette=["#d2a993", "#b6c1d0"])

    for s in ['left', 'right']:
        ax.spines[s].set_visible(False)

    ax.set_ylim([0.0, 0.12])  # xmin, xmax, ymin, ymax
    ax.set_yticks(np.linspace(0.0, 0.12, 3))
    ax.set_xticklabels(tool_legend, fontsize=16, rotation=340)

    plt.xlabel("")
    plt.ylabel("Twin discordancy", fontsize=16,)

    plt.legend([],[], frameon=False)

    plt.savefig("out_trio_twin.png", dpi=1500,bbox_inches='tight')


def plot_trio_mendelian(data_file):

    data_frame = pd.read_table(data_file)

    fig = pl.figure(figsize=(5, 3))
    ax = fig.add_subplot(1, 1, 1)
    ax.grid(b=True, color='grey', linestyle='-.', linewidth=0.5, alpha=0.2)

    sns.boxplot(data=data_frame, x="tool", y="ratio", hue="tech", palette=[ "#d2a993", "#b6c1d0"])

    for s in ['left', 'right']:
        ax.spines[s].set_visible(False)

    ax.set_ylim([0.8, 1.0])  # xmin, xmax, ymin, ymax
    ax.set_yticks(np.linspace(0.8, 1.0, 3))
    ax.set_xticklabels([], fontsize=6,)

    plt.xlabel("")
    plt.ylabel("Mendelian Consistency", fontsize=16,)


    plt.legend([],[], frameon=False)
    # plt.legend( frameon=False, loc=1, ncol=2,  fontsize=10,)

    plt.savefig("out_trio_mendelian.png", dpi=1500, bbox_inches='tight')
    plt.savefig("out_trio_mendelian.svg", dpi=1500, bbox_inches='tight')

if __name__ == '__main__':
    tool_legend = ["SVision-pro", "sniffles2", "SVision",  "cuteSV", "debreak"]

    # tool_legend = ["SVision-pro", "sniffles2", "SVision-SURVIVOR", "SVision-Jasmine", "cuteSV-SURVIVOR", "cuteSV-Jasmine", "debreak-SURVIVOR", "debreak-Jasmine"]
    # hg002_hifi = [0.957053897, 0.938466565, 0.942819441, 0.940449558, 0.930104568, 0.910244007, 0.874078151, 0.87313037]
    # hg002_ont = [0.974420902, 0.974365117, 0.946248476, 0.952910681, 0.962843296, 0.938426712, 0.92653329, 0.925302754]
    #
    # quartet_hifi = [0.009287421, 0.042985047, 0.066689072, 0.084997259, 0.059050325, 0.076366436, 0.08692053, 0.092397887]
    # quartet_ont = [0.009287421, 0.042985047, 0.066689072, 0.084997259, 0.059050325, 0.076366436, 0.08692053, 0.092397887]
    #
    #
    #
    # # plot_hg002_mendelian(hg002_hifi, hg002_ont)
    # #
    # #

    plot_trio_mendelian("data_trio_mendelian_hifi_s5_10_ont_s10.6trio.txt")
    #
    plot_quartet_identical("data_trio_twin_discordant_hifi_s10_ont_s10.quartet.txt")
