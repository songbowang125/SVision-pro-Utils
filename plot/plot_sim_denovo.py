import os
import pandas as pd
from matplotlib import pyplot as plt
import math
import os
import numpy as np
import scipy as sc
import pylab as pl
import seaborn as sns
from scipy.stats import ranksums
import random




# def plot_sim_denovo_accuracy(sim_hifi, sim_ont):
#
#     fig = pl.figure(figsize=(1.8, 5))
#     ax = fig.add_subplot(1, 1, 1)
#
#     survivor_res_hifi = [sim_hifi[1]]
#     jasmine_res_hifi = [sim_hifi[2]]
#
#     survivor_res_ont = [sim_ont[1]]
#     jasmine_res_ont = [sim_ont[2]]
#
#
#     # plt.plot(hg002_hifi, marker=".", color="#d2a993")
#
#     plt.scatter([0], sim_hifi[0], marker=".", color="#d2a993", s=200)
#     plt.scatter([1], survivor_res_hifi, marker="s", color="#d2a993", s=200)
#     plt.scatter([1], jasmine_res_hifi, marker="p", color="#d2a993", s=200)
#
#     # plt.plot(hg002_ont, marker=".", color="#b6c1d0")
#     plt.scatter([0], sim_ont[0], marker=".", color="#b6c1d0", s=200)
#     plt.scatter([1], jasmine_res_ont, marker="p", color="#b6c1d0", s=120)
#     plt.scatter([1], survivor_res_ont, marker="s", color="#b6c1d0", s=120)
#
#     for s in ['left', 'right']:
#         ax.spines[s].set_visible(False)
#
#     # ax.set_ylim([0.85, 1.0])  # xmin, xmax, ymin, ymax
#
#     ax.set_yticks(np.linspace(0.4, 1.0, 4))
#     ax.grid(b=True, color='grey', linestyle='-.', linewidth=0.5, alpha=0.2)
#
#     # ax.set_yticklabels([round(val, 2) for val in np.linspace(0.6, 1.0, 5)], fontsize=8)
#
#     plt.ylabel("CSV Mendelian Consistency", fontsize=16)
#     plt.xticks(range(0, len(tool_legend)), tool_legend, fontsize=16, rotation=340)
#
#
#     # ax.legend(handles=[Line2D([0], [0], label='HiFi', color='#d2a993', ls='-', lw=2), Line2D([0], [0], label='ONT', color='#b6c1d0', ls='-', lw=2)])
#     # plt.legend([Line2D([0], [0], label='HiFi', color='#d2a993', ls='-', lw=2)], [Line2D([0], [0], label='ONT', color='#b6c1d0', ls='-', lw=2)])
#     pl.savefig("out_sim_denovo.png", dpi=200, bbox_inches='tight',)



def plot_sim_denovo_accuracy(data_file):
    data_frame = pd.read_table(data_file)

    fig = pl.figure(figsize=(2, 5))
    ax = fig.add_subplot(1, 1, 1)
    ax.grid(b=True, color='grey', linestyle='-.', linewidth=0.5, alpha=0.2)

    # sns.boxplot(data=data_frame, x="tool", y="ratio", hue="tech",   whis=0, capprops={'color': 'none'}, medianprops={'color': 'none'}, boxprops={'facecolor': '#FFFFFF', 'edgecolor': 'none'}, whiskerprops={'color': 'none'}, showfliers=False)
    # sns.boxplot(data=data_frame, x="tool", y="ratio", hue="tech",whiskerprops={'color': 'none'}  )

    sns.swarmplot(data=data_frame, x="tool", y="ratio", hue="tech", dodge=True, color="black", size=10, palette=["#d2a993", "#b6c1d0"])

    for s in ['left', 'right']:
        ax.spines[s].set_visible(False)

    # ax.set_ylim([0.0, 0.12])  # xmin, xmax, ymin, ymax
    # ax.set_yticks(np.linspace(0.0, 0.12, 3))
    ax.set_xticklabels(["SVision-pro", "SVision"], fontsize=16, rotation=340)

    plt.xlabel("")
    plt.ylabel("CSV Mendelian Consistency", fontsize=16,)

    plt.legend([],[], frameon=False)

    plt.savefig("out_sim_denovo.png", dpi=1500,bbox_inches='tight')


if __name__ == '__main__':

    tool_legend = ["SVision-pro", "SVision",]

    # tool_legend = ["SVision-pro", "SVision-SURVIVOR", "SVision-Jasmine", ]
    # sim_hifi = [0.966690804, 0.532340116, 0.540864511]
    # sim_ont = [0.933986217, 0.335692996, 0.332960894]

    # plot_sim_denovo_accuracy(sim_hifi, sim_ont)

    plot_sim_denovo_accuracy("data_sim_trio.txt")