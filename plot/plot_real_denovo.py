import os
import pandas as pd
from matplotlib import pyplot as plt

import os
import numpy as np


def plot_denovo_number():

    val_dict = {
        "debreak_SURVIVOR": 8269,
        "debreak_Jasmine": 10597,
        "cuteSV_SURVIVOR": 5831,
        "cuteSV_Jasmine": 12074,
        "SVision_SURVIVOR": 6813,
        "SVision_Jasmine": 12468,
        "sniffles2": 90,
        "SVision-pro": 26,

    }


    fig, ax = plt.subplots(1, 1, figsize=(2, 5))
    # fig, ax = plt.subplots(1, 1, figsize=(4, 5))

    width = 0.5
    ax.barh(range(len(val_dict.keys())), val_dict.values(), width, color="#caab9b")

    plt.yticks(range(len(val_dict.keys())), val_dict.keys(), fontsize=16, rotation=0)

    # ax.grid(b=True, color='grey', linestyle='-.', linewidth=0.5, alpha=0.2)
    ax.grid( color='#bfbfbf', linestyle='-.', linewidth=0.2, alpha=0.2)


    for s in ['top', 'bottom', 'right']:
        ax.spines[s].set_visible(False)
    ax.set_xlim(0, 100)
    ax.set_xticks(np.linspace(0, 100, 3))



    # for s in ['top', 'bottom', 'right', "left"]:
    #     ax.spines[s].set_visible(False)
    # ax.set_xlim(4000, 13000)
    # ax.set_xticks(np.linspace(5000, 13000, 5))


    #
    plt.xlabel("")
    # plt.ylabel("Twin discordancy", fontsize=16,)

    plt.savefig("out_trio_denovo_num1.png", dpi=1500, bbox_inches='tight')
    plt.savefig("out_trio_denovo_num1.svg", dpi=1500, bbox_inches='tight')


def plt_denovo_fp_sniffles():

    fig, ax = plt.subplots(1, 1, figsize=(3, 3))
    size = 0.3

    wedges, _ = ax.pie([1, 59, 11, 22], counterclock=True, startangle=90, colors=['#a0aebf', '#b4aeb4', '#ebbca5', '#b4aeb4'], wedgeprops=dict(width=size, edgecolor='w'))

    # plt.text(0.5, 0.5, "Denovo SVs \nof sniffles (n=95)", horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=8)

    wedges[3].set_hatch('+')
    plt.savefig("out_trio_denovo_sniffles.png", dpi=1500, bbox_inches='tight')
    plt.savefig("out_trio_denovo_sniffles.svg", dpi=1500, bbox_inches='tight')


if __name__ == '__main__':
    # plot_denovo_fp()
    # plt_denovo_fp_sniffles()

    plot_denovo_number()