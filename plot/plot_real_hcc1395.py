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


def plot_hcc1395_recall():

    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    size = 0.2

    ax.pie([87, 13], counterclock=True, startangle=90, radius=1, colors=['#d2a993', '#ffffff'], wedgeprops=dict(width=size, edgecolor='w'))
    # ax.pie([90, 10], counterclock=True, startangle=90, radius=1, colors=['#d2a993', '#ffffff'], wedgeprops=dict(width=size, edgecolor='w'))
    # ax.pie([89, 11], counterclock=True, startangle=90, radius=1, colors=['#d2a993', '#ffffff'], wedgeprops=dict(width=size, edgecolor='w'))

    ax.pie([66, 34], counterclock=True, startangle=90, radius=1 - size, colors=['#a0aebf', '#ffffff'],  wedgeprops=dict(width=size, edgecolor='w'))
    # ax.pie([81, 19], counterclock=True, startangle=90, radius=1 - size, colors=['#a0aebf', '#ffffff'],  wedgeprops=dict(width=size, edgecolor='w'))
    # ax.pie([68, 32], counterclock=True, startangle=90, radius=1 - size, colors=['#a0aebf', '#ffffff'],  wedgeprops=dict(width=size, edgecolor='w'))


    ax.pie([6, 94], counterclock=True, startangle=90, radius=1 - size * 2, colors=['#abc7de', '#ffffff'],  wedgeprops=dict(width=size, edgecolor='w'))
    # ax.pie([0, 100], counterclock=True, startangle=90, radius=1 - size * 2, colors=['#abc7de', '#ffffff'],  wedgeprops=dict(width=size, edgecolor='w'))
    # ax.pie([29, 71], counterclock=True, startangle=90, radius=1 - size * 2, colors=['#abc7de', '#ffffff'],  wedgeprops=dict(width=size, edgecolor='w'))

    # ax.pie([30, 70], radius=1 - size,
    #        colors=['#cccccc', '#fdae6b'],
    #        labels=[0, 30], labeldistance=0.7, wedgeprops=dict(width=size, edgecolor='w'))
    # plt.show()

    # plt.savefig
    plt.savefig("out_hcc1395_somatic_accuracy.png", dpi=1000, bbox_inches='tight',)
    plt.savefig("out_hcc1395_somatic_accuracy.svg", dpi=1000, bbox_inches='tight',)


if __name__ == '__main__':


    plot_hcc1395_recall()

    # plot_vapor_validation()