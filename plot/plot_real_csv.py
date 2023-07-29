from matplotlib import pyplot as plt

import os
import numpy as np
import scipy as sc
import pylab as pl


def fmeasure_curve(f, p):
    return f * p / (2 * p - f)


def autolabel_line(x, y, ax):
    for a, b in zip(x, y):
        ax.annotate('n=%s' % b, xy=(b, a), xytext=(4, 0), textcoords="offset points")


def autolabel_rects(rects, ax):
    for rect in rects:
        width = rect.get_width()

        if width == 0:
            continue

        ax.annotate('%d' % width, xy=(width, rect.get_y()), xytext=(-6, 0), textcoords="offset points", ha='center')


def plot_csv_quartet_types_all(val_dict):
    type_legend = list(val_dict.keys())

    # # STEP: preprocess
    total_num_list = []
    match_num_list = []
    comp_missed_num_list = []
    comp_disorder_num_list = []
    mismatch_num_list = []

    for type in type_legend:

        total_num_list.append(val_dict[type][0])
        match_num_list.append(val_dict[type][1])
        comp_missed_num_list.append(val_dict[type][2])
        comp_disorder_num_list.append(val_dict[type][3])
        mismatch_num_list.append(val_dict[type][4])

    # # STEP: merge into one list
    all_num_list = [np.sum(match_num_list), np.sum(comp_missed_num_list), np.sum(comp_disorder_num_list), np.sum(mismatch_num_list)]
    all_name_list = ["Match", "Component missed", "Component disordered", "MisMatch"]
    all_color_lost = ["#be6839", "#579979", "#d6a250", "#bcbbdb"]

    # # STEP: plot a pie
    # First Ring (outside)
    fig, ax = plt.subplots(figsize=(6, 6))
    # fig = plt.figure(figsize=(6, 5))  # llc_x, llc_y, width, height

    ax.axis('equal')
    mypie, _ = ax.pie(all_num_list, radius=0.6, labels=["n={}".format(all_num_list[i]) for i in range(len(all_name_list))], colors=all_color_lost)
    plt.setp(mypie, width=0.3, edgecolor='white')

    plt.title("All CSV types")
    plt.savefig("out_csv_quartet_types_all.png", dpi=200)



def plot_csv_quartet_types(val_dict):
    type_legend = list(val_dict.keys())

    # # STEP: preprocess
    total_num_list = []
    match_num_list = []
    comp_missed_num_list = []
    comp_disorder_num_list = []
    mismatch_num_list = []

    for type in type_legend:

        total_num_list.append(val_dict[type][0])
        match_num_list.append(val_dict[type][1])
        comp_missed_num_list.append(val_dict[type][2])
        comp_disorder_num_list.append(val_dict[type][3])
        mismatch_num_list.append(val_dict[type][4])

    # # STEP: plot
    width = 0.16
    r1 = np.arange(len(type_legend))
    r2 = [x + width + 0.01 for x in r1]
    r3 = [x + width + 0.01 for x in r2]
    r4 = [x + width + 0.01 for x in r3]

    fig = pl.figure(figsize=(10, 6))  # llc_x, llc_y, width, height
    # fig = pl.figure()  # llc_x, llc_y, width, height

    ax = fig.add_subplot(1, 1, 1)

    ax.plot(total_num_list, [r + width for r in range(len(type_legend))], color="#7394bc", label="Total", marker='v')
    autolabel_line([r + width for r in range(len(type_legend))], total_num_list, ax)

    rects1 = ax.barh(r1, match_num_list, width, label="Match", color="#be6839")
    rects2 = ax.barh(r2, comp_missed_num_list, width, label="Component missed", color="#579979")
    rects3 = ax.barh(r3, comp_disorder_num_list, width, label="Component disordered", color="#d6a250")
    rects4 = ax.barh(r4, mismatch_num_list, width, label="MisMatch", color="#bcbbdb")

    autolabel_rects(rects1, ax)
    autolabel_rects(rects2, ax)
    autolabel_rects(rects3, ax)
    autolabel_rects(rects4, ax)

    ax.grid(b=True, color='grey', linestyle='-.', linewidth=0.5, alpha=0.2)

    ax.legend(loc='best')

    plt.yticks([r + width for r in range(len(type_legend))], type_legend, rotation=0)

    plt.title("Number of CSV", fontsize=10, fontweight='bold')
    plt.ylabel("CSV types", fontsize=10, fontweight='bold')

    plt.gcf().subplots_adjust(left=0.3)

    plt.savefig("out_csv_quartet_types.png", dpi=200)


if __name__ == '__main__':
    csv_quartet_types = {
        "Replacement": [16, 0, 14, 0, 2],
        "Inversion with deletion(s)": [10, 2, 7, 1, 0],
        "Inversion with insertion(s)": [6, 2, 3, 1, 0],
        "Replacement with inversion(s)": [6, 0, 6, 0, 0],
        "Complex duplication(s)": [5, 2, 1, 0, 2],
    }

    # plot_csv_quartet_types(csv_quartet_types)

    plot_csv_quartet_types_all(csv_quartet_types)