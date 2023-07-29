from matplotlib import pyplot as plt
import pandas as pd
import os
import numpy as np
import scipy as sc
import pylab as pl
import seaborn as sns


def fmeasure_curve(f, p):
    return f * p / (2 * p - f)


def plot_fmeasures(fstepsize=.05, stepsize=0.001):
    """Plots 10 fmeasure Curves into the current canvas."""
    p = sc.arange(0., 1., stepsize)[1:]
    for f in sc.arange(0., 1., fstepsize)[1:]:
        points = [(x, fmeasure_curve(f, x)) for x in p if 0 < fmeasure_curve(f, x) <= 1.5]
        xs, ys = zip(*points)
        curve, = pl.plot(xs, ys, "--", color="gray", linewidth=0.5)  # , label=r"$f=%.1f$"%f) # exclude labels, for legend
        # bad hack:
        # gets the 10th last datapoint, from that goes a bit to the left, and a bit down
        pl.annotate(r"$f1=%.2f$" % f, xy=(xs[-10], ys[-10]), xytext=(xs[-10] - 0.03, ys[-10] + 0.01), size="small", color="gray")


def plot_hg002_high_confident_precision_recall(hifi_val_dict, ont_val_dict, mode):

    seq_legend = ["HiFi", "ONT"]

    tool_legend = list(hifi_val_dict.keys())

    # # STEP: set up fig config
    fig = pl.figure(figsize=(6, 5))  # llc_x, llc_y, width, height
    ax = fig.add_subplot(1, 1, 1)

    # # STEP: plot f1 line
    # ax.xaxis.set_ticks_position('top')
    pl.title("HG002 high confident SSVs".format(mode), fontsize=15, fontweight='bold')
    pl.xlabel("Precision", fontsize=13, fontweight='bold')
    pl.ylabel("Recall", fontsize=13, fontweight='bold')
    plot_fmeasures()

    # # STEP: plot legends
    tool_legend_plot = [pl.scatter([], [], marker='s', s=100, color=tool_color[tool_legend[i]], label=tool_legend[i]) for i in range(len(tool_legend))]

    # seq_legend_plot = [pl.scatter([], [], marker=seq_symbol[seq_legend[i]], s=80, color='gray', label=seq_legend[i]) for i in range(len(seq_legend))]
    seq_legend_plot = [pl.scatter([], [], marker=seq_symbol[seq_legend[i]], s=80, facecolor='none', edgecolors='black', linewidths=1, label=seq_legend[i]) for i in range(len(seq_legend))]

    for patch in seq_legend_plot:
        tool_legend_plot.append(patch)
    pl.legend(handles=tool_legend_plot, loc=(0.02, 0.02), fancybox=True)

    # # STEP: plot hifi data points
    for tool, points in hifi_val_dict.items():
        if "NA" in points:
            continue
        ax.scatter(points[0], points[1], label=tool, s=90, linewidths=0.75, facecolor=tool_color[tool], alpha=0.75, marker=seq_symbol["HiFi"])

    # # STEP: plot ont data points
    for tool, points in ont_val_dict.items():
        if "NA" in points:
            continue
        ax.scatter(points[0], points[1], label=tool, s=90, linewidths=0.75, facecolor=tool_color[tool], alpha=0.75, marker=seq_symbol["ONT"])

    pl.tight_layout()
    pl.axis([0.7, 1.02, 0.7, 1.02])  # xmin, xmax, ymin, ymax
    # pl.show()
    pl.savefig("out_hg002_{}_high_confident.png".format(mode), dpi=1500)



def plot_hg002_high_confident_f1_gt(hifi_val_dict, ont_val_dict, mode):

    # # STEP: plot f1
    f1_hifi_values = [hifi_val_dict[tool][2] for tool in hifi_val_dict.keys()]
    f1_ont_values = [ont_val_dict[tool][2] for tool in hifi_val_dict.keys()]

    # gt_values = [hifi_val_dict[tool][3] for tool in hifi_val_dict.keys()]

    fig = plt.figure(figsize=(5, 2))  # llc_x, llc_y, width, height
    ax = fig.add_subplot(1, 1, 1)

    tool_legend = list(hifi_val_dict.keys())

    # # STEP: plot dash lines
    f1_hifi_values_without_na = [val for val in f1_hifi_values if val != "NA"]
    f1_ont_values_without_na = [val for val in f1_ont_values if val != "NA"]

    ax.plot(f1_hifi_values_without_na, color="gray", label="-", linewidth=0.25,)
    ax.plot(f1_ont_values_without_na, color="gray", label="-", linewidth=0.25,)

    # # STEP: plot data points for hifi
    for i in range(len(tool_legend)):
        tool = tool_legend[i]
        if f1_hifi_values[i] == "NA":
            continue
        ax.scatter(i, f1_hifi_values[i], label=tool, s=50, linewidths=0.75, facecolor=tool_color[tool], alpha=0.75, marker=seq_symbol["HiFi"])

    # # STEP: plot data points for ont
    for i in range(len(tool_legend)):
        tool = tool_legend[i]
        if f1_ont_values[i] == "NA" or i >= len(tool_legend):
            continue
        ax.scatter(i, f1_ont_values[i], label=tool, s=50, linewidths=0.75, facecolor=tool_color[tool], alpha=0.75, marker=seq_symbol["ONT"])

    pl.ylabel("F1", fontsize=13, fontweight='bold')
    pl.xticks([])

    # # STEP: plot GT
    gt_hifi_values = [hifi_val_dict[tool][3] for tool in hifi_val_dict.keys()]
    gt_ont_values = [ont_val_dict[tool][3] for tool in hifi_val_dict.keys()]

    # ax2 = fig.add_subplot(2, 1, 2)
    #
    # # # STEP: plot dash lines
    # gt_hifi_values_without_na = [val for val in gt_hifi_values if val != "NA"]
    # gt_ont_values_without_na = [val for val in gt_ont_values if val != "NA"]
    #
    # ax2.plot(gt_hifi_values_without_na, color="gray", label="-", linewidth=0.25,)
    # ax2.plot(gt_ont_values_without_na, color="gray", label="-", linewidth=0.25,)
    #
    # # # STEP: plot data points for hifi
    # for i in range(len(tool_legend)):
    #     tool = tool_legend[i]
    #     if gt_hifi_values[i] == "NA":
    #         continue
    #     ax2.scatter(i, gt_hifi_values[i], label=tool, s=50, linewidths=0.75, facecolor=tool_color[tool], alpha=0.75, marker=seq_symbol["HiFi"])
    #
    # # # STEP: plot data points for ont
    # for i in range(len(tool_legend)):
    #     tool = tool_legend[i]
    #     if gt_ont_values[i] == "NA" or i >= len(tool_legend):
    #         continue
    #     ax2.scatter(i, gt_ont_values[i], label=tool, s=50, linewidths=0.75, facecolor=tool_color[tool], alpha=0.75, marker=seq_symbol["ONT"])
    # pl.ylabel("GT", fontsize=10, fontweight='bold')
    # pl.xticks(range(0, len(tool_legend), 1), tool_legend)

    pl.savefig("out_hg002_{}_high_confident_f1_gt.png".format(mode), dpi=1500)






def plot_csv_simulated_precision_recall(hifi_val_dict, ont_val_dict, mode):

    seq_legend = ["HiFi", "ONT"]

    tool_legend = list(hifi_val_dict.keys())

    # # STEP: set up fig config
    fig = pl.figure(figsize=(6, 5))  # llc_x, llc_y, width, height
    ax = fig.add_subplot(1, 1, 1)

    # # STEP: plot f1 line
    # ax.xaxis.set_ticks_position('top')
    pl.title("Simualted CSVs (Region Match)".format(mode), fontsize=15, fontweight='bold')
    pl.xlabel("Precision", fontsize=13, fontweight='bold')
    pl.ylabel("Recall", fontsize=13, fontweight='bold')
    plot_fmeasures()

    # # STEP: plot legends
    tool_legend_plot = [pl.scatter([], [], marker='s', s=100, color=tool_color[tool_legend[i]], label=tool_legend[i]) for i in range(len(tool_legend))]

    # seq_legend_plot = [pl.scatter([], [], marker=seq_symbol[seq_legend[i]], s=80, color='gray', label=seq_legend[i]) for i in range(len(seq_legend))]
    seq_legend_plot = [pl.scatter([], [], marker=seq_symbol[seq_legend[i]], s=80, facecolor='none', edgecolors='black', linewidths=1, label=seq_legend[i]) for i in range(len(seq_legend))]

    for patch in seq_legend_plot:
        tool_legend_plot.append(patch)
    pl.legend(handles=tool_legend_plot, loc=(0.6, 0.05), fancybox=True)

    # # STEP: plot hifi data points
    for tool, points in hifi_val_dict.items():
        if "NA" in points:
            continue
        ax.scatter(points[0], points[1], label=tool, s=90, linewidths=0.75, facecolor=tool_color[tool], alpha=0.75, marker=seq_symbol["HiFi"])

    # # STEP: plot ont data points
    for tool, points in ont_val_dict.items():
        if "NA" in points:
            continue
        ax.scatter(points[0], points[1], label=tool, s=90, linewidths=0.75, facecolor=tool_color[tool], alpha=0.75, marker=seq_symbol["ONT"])

    pl.tight_layout()
    pl.axis([0.3, 1.02, 0.3, 1.02])  # xmin, xmax, ymin, ymax
    # pl.show()
    pl.savefig("out_csv_{}_simulated.png".format(mode), dpi=1500)


def plot_csv_simualted_f1_gt(hifi_val_dict, ont_val_dict, mode):

    # # STEP: plot f1
    f1_hifi_values = [hifi_val_dict[tool][2] for tool in hifi_val_dict.keys()]
    f1_ont_values = [ont_val_dict[tool][2] for tool in hifi_val_dict.keys()]

    # gt_values = [hifi_val_dict[tool][3] for tool in hifi_val_dict.keys()]

    fig = plt.figure(figsize=(5, 2))  # llc_x, llc_y, width, height
    ax = fig.add_subplot(1, 1, 1)

    tool_legend = list(hifi_val_dict.keys())

    # # STEP: plot dash lines
    f1_hifi_values_without_na = [val for val in f1_hifi_values if val != "NA"]
    f1_ont_values_without_na = [val for val in f1_ont_values if val != "NA"]

    ax.plot(f1_hifi_values_without_na, color="gray", label="-", linewidth=0.25,)
    ax.plot(f1_ont_values_without_na, color="gray", label="-", linewidth=0.25,)

    # # STEP: plot data points for hifi
    for i in range(len(tool_legend)):
        tool = tool_legend[i]
        if f1_hifi_values[i] == "NA":
            continue
        ax.scatter(i, f1_hifi_values[i], label=tool, s=50, linewidths=0.75, facecolor=tool_color[tool], alpha=0.75, marker=seq_symbol["HiFi"])

    # # STEP: plot data points for ont
    for i in range(len(tool_legend)):
        tool = tool_legend[i]
        if f1_ont_values[i] == "NA" or i >= len(tool_legend):
            continue
        ax.scatter(i, f1_ont_values[i], label=tool, s=50, linewidths=0.75, facecolor=tool_color[tool], alpha=0.75, marker=seq_symbol["ONT"])

    pl.ylabel("F1", fontsize=13, fontweight='bold')
    pl.xticks([])

    # # STEP: plot GT
    gt_hifi_values = [hifi_val_dict[tool][3] for tool in hifi_val_dict.keys()]
    gt_ont_values = [ont_val_dict[tool][3] for tool in hifi_val_dict.keys()]

    # ax2 = fig.add_subplot(2, 1, 2)
    #
    # # # STEP: plot dash lines
    # gt_hifi_values_without_na = [val for val in gt_hifi_values if val != "NA"]
    # gt_ont_values_without_na = [val for val in gt_ont_values if val != "NA"]
    #
    # ax2.plot(gt_hifi_values_without_na, color="gray", label="-", linewidth=0.25,)
    # ax2.plot(gt_ont_values_without_na, color="gray", label="-", linewidth=0.25,)
    #
    # # # STEP: plot data points for hifi
    # for i in range(len(tool_legend)):
    #     tool = tool_legend[i]
    #     if gt_hifi_values[i] == "NA":
    #         continue
    #     ax2.scatter(i, gt_hifi_values[i], label=tool, s=50, linewidths=0.75, facecolor=tool_color[tool], alpha=0.75, marker=seq_symbol["HiFi"])
    #
    # # # STEP: plot data points for ont
    # for i in range(len(tool_legend)):
    #     tool = tool_legend[i]
    #     if gt_ont_values[i] == "NA" or i >= len(tool_legend):
    #         continue
    #     ax2.scatter(i, gt_ont_values[i], label=tool, s=50, linewidths=0.75, facecolor=tool_color[tool], alpha=0.75, marker=seq_symbol["ONT"])
    # pl.ylabel("GT", fontsize=10, fontweight='bold')
    # pl.xticks(range(0, len(tool_legend), 1), tool_legend)

    pl.savefig("out_csv_{}_simualted_f1_gt.png".format(mode), dpi=1500)



def plot_csv_structure_boxplot(csv_structure):
    csv_structure_pd = pd.DataFrame(csv_structure)

    print(csv_structure_pd)

    fig = pl.figure(figsize=(2, 2.5))

    ax = fig.add_subplot(1, 1, 1)

    sns.boxplot(data=csv_structure_pd, fliersize=2, color="#fc956c", width=0.7,  linewidth=1)

    # Remove axes splines
    for s in ['left', 'right']:
        ax.spines[s].set_visible(False)

    plt.ylabel("structure concordance")

    plt.xticks(fontsize=8, rotation=340)
    plt.yticks(fontsize=8,)

    ax.grid(b=True, color='grey', linestyle='-.', linewidth=0.5, alpha=0.2)

    pl.title("Simualted CSVs (Structure Match)".format(mode), fontsize=9, fontweight='bold')

    pl.savefig("out_csv_{}_structure_concordance.png".format(mode), dpi=1500, bbox_inches='tight',)



if __name__ == '__main__':
    seq_symbol = {"HiFi": "o", "ONT": "^"}  # # for hifi and ont, use different symbles
    tool_color = {'SVision-pro': '#4386bc', 'SVision': '#91a3cd', 'cuteSV': '#ff8710', 'Sniffles2': '#f889c3', 'pbsv': '#ac6036', "SVDSS": "#9bd4d7", "PAV": "#cad886", "SVIM-asm": "#d88686", "debreak": "#579979"}

    # hg002_read_hifi_high_confident = data.hg002_read_hifi_high_confident
    # hg002_read_ont_high_confident = data.hg002_read_ont_high_confident
    # hg002_asm_hifi_high_confident = data.hg002_asm_hifi_high_confident
    # hg002_asm_ont_high_confident = data.hg002_asm_ont_high_confident
    #
    mode = "read"
    # # mode = "assembly"
    #
    # plot_hg002_high_confident_precision_recall(hg002_read_hifi_high_confident, hg002_read_ont_high_confident, mode=mode)
    #
    # plot_hg002_high_confident_f1_gt(hg002_read_hifi_high_confident, hg002_read_ont_high_confident, mode=mode)

    hg002_hifi_high_confident = {
        "SVision-pro": [0.9663, 0.9689, 0.9676, 0.9857],
        "SVision": [0.923, 0.9591, 0.9407, 0.9741],
        "Sniffles2": [0.9447, 0.9654, 0.955, 0.978],
        "cuteSV": [0.9607, 0.9545, 0.9576, 0.9885],
        "debreak": [0.9597, 0.971, 0.9653, 0.9634],
        "pbsv": [0.9296, 0.9537, 0.9415, 0.9765],
        "SVDSS": ["NA", "NA", "NA", 'NA']
    }

    hg002_ont_high_confident = {
        "SVision-pro": [0.9509, 0.9406, 0.9457, 0.9372],
        "SVision": [0.8065, 0.9187, 0.8589, 0.8735],
        "Sniffles2": [0.9237, 0.9304, 0.927, 0.9575],
        "cuteSV": [0.9338, 0.9022, 0.9178, 0.988],
        "debreak": [0.9478, 0.9474, 0.9476, 0.9379],
        "pbsv": [0.7366, 0.876, 0.8003, 0.9187],
        "SVDSS": ["NA", "NA", "NA", 'NA']
    }

    csv_hifi_simulated = {
        "SVision-pro": [0.9841, 0.951, 0.9672, 0.94],
        "SVision": [0.8391, 0.9456, 0.8892, 0.8759],
        "Sniffles2": [0.418, 0.374, 0.3947, 0.237],
        "cuteSV": [0.5308, 0.459, 0.4923, 0.7879],
        "debreak": [0.5203, 0.7553, 0.6161, 0.139],
    }

    csv_ont_simulated = {
        "SVision-pro": [0.9881, 0.9483, 0.9678, 0.8564],
        "SVision": [0.8027, 0.898, 0.8477, 0.693],
        "Sniffles2": [0.4716, 0.582, 0.521, 0.2256],
        "cuteSV": [0.5145, 0.46, 0.4857, 0.692],
        "debreak": [0.5314, 0.7213, 0.6119, 0.1238],
    }


    csv_structure = {
        "SVision-pro": [0.985309549, 0.96173913, 0.985278654, 0.964850615],
        "SVision": [0.827513966, 0.682678822, 0.91716602, 0.862657758],
        "Sniffles2": [0.258178603, 0.161904762, 0.254901961, 0.159793814],
        "cuteSV": [0.182591623, 0.167369055, 0.179375454, 0.166666667],
        "debreak": [0.118959108, 0.112253289, 0.112974404, 0.091497227],
    }


    # plot_hg002_high_confident_precision_recall(hg002_hifi_high_confident, hg002_ont_high_confident, mode)
    # plot_hg002_high_confident_f1_gt(hg002_hifi_high_confident, hg002_ont_high_confident, mode)

    # plot_csv_simulated_precision_recall(csv_hifi_simulated, csv_ont_simulated, mode)
    # plot_csv_simualted_f1_gt(csv_hifi_simulated, csv_ont_simulated, mode)

    plot_csv_structure_boxplot(csv_structure)