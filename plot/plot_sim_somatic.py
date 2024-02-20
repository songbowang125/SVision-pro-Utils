
from matplotlib import pyplot as plt
import math
import os
import numpy as np
import scipy as sc
import pylab as pl
import seaborn as sns
from scipy.stats import ranksums
import random

def autolabel(rects, ax):
    for rect in rects:
        height = rect.get_height()
        ax.annotate('%.2f' % height,
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3点垂直偏移
                    textcoords="offset points",
                    ha='center', va='bottom')


def plot_sim_somatic_accuracy(val_dict, seq='hifi'):

    if seq == "hifi":
        marker = "."
        line = "-"
    else:
        marker = '.'
        line = "--"

    # af_total = np.array([1090, 1088, 1099, 1058, 1086, 1045, 1024, 1055, 1096])
    # af_legend = [0.25, 0.2, 0.15, 0.10, 0.08, 0.05, 0.03, 0.02, 0.01]
    # af_legend = [0.10, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01]
    af_legend = [0.10,  0.08, 0.05, 0.04, 0.03, 0.02, 0.01]

    # # STEP: plot
    r1 = np.arange(len(af_legend))

    fig = pl.figure(figsize=(2.5, 3))  # llc_x, llc_y, width, height
    # fig = pl.figure(figsize=(2.5, 1))  # llc_x, llc_y, width, height
    # fig = pl.figure(figsize=(2.5, 1.5))  # llc_x, llc_y, width, height

    ax = fig.add_subplot(1, 1, 1)
    svision_pro_256_na = np.isfinite(np.array(val_dict["SVision-pro-256"]).astype(np.double))
    svision_pro_512_na = np.isfinite(np.array(val_dict["SVision-pro-512"]).astype(np.double))
    svision_pro_1024_na = np.isfinite(np.array(val_dict["SVision-pro-1024"]).astype(np.double))
    sniffles2_na = np.isfinite(np.array(val_dict["sniffles2"]).astype(np.double))
    nanomonsv_na = np.isfinite(np.array(val_dict["nanomonsv"]).astype(np.double))

    line1 = ax.plot(r1[svision_pro_256_na], np.array(val_dict["SVision-pro-256"])[svision_pro_256_na], label="SVision-pro-256", color="#4386bc", marker=marker, ls=line)
    line2 = ax.plot(r1[svision_pro_512_na], np.array(val_dict["SVision-pro-512"])[svision_pro_512_na], label="SVision-pro-512", color="#579979", marker=marker, ls=line)
    line3 = ax.plot(r1[svision_pro_1024_na], np.array(val_dict["SVision-pro-1024"])[svision_pro_1024_na], label="SVision-pro-1024", color="#ebbca5", marker=marker, ls=line)
    line4 = ax.plot(r1[sniffles2_na], np.array(val_dict["sniffles2"])[sniffles2_na], label="sniffles2", color="#a0aebf", marker=marker, ls=line)
    line5 = ax.plot(r1[nanomonsv_na], np.array(val_dict["nanomonsv"])[nanomonsv_na], label="nanomonsv", color="#abc7de", marker=marker, ls=line)


    ax.grid(b=True, color='grey', linestyle='-.', linewidth=0.5, alpha=0.2)
    # ax.legend(loc='best', fontsize=7,)

    # Remove axes splines
    for s in [ 'left', 'right']:
        ax.spines[s].set_visible(False)

    plt.xticks([r for r in range(len(af_legend))], ["{}".format(af) for af in af_legend], fontsize=8,)

    # # for SSV
    ax.set_ylim([0.65, 1.02])  # xmin, xmax, ymin, ymax
    ax.set_yticks(np.linspace(0.7, 1.0, 4))
    ax.set_yticklabels([], fontsize=8)
    ax.set_yticklabels([round(val, 2) for val in np.linspace(0.7, 1.0, 4)], fontsize=8)

    # # for SSV ONT
    ax.set_ylim([0.65, 1.02])  # xmin, xmax, ymin, ymax
    ax.set_yticks(np.linspace(0.7, 1.0, 4))
    ax.set_yticklabels([], fontsize=8)
    # ax.set_yticklabels([round(val, 2) for val in np.linspace(0.7, 1.0, 4)], fontsize=8)

    # # for CSV hifi
    # ax.set_ylim([0.9, 1.02])  # xmin, xmax, ymin, ymax
    # ax.set_yticks(np.linspace(0.9, 1.0, 2))
    # ax.set_yticklabels([], fontsize=8)
    # ax.set_yticklabels([round(val, 2) for val in np.linspace(0.9, 1.0, 2)], fontsize=8)
    # plt.xticks([r for r in range(len(af_legend))], [], fontsize=8,)

    # # # for CSV ont
    # ax.set_ylim([0.8, 1.02])  # xmin, xmax, ymin, ymax
    # ax.set_yticks(np.linspace(0.8, 1.0, 3))
    # ax.set_yticklabels([], fontsize=8)
    # ax.set_yticklabels([round(val, 2) for val in np.linspace(0.8, 1.0, 3)], fontsize=8)

    # plt.ylabel("SSV HiFi")

    plt.savefig("out_sim_somatic_accuracy.png", dpi=200, bbox_inches='tight',)
    # plt.savefig("out_sim_somatic_accuracy.svg")



def plot_sim_somatic_accuracy_1024(val_dict, val_dict_ont):


    # af_legend = [0.10,  0.08, 0.05, 0.04, 0.03, 0.02, 0.01]
    af_legend = [ 0.04, 0.03, 0.02, 0.01]

    for key in val_dict:
        val_dict[key][-4] = val_dict[key][0]
        val_dict[key][-3] = val_dict[key][3]

        val_dict[key] = val_dict[key][-4: ]

    for key in val_dict_ont:
        val_dict[key][-4] = val_dict[key][0]
        val_dict[key][-3] = val_dict[key][3]

        val_dict_ont[key] = val_dict_ont[key][-4: ]

    # # STEP: plot
    r1 = np.arange(len(af_legend))

    # fig = pl.figure(figsize=(2.5, 3))  # llc_x, llc_y, width, height
    fig = pl.figure(figsize=(1.2, 3))  # llc_x, llc_y, width, height


    ax = fig.add_subplot(1, 1, 1)

    marker = "."
    line = "-"
    svision_pro_1024_na = np.isfinite(np.array(val_dict["SVision-pro-1024"]).astype(np.double))
    sniffles2_na = np.isfinite(np.array(val_dict["sniffles2"]).astype(np.double))
    nanomonsv_na = np.isfinite(np.array(val_dict["nanomonsv"]).astype(np.double))

    ax.plot(r1[svision_pro_1024_na], np.array(val_dict["SVision-pro-1024"])[svision_pro_1024_na], label="SVision-pro", color="#ebbca5", marker=marker, ls=line)
    ax.plot(r1[sniffles2_na], np.array(val_dict["sniffles2"])[sniffles2_na], label="sniffles2", color="#a0aebf", marker=marker, ls=line)
    ax.plot(r1[nanomonsv_na], np.array(val_dict["nanomonsv"])[nanomonsv_na], label="nanomonsv", color="#abc7de", marker=marker, ls=line)

    marker = '.'
    line = "--"
    svision_pro_1024_na_ont = np.isfinite(np.array(val_dict_ont["SVision-pro-1024"]).astype(np.double))
    sniffles2_na_ont = np.isfinite(np.array(val_dict_ont["sniffles2"]).astype(np.double))
    nanomonsv_na_ont = np.isfinite(np.array(val_dict_ont["nanomonsv"]).astype(np.double))

    ax.plot(r1[svision_pro_1024_na_ont], np.array(val_dict_ont["SVision-pro-1024"])[svision_pro_1024_na], label="SVision-pro", color="#ebbca5", marker=marker, ls=line)
    ax.plot(r1[sniffles2_na_ont], np.array(val_dict_ont["sniffles2"])[sniffles2_na], label="sniffles2", color="#a0aebf", marker=marker, ls=line)
    ax.plot(r1[nanomonsv_na_ont], np.array(val_dict_ont["nanomonsv"])[nanomonsv_na], label="nanomonsv", color="#abc7de", marker=marker, ls=line)

    # ax.grid(b=True, color='grey', linestyle='-.', linewidth=0.5, alpha=0.2)
    ax.grid(color='grey', linestyle='-.', linewidth=0.5, alpha=0.2)

    # ax.legend(loc='best', fontsize=7,)

    # Remove axes splines
    for s in [ 'left', 'right']:
        ax.spines[s].set_visible(False)

    plt.xticks([r for r in range(len(af_legend))], ["0.10", "0.04", "0.02", "0.01"], fontsize=8,rotation=340)

    # # for SSV
    ax.set_ylim([0.55, 1.02])  # xmin, xmax, ymin, ymax
    ax.set_yticks(np.linspace(0.6, 1.0, 4))
    ax.set_yticklabels([], fontsize=10)
    ax.set_yticklabels([round(val, 2) for val in np.linspace(0.6, 1.0, 4)], fontsize=8)

    # # for CSV
    # ax.set_ylim([0.9, 1.02])  # xmin, xmax, ymin, ymax
    # ax.set_yticks(np.linspace(0.9, 1.0, 2))
    # ax.set_yticklabels([], fontsize=8)
    # ax.set_yticklabels([round(val, 2) for val in np.linspace(0.9, 1.0, 2)], fontsize=8)
    # plt.xticks([r for r in range(len(af_legend))], [], fontsize=8,)

    # # # for CSV ont
    # ax.set_ylim([0.8, 1.02])  # xmin, xmax, ymin, ymax
    # ax.set_yticks(np.linspace(0.8, 1.0, 3))
    # ax.set_yticklabels([], fontsize=8)
    # ax.set_yticklabels([round(val, 2) for val in np.linspace(0.8, 1.0, 3)], fontsize=8)

    # plt.ylabel("SSV HiFi")
    # plt.ylabel("Somatic SSV accuracy", fontsize=16)

    plt.savefig("out_sim_somatic_accuracy.png", dpi=1500, bbox_inches='tight',)
    # plt.savefig("out_sim_somatic_accuracy.svg")


def plot_sim_somatic_accuracy_1024_csv(val_dict, val_dict_ont):

    # af_total = np.array([1090, 1088, 1099, 1058, 1086, 1045, 1024, 1055, 1096])
    # af_legend = [0.25, 0.2, 0.15, 0.10, 0.08, 0.05, 0.03, 0.02, 0.01]
    # af_legend = [0.10, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01]
    af_legend = [0.04, 0.03, 0.02, 0.01]

    for key in val_dict:
        val_dict[key][-4] = val_dict[key][0]
        val_dict[key][-3] = val_dict[key][-4]

        val_dict[key] = val_dict[key][-4: ]



    for key in val_dict_ont:
        val_dict[key][-4] = val_dict[key][0]
        val_dict[key][-3] = val_dict[key][-4]

        val_dict_ont[key] = val_dict_ont[key][-4: ]

    # # STEP: plot
    r1 = np.arange(len(af_legend))

    fig = pl.figure(figsize=(1.2, 3))  # llc_x, llc_y, width, height
    # fig = pl.figure(figsize=(2.5, 1))  # llc_x, llc_y, width, height
    # fig = pl.figure(figsize=(2.5, 1.5))  # llc_x, llc_y, width, height

    ax = fig.add_subplot(1, 1, 1)

    marker = "."
    line = "-"
    svision_pro_1024_na = np.isfinite(np.array(val_dict["SVision-pro-1024"]).astype(np.double))
    sniffles2_na = np.isfinite(np.array(val_dict["sniffles2"]).astype(np.double))
    nanomonsv_na = np.isfinite(np.array(val_dict["nanomonsv"]).astype(np.double))

    ax.plot(r1[svision_pro_1024_na], np.array(val_dict["SVision-pro-1024"])[svision_pro_1024_na], label="SVision-pro", color="#ebbca5", marker=marker, ls=line)
    ax.plot(r1[sniffles2_na], np.array(val_dict["sniffles2"])[sniffles2_na], label="sniffles2", color="#a0aebf", marker=marker, ls=line)
    ax.plot(r1[nanomonsv_na], np.array(val_dict["nanomonsv"])[nanomonsv_na], label="nanomonsv", color="#abc7de", marker=marker, ls=line)

    marker = '.'
    line = "--"
    svision_pro_1024_na_ont = np.isfinite(np.array(val_dict_ont["SVision-pro-1024"]).astype(np.double))
    sniffles2_na_ont = np.isfinite(np.array(val_dict_ont["sniffles2"]).astype(np.double))
    nanomonsv_na_ont = np.isfinite(np.array(val_dict_ont["nanomonsv"]).astype(np.double))

    ax.plot(r1[svision_pro_1024_na_ont], np.array(val_dict_ont["SVision-pro-1024"])[svision_pro_1024_na], label="SVision-pro", color="#ebbca5", marker=marker, ls=line)
    ax.plot(r1[sniffles2_na_ont], np.array(val_dict_ont["sniffles2"])[sniffles2_na], label="sniffles2", color="#a0aebf", marker=marker, ls=line)
    ax.plot(r1[nanomonsv_na_ont], np.array(val_dict_ont["nanomonsv"])[nanomonsv_na], label="nanomonsv", color="#abc7de", marker=marker, ls=line)

    # ax.grid(b=True, color='grey', linestyle='-.', linewidth=0.5, alpha=0.2)
    ax.grid( color='grey', linestyle='-.', linewidth=0.5, alpha=0.2)

    # ax.legend(loc='best', fontsize=7,)

    # Remove axes splines
    for s in [ 'left', 'right']:
        ax.spines[s].set_visible(False)

    plt.xticks([r for r in range(len(af_legend))], ["0.10", "0.04", "0.02", "0.01"], fontsize=8,rotation=340)

    ax.set_ylim([0.65, 1.02])  # xmin, xmax, ymin, ymax
    ax.set_yticks(np.linspace(0.7, 1.0, 4))
    ax.set_yticklabels([], fontsize=8)
    # ax.set_yticklabels([round(val, 2) for val in np.linspace(0.7, 1.0, 4)], fontsize=8)

    # # for CSV
    # ax.set_ylim([0.9, 1.02])  # xmin, xmax, ymin, ymax
    # ax.set_yticks(np.linspace(0.9, 1.0, 2))
    # ax.set_yticklabels([], fontsize=8)
    # ax.set_yticklabels([round(val, 2) for val in np.linspace(0.9, 1.0, 2)], fontsize=8)
    # plt.xticks([r for r in range(len(af_legend))], [], fontsize=8,)

    # # # for CSV ont
    # ax.set_ylim([0.8, 1.02])  # xmin, xmax, ymin, ymax
    # ax.set_yticks(np.linspace(0.8, 1.0, 3))
    # ax.set_yticklabels([], fontsize=8)
    # ax.set_yticklabels([round(val, 2) for val in np.linspace(0.8, 1.0, 3)], fontsize=8)

    # plt.ylabel("SSV HiFi")

    plt.savefig("out_sim_somatic_accuracy.png", dpi=1500, bbox_inches='tight',)
    # plt.savefig("out_sim_somatic_accuracy.svg")


def plot_sim_somatic_precision_recall():
    import matplotlib.pyplot as plt

    # 示例数据
    thresholds = [0.2, 0.4, 0.6, 0.8]  # 示例阈值
    precision_values = [0.6, 0.7, 0.8, 0.75]  # 示例精确率值
    recall_values = [0.5, 0.6, 0.7, 0.65]  # 示例召回率值

    # 计算F1值
    f1_values = [2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
                 for precision, recall in zip(precision_values, recall_values)]

    # 绘图
    fig, ax1 = plt.subplots()

    # 绘制精确率和召回率的柱状图
    ax1.bar(thresholds, precision_values, width=0.1, color='b', align='center', label='Precision')
    ax1.bar([t + 0.1 for t in thresholds], recall_values, width=0.1, color='g', align='center', label='Recall')
    ax1.set_xlabel('Thresholds')
    ax1.set_ylabel('Precision/Recall')
    ax1.legend(loc='upper left')

    # 添加第二个y轴并绘制F1的折线图
    ax2 = ax1.twinx()
    ax2.plot(thresholds, f1_values, marker='o', linestyle='-', color='r', label='F1')
    ax2.set_ylabel('F1', color='r')
    ax2.tick_params(axis='y', labelcolor='r')
    ax2.legend(loc='upper right')

    plt.title('Precision, Recall, and F1')
    plt.show()

if __name__ == '__main__':
    # hifi_accurate = data.sim_somatic_hifi_accurate_values
    # ont_accurate = data.sim_somatic_ont_accurate_values
    #
    # plot_sim_somatic_accuracy(hifi_accurate, ont_accurate)

    #
    # hifi_accurate_runtime = data.sim_somatic_hifi_allaccurate_runtime
    # ont_accurate_runtime = data.sim_somatic_ont_allaccurate_runtime
    #
    # plot_sim_somatic_allaccuracy_runtime(hifi_accurate_runtime, ont_accurate_runtime)



    sim_ssv_hifi_val_dict = {"SVision-pro-256": [None, None, None, None, None, None, None],
                     "SVision-pro-512": [None, None, None, None, None, None, None],
                     "SVision-pro-1024": [0.982, 0.981, 0.977, 0.979, 0.973, 0.971, 0.965],
                     "sniffles2": [0.938, 0.929, 0.934, 0.938, 0.921, 0.933, 0.919],
                     "nanomonsv": [0.728, 0.715, 0.702, 0.705, 0.7, 0.72, 0.687]}

    sim_ssv_ont_val_dict = {"SVision-pro-256": [None, None, None, None, None, None, None],
                    "SVision-pro-512": [None, None, None, None, None, None, None],
                    "SVision-pro-1024": [0.989, 0.992, 0.981, 0.981, 0.974, 0.977, 0.889],
                    "sniffles2": [0.907, 0.925, 0.916, 0.916, 0.924, 0.927, 0.843],
                    "nanomonsv": [0.966, 0.966, 0.967, 0.958, 0.945, 0.931, 0.816]}

    sim_csv_hifi_val_dict = {
        "SVision-pro-256": [None, None, None, None, None, None, None],
        "SVision-pro-512": [None, None, None, None, None, None, None],
        "SVision-pro-1024": [0.959, 0.953, 0.962, 0.969, 0.958, 0.952, 0.941],
        "sniffles2": [None, None, None, None, None, None, None],
        "nanomonsv": [None, None, None, None, None, None, None],
    }

    sim_csv_ont_val_dict = {
        "SVision-pro-256": [None, None, None, None, None, None, None],
        "SVision-pro-512": [None, None, None, None, None, None, None],
        "SVision-pro-1024": [0.983, 0.96, 0.959, 0.954, 0.936, 0.947, 0.92],
        "sniffles2": [None, None, None, None, None, None, None],
        "nanomonsv": [None, None, None, None, None, None, None],
    }
    # plot_sim_somatic_accuracy(sim_ssv_ont_val_dict, "ont")
    # plot_sim_somatic_accuracy_1024(sim_ssv_hifi_val_dict, sim_ssv_ont_val_dict)
    plot_sim_somatic_accuracy_1024_csv(sim_csv_hifi_val_dict, sim_csv_ont_val_dict)

    # plot_sim_somatic_round2("/mnt/c/workspace/test/sim_somatic_ccs/HG002_SVs_Tier1_v0.6_high_confident.bamsurgeon.round2.vcf", "/mnt/c/workspace/test/sim_somatic_ont/sv_call_round2/svision-pro-256_bench_region")
