import os
import pysam
import numpy as np


def compare_csv_subcomponents(truvari_res_folder):

    base_csv_list = []
    call_csv_list = []

    # # collect from base
    for record in pysam.VariantFile(os.path.join(truvari_res_folder, "tp-base.vcf")):

        csv_info = {}

        csv_type = record.info["SVTYPE"]

        for i in range(len(csv_type.split("+"))):
            ssv_type = csv_type.split("+")[i]
            ssv_bkp = record.info["BKPS"][i]

            if ssv_type in ["dDUP", "idDUP", "tDUP", "itDUP"]:
                ssv_type = "DUP"

            ssv_start = int(ssv_bkp.split("_")[3])

            ssv_end = int(ssv_bkp.split("_")[4])

            csv_info[ssv_type] = [ssv_start, ssv_end]

        base_csv_list.append(csv_info)

    # # collect from target
    # # collect from base
    for record in pysam.VariantFile(os.path.join(truvari_res_folder, "tp-call.vcf")):

        csv_info = {}

        csv_type = record.info["SVTYPE"]

        for i in range(len(csv_type.split("+"))):
            ssv_type = csv_type.split("+")[i]

            if ssv_type in ["dDUP", "idDUP", "tDUP", "itDUP"]:
                ssv_type = "DUP"

            if "svision_pro" in truvari_res_folder:
                ssv_bkp = record.info["BKPS"][i]

                ssv_start = int(ssv_bkp.split("_")[3])

                ssv_end = int(ssv_bkp.split("_")[4])
            elif "svision" in truvari_res_folder:
                ssv_bkp = record.info["BKPS"][i]

                ssv_start = int(ssv_bkp.split("-")[1])

                ssv_end = int(ssv_bkp.split("-")[2])

            else:
                ssv_start = record.start + 1
                ssv_end = record.stop

            csv_info[ssv_type] = [ssv_start, ssv_end]

        call_csv_list.append(csv_info)


    total_event = 0
    match_event = 0
    component_miss_event = 0
    component_disorder_event = 0

    bkp_diff = []
    for i in range(len(base_csv_list)):
        total_event += 1

        base_csv_info = base_csv_list[i]
        call_csv_info = call_csv_list[i]

        component_miss_flag = False
        for base_ssv in base_csv_info:

            if base_ssv not in call_csv_info:
                component_miss_flag = True
            else:
                base_ssv_length = abs(base_csv_info[base_ssv][1] - base_csv_info[base_ssv][0])
                call_ssv_length = abs(call_csv_info[base_ssv][1] - call_csv_info[base_ssv][0])

                # if abs(base_ssv_length - call_ssv_length) > 1 :
                bkp_diff.append(abs(base_ssv_length - call_ssv_length))

        if component_miss_flag is False:
            # if list(base_csv_info.keys()) != list(call_csv_info.keys()):
            #     component_disorder_event += 1
            #     print(base_csv_info, call_csv_info)
            # else:
                match_event += 1
        else:
            component_miss_event += 1

    print(truvari_res_folder)
    print(total_event)
    print(match_event)
    # print(component_disorder_event)
    print(component_miss_event)

    print(np.average(bkp_diff), np.std(bkp_diff))

    bkp_collect = [[], [], [], []]

    for i in bkp_diff:
        # print(i)
        if i < 10:
            bkp_collect[0].append(i)
        elif i < 100:
            bkp_collect[1].append(i)
        elif i < 500:
            bkp_collect[2].append(i)
        else:
            bkp_collect[3].append(i)

    print("bkp shift")
    for i in range(len(bkp_collect)):
        bkp_collect[i] = len(bkp_collect[i])
        print(bkp_collect[i])


if __name__ == '__main__':

    compare_csv_subcomponents("/mnt/c/workspace/test_new/sim_csv_ccs/release.svision_pro_v1.6.s5_bench_region")

    #
    # compare_csv_subcomponents("/mnt/c/workspace/test_new/sim_csv_ont/sim_csv.svision.s10_bench_region")
    #
    # compare_csv_subcomponents("/mnt/c/workspace/test_new/sim_csv_ont/sim_csv.sniffles2.s5.m50_bench_region")
    #
    # compare_csv_subcomponents("/mnt/c/workspace/test_new/sim_csv_ccs/sim_csv.cutesv.s10.m50_bench_region")
    #
    #
    # compare_csv_subcomponents("/mnt/c/workspace/test_new/sim_csv_ont/sim_csv.debreak.s10.m50_bench_region")
