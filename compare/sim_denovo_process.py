import os
import pysam

def determine_mendelian_error(child_gt, mother_gt, father_gt):

    if father_gt == "1/1" and mother_gt == "1/1" and child_gt != "1/1":
        return "bad"

    if father_gt == "0/0" and mother_gt == "0/0" and child_gt != "0/0":
        return "bad"

    if ((father_gt == "1/1" and mother_gt == "0/0") or (father_gt == "0/0" and mother_gt == "1/1")) and child_gt != "0/1":
        return "bad"

    if ((father_gt == "1/1" and mother_gt == "0/1") or (father_gt == "0/1" and mother_gt == "1/1")) and child_gt == "0/0":
        return "bad"

    if ((father_gt == "0/1" and mother_gt == "0/0") or (father_gt == "0/0" and mother_gt == "0/1")) and child_gt == "1/1":
        return "bad"

    return "ok"



def add_trio_gt_to_denovo_benchmark_vcf(benchmark_vcf, out_vcf, denovo_log, sv_class="SSV"):

    out_vcf = open(out_vcf, "w")

    benchmark_vcf = pysam.VariantFile(benchmark_vcf)

    out_vcf.write(str(benchmark_vcf.header))

    for record in benchmark_vcf:
        record_chrom = record.contig
        record_start = record.start + 1
        record_svtype = record.info["SVTYPE"]

        if "dDUP" == record_svtype or "idDUP" == record_svtype:
            search_info = "{}_{}".format(record_chrom, record.info["BKPS"][0].split("_")[3])

        elif "dDUP" in record_svtype :
            all_bkps = []

            for bkp_str in record.info["BKPS"]:
                bkp_str_split= bkp_str.split("_")
                all_bkps.append(int(bkp_str_split[3]))
                all_bkps.append(int(bkp_str_split[4]))

            search_info = "{}_{}".format(record_chrom, min(all_bkps))
        else:
            search_info = "{}_{}".format(record_chrom, record_start)

        found_flag = False

        # # search in log file
        for line in open(denovo_log):
            line_split = line.strip().split("\t")

            target_info = "{}_{}".format(line_split[0], line_split[1])

            if search_info == target_info:

                if sv_class == "SSV":
                    father_allele = float(line.strip().split("\t")[6].split(":")[1].split("-")[0])
                    mother_allele = float(line.strip().split("\t")[6].split(":")[1].split("-")[1])
                    child_allele = float(line.strip().split("\t")[6].split(":")[1].split("-")[2])
                else:
                    father_allele = float(line.strip().split("\t")[5].split(":")[1].split("-")[0])
                    mother_allele = float(line.strip().split("\t")[5].split(":")[1].split("-")[1])
                    child_allele = float(line.strip().split("\t")[5].split(":")[1].split("-")[2])


                if child_allele == 0.0:
                    continue

                if father_allele == 0:
                    father_gt = "0/0"
                elif father_allele == 0.5:
                    father_gt = "0/1"
                elif father_allele == 1:
                    father_gt = "1/1"
                else:
                    print("No this allele", father_allele)
                    exit()

                if mother_allele == 0:
                    mother_gt = "0/0"
                elif mother_allele == 0.5:
                    mother_gt = "0/1"
                elif mother_allele == 1:
                    mother_gt = "1/1"
                else:
                    print("No this allele", mother_allele)
                    exit()

                if child_allele == 0:
                    child_gt = "0/0"
                elif child_allele == 0.5:
                    child_gt = "0/1"
                elif child_allele == 1:
                    child_gt = "1/1"
                else:
                    print("No this allele", child_allele)
                    exit()

                found_flag = True

                itbkp_str = []

                for simple_sv_type in record_svtype.split("+"):
                    itbkp_str.append("{}_{}_{}_{}".format(simple_sv_type, child_gt, father_gt, mother_gt))

                itbkp_str = ",".join(itbkp_str)
                break

        if not found_flag:
            print(search_info)

        else:
            record_split = str(record).split("\t")

            record_split[7] = "Trio={};".format(itbkp_str) + record_split[7]

            out_vcf.write("\t".join(record_split))


def evaluate_trio_gt_concordance_from_truvari_res(truvari_folder, mode="denovo"):
    # # calculate component concordance

    base_record_list = []
    base_trio_gt_list = []
    call_trio_gt_list = []

    # collect from tp base
    for record in pysam.VariantFile(os.path.join(truvari_folder, "tp-base.vcf")):
        base_trio_gt_list.append("_".join(record.info["Trio"][0].split("_")[1: ]))

        base_record_list.append("{}_{}_{}".format(record.contig, record.start + 1, record.info["SVTYPE"]))

    # collect from tp call
    for record in pysam.VariantFile(os.path.join(truvari_folder, "tp-call.vcf")):

        father_gt = str(record).strip().split("\t")[-2].split(":")[0]
        mother_gt = str(record).strip().split("\t")[-1].split(":")[0]
        child_gt = str(record).strip().split("\t")[-3].split(":")[0]

        if "survivor" in truvari_folder or "jasmine" in truvari_folder:

            # convert_dict = {"1/1": "1/1", "0/1": "0/1", "0/0": "./.", "./.": "0/0"}
            convert_dict = {"1/1": "1/1", "0/1": "0/1", "0/0": "0/0", "./.": "0/0"}

            father_gt = convert_dict[father_gt]
            mother_gt = convert_dict[mother_gt]
            child_gt = convert_dict[child_gt]


        call_trio_gt_list.append("{}_{}_{}".format(child_gt, father_gt, mother_gt))


    # calculate concordance
    denovo_tot_events = 0
    denovo_match_events = 0

    for i in range(len(base_trio_gt_list)):

        # # skpi denovo event
        if not base_trio_gt_list[i].endswith("_0/0_0/0"):
            continue

        denovo_tot_events += 1

        if call_trio_gt_list[i] == base_trio_gt_list[i]:
            denovo_match_events += 1

        # else:
        #     print(base_record_list[i], base_trio_gt_list[i], call_trio_gt_list[i])
    # print("denovo", denovo_tot_events, denovo_match_events, round(denovo_match_events / denovo_tot_events, 3))
    # print(denovo_tot_events)
    # print(denovo_match_events)
    # print(round(denovo_match_events / denovo_tot_events, 3))

    inherited_tot_events = 0
    inherited_match_events = 0
    inherited_mendelian_dict = {"total": 0, "ok": 0, "bad": 0}

    for i in range(len(base_trio_gt_list)):

        # # skpi denovo event
        if base_trio_gt_list[i].endswith("_0/0_0/0"):
            continue

        inherited_tot_events += 1

        if call_trio_gt_list[i] == base_trio_gt_list[i]:
            inherited_match_events += 1
        # else:
        #     print(base_record_list[i], base_trio_gt_list[i], call_trio_gt_list[i])

        inherited_mendelian_dict[determine_mendelian_error(call_trio_gt_list[i].split("_")[0], call_trio_gt_list[i].split("_")[1], call_trio_gt_list[i].split("_")[2])] += 1
    # print(inherited_mendelian_dict)

    # print("inherited", inherited_tot_events, inherited_match_events, round(inherited_match_events / inherited_tot_events, 3))
    # print(inherited_tot_events)
    # print(inherited_match_events)
    # print(round(inherited_match_events / inherited_tot_events, 3))


    match_events = denovo_match_events + inherited_match_events
    total_events = denovo_tot_events + inherited_tot_events

    print(match_events, total_events - match_events, match_events / total_events, truvari_folder)
    print(denovo_match_events / denovo_tot_events, inherited_match_events / inherited_tot_events)


if __name__ == '__main__':
    # add_trio_gt_to_denovo_benchmark_vcf("/mnt/c/workspace/test_new/sim_denovo_ccs/ground_truth/HG002_SVs_Tier1_v0.6.final.vcf",
    #                                     "/mnt/c/workspace/test_new/sim_denovo_ccs/ground_truth/HG002_SVs_Tier1_v0.6.final.addTrio.vcf",
    #                                     "/mnt/c/workspace/test_new/sim_denovo_ccs/ground_truth/HG002_SVs_Tier1_v0.6.final.denovo_allele.log",
    #                                     sv_class="SSV")


    # evaluate_trio_gt_concordance_from_truvari_res("/mnt/c/workspace/test_new/sim_denovo_ccs/complex_old/svision_new.s10.survivor_bench_region")

    print("--------------------hifi")
    evaluate_trio_gt_concordance_from_truvari_res("/mnt/c/workspace/test_new/sim_denovo_ccs/complex/release_256.svision_pro_v1.6.s10_bench_region")

    evaluate_trio_gt_concordance_from_truvari_res("/mnt/c/workspace/test_new/sim_denovo_ccs/complex/svision_survivor.child-father-mother.s10_bench_region")
    evaluate_trio_gt_concordance_from_truvari_res("/mnt/c/workspace/test_new/sim_denovo_ccs/complex/svision_jasmine.child-father-mother.s10_bench_region")


    print("--------------------ont")

    evaluate_trio_gt_concordance_from_truvari_res("/mnt/c/workspace/test_new/sim_denovo_ont/complex/release_256.svision_pro_v1.6.s10_bench_region")
    evaluate_trio_gt_concordance_from_truvari_res("/mnt/c/workspace/test_new/sim_denovo_ont/complex/svision_survivor.child-father-mother.s10_bench_region")
    evaluate_trio_gt_concordance_from_truvari_res("/mnt/c/workspace/test_new/sim_denovo_ont/complex/svision_jasmine.child-father-mother.s10_bench_region")


