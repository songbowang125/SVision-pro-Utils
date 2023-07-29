import os.path

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


def calculate_quartet_twin(vcf_path, min_support=10, exist_in_child=True, convert_gt=True, output_discordant=False):

    cnt_dict = {"total": 0, "identical": 0, "not_identical": 0}

    in_vcf = pysam.VariantFile(vcf_path)

    if output_discordant is True:
        discordant_file = open(vcf_path.replace(".vcf", ".twin_discordant.vcf"), "w")
        cordant_file = open(vcf_path.replace(".vcf", ".twin_cordant.vcf"), "w")

        discordant_file.write(str(in_vcf.header))
        cordant_file.write(str(in_vcf.header))

    for record in in_vcf:

        line = str(record)

        line_split = line.strip().split("\t")

        # # STEP: only keep autosome
        if record.contig not in ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", ]:
            continue

        if "PASS" not in line:
            continue

        # # STEP: must exist in child sample
        if exist_in_child:
            if "svision_pro" not in vcf_path:
                if "SUPP_VEC=1" not in str(record):
                    continue

        if "BND" in line:
            continue

        # # STEP: extract info
        child_gt_info = line_split[-4]
        child2_gt_info = line_split[-3]
        mother_gt_info = line_split[-2]
        father_gt_info = line_split[-1]

        child_gt = child_gt_info.split(":")[0]
        child2_gt = child2_gt_info.split(":")[0]
        mother_gt = mother_gt_info.split(":")[0]
        father_gt = father_gt_info.split(":")[0]

        cnt_dict["total"] += 1


        if output_discordant is False:
            if child_gt != child2_gt:
                cnt_dict["not_identical"] += 1
                # print(record, end="")

            else:
                cnt_dict["identical"] += 1

        else:
            if child_gt != child2_gt:
                cnt_dict["not_identical"] += 1
                discordant_file.write(str(record))

            else:
                cnt_dict["identical"] += 1
                cordant_file.write(str(record))

    print(cnt_dict["not_identical"] / cnt_dict["total"], cnt_dict, vcf_path)


def calculate_mendelian(vcf_path, min_support=10, exist_in_child=True, convert_gt=True):

    cnt_dict = {"total": 0, "ok": 0, "bad": 0}

    for record in pysam.VariantFile(vcf_path):

        line = str(record)

        if "PASS" not in str(record):
            continue

        line_split = line.strip().split("\t")

        # # STEP: only keep autosome
        if record.contig not in ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", ]:
            continue

        # # STEP: must exist in child sample
        if exist_in_child:
            if "svision_pro" not in vcf_path:
                if "SUPP_VEC=1" not in str(record):
                    continue

        # # STEP: extract info
        child_gt_info = line_split[-3]
        mother_gt_info = line_split[-2]
        father_gt_info = line_split[-1]

        child_gt = child_gt_info.split(":")[0]
        mother_gt = mother_gt_info.split(":")[0]
        father_gt = father_gt_info.split(":")[0]

        # #  STEP: filter by support fix sniffles2 min support bug
        if "sniffles2" in vcf_path:
            child_support = int(child_gt_info.split(":")[-2])
            mother_support = int(mother_gt_info.split(":")[-2])
            father_support = int(father_gt_info.split(":")[-2])

        elif "svision_pro" in vcf_path:
            child_support = int(child_gt_info.split(":")[-1]) if child_gt_info.split(":")[-1] != "NA" else 0
            mother_support = int(mother_gt_info.split(":")[-1]) if mother_gt_info.split(":")[-1] != "NA" else 0
            father_support = int(father_gt_info.split(":")[-1]) if father_gt_info.split(":")[-1] != "NA" else 0

        else:
            # # for merge strategy, we already set the min support to 10 when calling each separate sample
            child_support = min_support + 1
            mother_support = min_support + 1
            father_support = min_support + 1

        if child_support < min_support:
            child_gt = "0/0"

            if exist_in_child:
                continue

        if mother_support < min_support:
            mother_gt = "0/0"
            # continue

        if father_support < min_support:
            father_gt = "0/0"
            # continue

        # # STEP: convert ./. to 0/0 to unify gt representation
        if convert_gt:
            gt_convert = {"1/1": "1/1", "0/1": "0/1", "0/0": "0/0", "./.": "0/0"}

            child_gt = gt_convert[child_gt]
            mother_gt = gt_convert[mother_gt]
            father_gt = gt_convert[father_gt]

        cnt_dict["total"] += 1

        mendelian_res = determine_mendelian_error(child_gt, mother_gt, father_gt)

        # if mendelian_res == "bad":
        #     print(str(record).strip())

        cnt_dict[mendelian_res] += 1

    # print(cnt_dict["ok"] / (cnt_dict["total"]) , cnt_dict, vcf_path)
    print(round(cnt_dict["ok"] / cnt_dict["total"], 4))

if __name__ == '__main__':

    # print("------------------HiFi mendelian")
    # # #
    trios = ["hg002-3-4", "hg005-6-7", "na19240-38-39", "hg00514-2-3", "hg00733-1-2", "lcl5-6-7-8"]
    # # #
    print("------------SVision-pro")
    for trio in trios:
        # calculate_mendelian("/mnt/c/workspace/test_new/real_trio/release_{}_256.svision_pro_v1.6.s10.vcf".format(trio), min_support=1)
        calculate_mendelian("/mnt/c/workspace/test_new/real_trio/supp5/release_{}_256.svision_pro_v1.6.s5.vcf".format(trio), min_support=1)
    #
    print("------------sniffles2")
    for trio in trios:
       # calculate_mendelian("/mnt/c/workspace/test_new/real_trio/sniffles2.{}.s10.vcf".format(trio), min_support=1)
       calculate_mendelian("/mnt/c/workspace/test_new/real_trio/supp5/sniffles2.{}.s5.vcf".format(trio), min_support=1)

    print("------------svision_jasmine")
    for trio in trios:
        calculate_mendelian("/mnt/c/workspace/test_new/real_trio/supp5/svision_jasmine.{}.s5.vcf".format(trio), min_support=1)
    #
    print("------------svision_survivor")
    for trio in trios:
        calculate_mendelian("/mnt/c/workspace/test_new/real_trio/supp5/svision_survivor.{}.s5.vcf".format(trio), min_support=1)

    print("------------cutesv_jasmine")
    for trio in trios:
        calculate_mendelian("/mnt/c/workspace/test_new/real_trio/supp5/cutesv_jasmine.{}.s5.vcf".format(trio), min_support=1)

    print("------------cutesv_survivor")
    for trio in trios:
        calculate_mendelian("/mnt/c/workspace/test_new/real_trio/supp5/cutesv_survivor.{}.s5.vcf".format(trio), min_support=1)

    print("------------debreak_jasmine")
    for trio in trios:
        calculate_mendelian("/mnt/c/workspace/test_new/real_trio/supp5/debreak_jasmine.{}.s5.vcf".format(trio), min_support=1)

    print("------------debreak_survivor")
    for trio in trios:
        calculate_mendelian("/mnt/c/workspace/test_new/real_trio/supp5/debreak_survivor.{}.s5.vcf".format(trio), min_support=1)


    print("------------------ONT mendelian")

    trios = ["hg002-3-4", "hg005-6-7", "lcl5-6-7-8"]
    # print("------------SVision-pro")
    for trio in trios:
        calculate_mendelian("/mnt/c/workspace/test_new/real_trio/ont/release_{}_ONT_256.svision_pro_v1.6.s10.vcf".format(trio), min_support=10)

    print("------------sniffles2")
    for trio in trios:
        calculate_mendelian("/mnt/c/workspace/test_new/real_trio/ont/sniffles2.{}_ONT.s10.vcf".format(trio), min_support=10)
    # #
    print("------------svision_jasmine")
    for trio in trios:
        calculate_mendelian("/mnt/c/workspace/test_new/real_trio/ont/svision_jasmine.{}_ONT.s10.vcf".format(trio), min_support=10)
    # #
    print("------------svision_survivor")
    for trio in trios:
        calculate_mendelian("/mnt/c/workspace/test_new/real_trio/ont/svision_survivor.{}_ONT.s10.vcf".format(trio), min_support=10)
    #
    print("------------cutesv_jasmine")
    for trio in trios:
        calculate_mendelian("/mnt/c/workspace/test_new/real_trio/ont/cutesv_jasmine.{}_ONT.s10.vcf".format(trio), min_support=10)
    #
    print("------------cutesv_survivor")
    for trio in trios:
        calculate_mendelian("/mnt/c/workspace/test_new/real_trio/ont/cutesv_survivor.{}_ONT.s10.vcf".format(trio), min_support=10)

    print("------------debreak_jasmine")
    for trio in trios:
        calculate_mendelian("/mnt/c/workspace/test_new/real_trio/ont/debreak_jasmine.{}_ONT.s10.vcf".format(trio), min_support=10)

    print("------------debreak_survivor")
    for trio in trios:
        calculate_mendelian("/mnt/c/workspace/test_new/real_trio/ont/debreak_survivor.{}_ONT.s10.vcf".format(trio), min_support=10)


    #

    print("------------------quartet twin hifi")

    calculate_quartet_twin("/mnt/c/workspace/test_new/real_trio/release_lcl5-6-7-8_256.svision_pro_v1.6.s10.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/c/workspace/test_new/real_trio/supp5/release_lcl5-6-7-8_256.svision_pro_v1.6.s5.vcf", output_discordant=False)

    calculate_quartet_twin("/mnt/c/workspace/test_new/real_trio/sniffles2.lcl5-6-7-8.s10.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/c/workspace/test_new/real_trio/supp5/sniffles2.lcl5-6-7-8.s5.vcf", output_discordant=False)

    calculate_quartet_twin("/mnt/c/workspace/test_new/real_trio/svision_survivor.lcl5-6-7-8.s10.vcf", output_discordant=False)
    calculate_quartet_twin("/mnt/c/workspace/test_new/real_trio/svision_jasmine.lcl5-6-7-8.s10.vcf", output_discordant=False)
    calculate_quartet_twin("/mnt/c/workspace/test_new/real_trio/cutesv_survivor.lcl5-6-7-8.s10.vcf", output_discordant=False)
    calculate_quartet_twin("/mnt/c/workspace/test_new/real_trio/cutesv_jasmine.lcl5-6-7-8.s10.vcf", output_discordant=False)
    calculate_quartet_twin("/mnt/c/workspace/test_new/real_trio/debreak_survivor.lcl5-6-7-8.s10.vcf", output_discordant=False)
    calculate_quartet_twin("/mnt/c/workspace/test_new/real_trio/debreak_jasmine.lcl5-6-7-8.s10.vcf", output_discordant=False)

    #
    print("------------------quartet twin ont")
    calculate_quartet_twin("/mnt/c/workspace/test_new/real_trio/ont/release_lcl5-6-7-8_ONT_256.svision_pro_v1.6.s10.vcf", output_discordant=False)

    calculate_quartet_twin("/mnt/c/workspace/test_new/real_trio/ont/sniffles2.lcl5-6-7-8_ONT.s10.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/c/workspace/test_new/real_trio/ont/svision_survivor.lcl5-6-7-8_ONT.s10.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/c/workspace/test_new/real_trio/ont/svision_jasmine.lcl5-6-7-8_ONT.s10.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/c/workspace/test_new/real_trio/ont/cutesv_survivor.lcl5-6-7-8_ONT.s10.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/c/workspace/test_new/real_trio/ont/cutesv_jasmine.lcl5-6-7-8_ONT.s10.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/c/workspace/test_new/real_trio/ont/debreak_survivor.lcl5-6-7-8_ONT.s10.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/c/workspace/test_new/real_trio/ont/debreak_jasmine.lcl5-6-7-8_ONT.s10.vcf", output_discordant=False)
