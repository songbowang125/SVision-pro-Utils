import os.path

import pysam
from pybedtools import BedTool

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

    print(cnt_dict["ok"] / (cnt_dict["total"]) , cnt_dict, vcf_path)
    # print(round(cnt_dict["ok"] / cnt_dict["total"], 4))


def extract_high_confident_by_bedtools(input_path, input_vcf, output_path):

    # # STEP 1: extract child records and remove bnds and remove unresoved contigs
    simplifed_bed = os.path.join(output_path, input_vcf + ".simplified.bed")

    with open(simplifed_bed, "w") as fout, pysam.VariantFile(os.path.join(input_path, input_vcf)) as fin:
        for record in fin:

            if "PASS" not in str(record):
                continue

            if "BND" in str(record):
                continue

            # # STEP: only keep autosome
            if record.contig not in ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", ]:
                continue

            # # STEP: must exist in child sample
            if "svision_pro" not in input_vcf:
                if "SUPP_VEC=1" not in str(record):
                    continue

            fout.write("{}\t{}\t{}\n".format(record.contig, record.start, record.start + 1))

    # # STEP 2: run bedtools intersect
    hg002_base_bed = '/mnt/d/workspace/svision-pro/real_trio_highconf/Tier1_v0.6.hg38.bed'

    overlap_bed = os.path.join(output_path, input_vcf + ".simplified.high_conf.bed")

    os.system("bedtools intersect -wo -a {} -b {} > {}".format(simplifed_bed, hg002_base_bed, overlap_bed))

    with open(overlap_bed) as bed_in, pysam.VariantFile(os.path.join(input_path, input_vcf)) as vcf_in, open(os.path.join(output_path, input_vcf), "w") as vcf_out:

        vcf_out.write(str(vcf_in.header))

        # # load from bed
        high_conf_set = []

        for line in bed_in:
            line_split = line.strip().split()

            high_conf_set.append("{}_{}".format(line_split[0], line_split[1]))

        # # load from vcf and output high conf record
        for record in vcf_in:

            if "PASS" not in str(record):
                continue

            if "BND" in str(record):
                continue

            record_id = "{}_{}".format(record.contig, record.start)

            if record_id in high_conf_set:
                vcf_out.write(str(record))

    os.remove(overlap_bed)
    os.remove(simplifed_bed)


def extract_high_confident_by_truvari(input_path, input_vcf, output_path):

    # # STEP1: extract child records and remove bnds and remove unresoved contigs
    simplifed_vcf = os.path.join(output_path, input_vcf + ".simplified.vcf")

    with open(simplifed_vcf, "w") as fout, pysam.VariantFile(os.path.join(input_path, input_vcf)) as fin:
        fout.write(str(fin.header))

        for record in fin:

            line = str(record)

            if "PASS" not in str(record):
                continue

            if "BND" in str(record):
                continue

            # # STEP: only keep autosome
            if record.contig not in ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", ]:
                continue

            # # STEP: must exist in child sample
            if "svision_pro" not in input_vcf:
                if "SUPP_VEC=1" not in str(record):
                    continue

            fout.write(str(record))

    # # STEP2: run truvari to get high confident calls
    hg002_base_vcf = '/mnt/d/workspace/svision-pro//real_ccs/ground_truth/HG002_SVs_Tier1_v0.6.vcf.gz'
    hg002_base_bed = '/mnt/d/workspace/svision-pro/real_trio_highconf/Tier1_v0.6.hg38.bed'
    ref = "/mnt/h/data/ref/grch38/GRCh38.d1.vd1.fa"

    sorted_vcf = os.path.join(output_path, input_vcf + ".srt.vcf")
    cmd_str = 'bcftools sort -o {0} {1}'.format(sorted_vcf, simplifed_vcf)
    os.system(cmd_str)

    gz_vcf = sorted_vcf + '.gz'

    cmd_str = 'bgzip -c {0} > {1}'.format(sorted_vcf, gz_vcf)
    os.system(cmd_str)

    cmd_str = 'tabix -p vcf {0}'.format(gz_vcf)
    os.system(cmd_str)

    truvari_output_path = os.path.join(output_path, input_vcf.replace(".vcf", ""))
    cmd_str = 'truvari bench -b {} -c {} --includebed {} --typeignore  -r 1000 -p 0 -f {} -o {}'.format(hg002_base_vcf, gz_vcf, hg002_base_bed, ref, truvari_output_path)
    os.system(cmd_str)

    os.system("mv {} {}".format(os.path.join(truvari_output_path, "fp.vcf"), os.path.join(output_path, input_vcf)))

    os.remove(simplifed_vcf)
    os.remove(sorted_vcf)
    os.remove(gz_vcf)
    os.remove(gz_vcf + ".tbi")


if __name__ == '__main__':

    print("------------------HiFi mendelian")

    trios = ["hg002-3-4", "hg005-6-7", "na19240-38-39", "hg00514-2-3", "hg00733-1-2", "lcl5-6-7-8"]
    #
    # trios = ["na19240-38-39", "hg00514-2-3", "hg00733-1-2",]
    #
    print("------------SVision-pro")
    for trio in trios:
    #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio/release_{}_256.svision_pro_v1.7.s10.vcf".format(trio), min_support=1)

        # extract_high_confident_by_bedtools("/mnt/d/workspace/svision-pro/real_trio/supp5/", "release_{}_256.svision_pro_v1.6.s5.vcf".format(trio), "/mnt/d/workspace/svision-pro/real_trio_highconf/hifi/")

        calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio/release_{}_256.svision_pro_v1.6.s10.vcf".format(trio), min_support=1)
    # # # #
    # print("------------sniffles2")
    # for trio in trios:
    # #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio/sniffles2.{}.s10.vcf".format(trio), min_support=1)
    #
    #     extract_high_confident_by_bedtools("/mnt/d/workspace/svision-pro/real_trio/supp5/", "sniffles2.{}.s5.vcf".format(trio), "/mnt/d/workspace/svision-pro/real_trio_highconf/hifi/")
    #     #
    #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio_highconf/hifi/sniffles2.{}.s5.vcf".format(trio), min_support=1)
    #
    # print("------------svision_jasmine")
    # for trio in trios:
    # #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio/supp10_supp2/svision_jasmine.{}.s10.s2.vcf".format(trio), min_support=1)
    # #
    #     extract_high_confident_by_bedtools("/mnt/d/workspace/svision-pro/real_trio/supp5/", "svision_jasmine.{}.s5.vcf".format(trio), "/mnt/d/workspace/svision-pro/real_trio_highconf/hifi/")
    #     #
    #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio_highconf/hifi/svision_jasmine.{}.s5.vcf".format(trio), min_support=1)
    #
    #
    # print("------------svision_survivor")
    # for trio in trios:
    # #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio/supp10_supp2/svision_survivor.{}.s10.s2.vcf".format(trio), min_support=1)
    #
    #     extract_high_confident_by_bedtools("/mnt/d/workspace/svision-pro/real_trio/supp5/", "svision_survivor.{}.s5.vcf".format(trio), "/mnt/d/workspace/svision-pro/real_trio_highconf/hifi/")
    #     #
    #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio_highconf/hifi/svision_survivor.{}.s5.vcf".format(trio), min_support=1)
    #
    #
    # print("------------cutesv_jasmine")
    # for trio in trios:
    # #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio/supp10_supp2/cutesv_jasmine.{}.s10.s2.vcf".format(trio), min_support=1)
    #
    #     extract_high_confident_by_bedtools("/mnt/d/workspace/svision-pro/real_trio/supp5/", "cutesv_jasmine.{}.s5.vcf".format(trio), "/mnt/d/workspace/svision-pro/real_trio_highconf/hifi/")
    #
    #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio_highconf/hifi/cutesv_jasmine.{}.s5.vcf".format(trio), min_support=1)
    #
    #
    # print("------------cutesv_survivor")
    # for trio in trios:
    # #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio/supp10_supp2/cutesv_survivor.{}.s10.s2.vcf".format(trio), min_support=1)
    #
    #     extract_high_confident_by_bedtools("/mnt/d/workspace/svision-pro/real_trio/supp5/", "cutesv_survivor.{}.s5.vcf".format(trio), "/mnt/d/workspace/svision-pro/real_trio_highconf/hifi/")
    #     #
    #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio_highconf/hifi/cutesv_survivor.{}.s5.vcf".format(trio), min_support=1)
    #
    #
    # print("------------debreak_jasmine")
    # for trio in trios:
    # #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio/supp10_supp2/debreak_jasmine.{}.s10.s2.vcf".format(trio), min_support=1)
    # #
    #     extract_high_confident_by_bedtools("/mnt/d/workspace/svision-pro/real_trio/supp5/", "debreak_jasmine.{}.s5.vcf".format(trio), "/mnt/d/workspace/svision-pro/real_trio_highconf/hifi/")
    #
    #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio_highconf/hifi/debreak_jasmine.{}.s5.vcf".format(trio), min_support=1)
    # #
    # #
    # print("------------debreak_survivor")
    # for trio in trios:
    # #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio/supp10_supp2/debreak_survivor.{}.s10.s2.vcf".format(trio), min_support=1)
    # #
    #     extract_high_confident_by_bedtools("/mnt/d/workspace/svision-pro/real_trio/supp5/", "debreak_survivor.{}.s5.vcf".format(trio), "/mnt/d/workspace/svision-pro/real_trio_highconf/hifi/")
    #     #
    #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio_highconf/hifi/debreak_survivor.{}.s5.vcf".format(trio), min_support=1)


    # #
    print("------------------ONT mendelian")
    # #
    trios = ["hg002-3-4", "hg005-6-7", "lcl5-6-7-8"]
    print("------------SVision-pro")
    for trio in trios:
    # #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio/ont/release_{}_ONT_256.svision_pro_v1.7.s10.vcf".format(trio), min_support=10)
    # #
    #     extract_high_confident_by_bedtools("/mnt/d/workspace/svision-pro/real_trio/ont", "release_{}_ONT_256.svision_pro_v1.7.s10.vcf".format(trio), "/mnt/d/workspace/svision-pro/real_trio_highconf/ont/")
    #     #
        calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio/ont/release_{}_ONT_256.svision_pro_v1.7.s10.vcf".format(trio), min_support=10)
    # #
    # print("------------sniffles2")
    # for trio in trios:
    # #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio/ont/sniffles2.{}_ONT.s10.vcf".format(trio), min_support=10)
    # #
    #     extract_high_confident_by_bedtools("/mnt/d/workspace/svision-pro/real_trio/ont", "sniffles2.{}_ONT.s10.vcf".format(trio), "/mnt/d/workspace/svision-pro/real_trio_highconf/ont/")
    #
    #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio_highconf/ont/sniffles2.{}_ONT.s10.vcf".format(trio), min_support=10)
    # #
    # # # # #
    # print("------------svision_jasmine")
    # for trio in trios:
    # #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio/ont/supp10_supp5/svision_jasmine.{}.s10.s5.vcf".format(trio), min_support=1)
    # # #
    #     extract_high_confident_by_bedtools("/mnt/d/workspace/svision-pro/real_trio/ont", "svision_jasmine.{}_ONT.s10.vcf".format(trio), "/mnt/d/workspace/svision-pro/real_trio_highconf/ont/")
    # #
    #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio_highconf/ont/svision_jasmine.{}_ONT.s10.vcf".format(trio), min_support=10)
    # #
    # #
    # # # #
    # print("------------svision_survivor")
    # for trio in trios:
    # #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio/ont/supp10_supp5/svision_survivor.{}.s10.s5.vcf".format(trio), min_support=1)
    #
    #     extract_high_confident_by_bedtools("/mnt/d/workspace/svision-pro/real_trio/ont", "svision_survivor.{}_ONT.s10.vcf".format(trio), "/mnt/d/workspace/svision-pro/real_trio_highconf/ont/")
    #     #
    #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio_highconf/ont/svision_survivor.{}_ONT.s10.vcf".format(trio), min_support=10)
    # #
    # #
    # # #
    # print("------------cutesv_jasmine")
    # for trio in trios:
    # #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio/ont/supp10_supp5/cutesv_jasmine.{}.s10.s5.vcf".format(trio), min_support=1)
    # #
    #     extract_high_confident_by_bedtools("/mnt/d/workspace/svision-pro/real_trio/ont/", "cutesv_jasmine.{}_ONT.s10.vcf".format(trio), "/mnt/d/workspace/svision-pro/real_trio_highconf/ont/")
    #     #
    #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio_highconf/ont/cutesv_jasmine.{}_ONT.s10.vcf".format(trio), min_support=10)
    # #
    # #
    # # #
    # print("------------cutesv_survivor")
    # for trio in trios:
    # #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio/ont/supp10_supp5/cutesv_survivor.{}.s10.s5.vcf".format(trio), min_support=1)
    #
    #     extract_high_confident_by_bedtools("/mnt/d/workspace/svision-pro/real_trio/ont/", "cutesv_survivor.{}_ONT.s10.vcf".format(trio), "/mnt/d/workspace/svision-pro/real_trio_highconf/ont/")
    #     #
    #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio_highconf/ont/cutesv_survivor.{}_ONT.s10.vcf".format(trio), min_support=10)
    #
    #
    # print("------------debreak_jasmine")
    # for trio in trios:
    # #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio/ont/supp10_supp5/debreak_jasmine.{}.s10.s5.vcf".format(trio), min_support=1)
    #
    #     extract_high_confident_by_bedtools("/mnt/d/workspace/svision-pro/real_trio/ont", "debreak_jasmine.{}_ONT.s10.vcf".format(trio), "/mnt/d/workspace/svision-pro/real_trio_highconf/ont/")
    #     #
    #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio_highconf/ont/debreak_jasmine.{}_ONT.s10.vcf".format(trio), min_support=10)
    #
    #
    # print("------------debreak_survivor")
    # for trio in trios:
    # #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio/ont/supp10_supp5/debreak_survivor.{}.s10.s5.vcf".format(trio), min_support=1)
    # #
    #     extract_high_confident_by_bedtools("/mnt/d/workspace/svision-pro/real_trio/ont", "debreak_survivor.{}_ONT.s10.vcf".format(trio), "/mnt/d/workspace/svision-pro/real_trio_highconf/ont/")
    #     #
    #     calculate_mendelian("/mnt/d/workspace/svision-pro/real_trio_highconf/ont/debreak_survivor.{}_ONT.s10.vcf".format(trio), min_support=10)

    # print("------------------quartet twin hifi")

    # calculate_quartet_twin("/mnt/d/workspace/svision-pro/real_trio//svision_survivor.lcl5-6-7-8.s10.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/d/workspace/svision-pro/real_trio//svision_jasmine.lcl5-6-7-8.s10.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/d/workspace/svision-pro/real_trio//cutesv_survivor.lcl5-6-7-8.s10.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/d/workspace/svision-pro/real_trio//cutesv_jasmine.lcl5-6-7-8.s10.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/d/workspace/svision-pro/real_trio//debreak_survivor.lcl5-6-7-8.s10.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/d/workspace/svision-pro/real_trio//debreak_jasmine.lcl5-6-7-8.s10.vcf", output_discordant=False)

    # calculate_quartet_twin("/mnt/d/workspace/svision-pro/real_trio/release_lcl5-6-7-8_256.svision_pro_v1.7.s10.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/d/workspace/svision-pro/real_trio/sniffles2.lcl5-6-7-8.s10.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/d/workspace/svision-pro/real_trio/supp10_supp2/svision_survivor.lcl5-6-7-8.s10.s2.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/d/workspace/svision-pro/real_trio/supp10_supp2/svision_jasmine.lcl5-6-7-8.s10.s2.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/d/workspace/svision-pro/real_trio/supp10_supp2/cutesv_survivor.lcl5-6-7-8.s10.s2.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/d/workspace/svision-pro/real_trio/supp10_supp2/cutesv_jasmine.lcl5-6-7-8.s10.s2.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/d/workspace/svision-pro/real_trio/supp10_supp2/debreak_survivor.lcl5-6-7-8.s10.s2.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/d/workspace/svision-pro/real_trio/supp10_supp2/debreak_jasmine.lcl5-6-7-8.s10.s2.vcf", output_discordant=False)

    # calculate_quartet_twin("/mnt/d/workspace/svision-pro/real_trio_highconf/hifi/release_lcl5-6-7-8_256.svision_pro_v1.7.s10.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/d/workspace/svision-pro/real_trio_highconf/hifi/sniffles2.lcl5-6-7-8.s10.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/d/workspace/svision-pro/real_trio_highconf/hifi/svision_survivor.lcl5-6-7-8.s10.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/d/workspace/svision-pro/real_trio_highconf/hifi/svision_jasmine.lcl5-6-7-8.s10.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/d/workspace/svision-pro/real_trio_highconf/hifi/cutesv_survivor.lcl5-6-7-8.s10.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/d/workspace/svision-pro/real_trio_highconf/hifi/cutesv_jasmine.lcl5-6-7-8.s10.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/d/workspace/svision-pro/real_trio_highconf/hifi/debreak_survivor.lcl5-6-7-8.s10.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/d/workspace/svision-pro/real_trio_highconf/hifi/debreak_jasmine.lcl5-6-7-8.s10.vcf", output_discordant=False)


    #
    # print("------------------quartet twin ont")
    # calculate_quartet_twin("/mnt/c/workspace/test_new/real_trio/ont/release_lcl5-6-7-8_ONT_256.svision_pro_v1.6.s10.vcf", output_discordant=False)
    #
    # calculate_quartet_twin("/mnt/c/workspace/test_new/real_trio/ont/sniffles2.lcl5-6-7-8_ONT.s10.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/c/workspace/test_new/real_trio/ont/svision_survivor.lcl5-6-7-8_ONT.s10.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/c/workspace/test_new/real_trio/ont/svision_jasmine.lcl5-6-7-8_ONT.s10.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/c/workspace/test_new/real_trio/ont/cutesv_survivor.lcl5-6-7-8_ONT.s10.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/c/workspace/test_new/real_trio/ont/cutesv_jasmine.lcl5-6-7-8_ONT.s10.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/c/workspace/test_new/real_trio/ont/debreak_survivor.lcl5-6-7-8_ONT.s10.vcf", output_discordant=False)
    # calculate_quartet_twin("/mnt/c/workspace/test_new/real_trio/ont/debreak_jasmine.lcl5-6-7-8_ONT.s10.vcf", output_discordant=False)
