import pysam

import os

def overlap_svision_pro_with_others(input_vcf, output_path, trio):

    other_call_path = {
        "sniffles2": "/mnt/d/workspace/svision-pro/real_trio/sniffles2.{}.s10.vcf".format(trio),
        "svision_jasmine": "/mnt/d/workspace/svision-pro/real_trio/supp10_supp2/svision_jasmine.{}.s10.s2.vcf".format(trio),
        "svision_survivor": "/mnt/d/workspace/svision-pro/real_trio/supp10_supp2/svision_survivor.{}.s10.s2.vcf".format(trio),
        "cutesv_jasmine": "/mnt/d/workspace/svision-pro/real_trio/supp10_supp2/cutesv_jasmine.{}.s10.s2.vcf".format(trio),
        "cutesv_survivor": "/mnt/d/workspace/svision-pro/real_trio/supp10_supp2/cutesv_survivor.{}.s10.s2.vcf".format(trio),
        "debreak_jasmine": "/mnt/d/workspace/svision-pro/real_trio/supp10_supp2/debreak_jasmine.{}.s10.s2.vcf".format(trio),
        "debreak_survivor": "/mnt/d/workspace/svision-pro/real_trio/supp10_supp2/debreak_survivor.{}.s10.s2.vcf".format(trio),

    }

    svision_pro_bed = os.path.join(output_path, "svision_pro.simplified.bed")

    with open(svision_pro_bed, "w") as fout, pysam.VariantFile(input_vcf) as fin:
        for record in fin:

            if "PASS" not in str(record):
                continue

            if "BND" in str(record):
                continue


            # # STEP: only keep autosome
            if record.contig not in ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", ]:
                continue

            if record.info["SVTYPE"] == "INS":
                fout.write("{}\t{}\t{}\n".format(record.contig, record.start, record.start + record.info["SVLEN"]))

            else:
                fout.write("{}\t{}\t{}\n".format(record.contig, record.start, record.start + 1))

    for other_caller in other_call_path.keys():

        other_vcf = other_call_path[other_caller]

        simplifed_bed = os.path.join(output_path, other_caller + ".simplified.bed")

        with open(simplifed_bed, "w") as fout, pysam.VariantFile(other_vcf) as fin:
            for record in fin:

                if "PASS" not in str(record):
                    continue

                if "BND" in str(record):
                    continue

                # # STEP: only keep autosome
                if record.contig not in ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", ]:
                    continue

                # # STEP: must exist in child sample
                if "SUPP_VEC=1" not in str(record):
                    continue

                if record.info["SVTYPE"] == "INS":
                    fout.write("{}\t{}\t{}\n".format(record.contig, record.start, record.start + record.info["SVLEN"]))

                else:
                    fout.write("{}\t{}\t{}\n".format(record.contig, record.start, record.start + 1))
        overlap_bed = os.path.join(output_path, "svision-pro.{}.simplified.paralogous.bed".format(other_caller))

        os.system("bedtools intersect -wo -a {} -b {} > {}".format(svision_pro_bed, simplifed_bed, overlap_bed))

        no_overlap_cnt = 0
        with open(overlap_bed) as bed_in, pysam.VariantFile(input_vcf) as vcf_in:

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

                if record_id not in high_conf_set:

                    print(str(record).strip())
                    no_overlap_cnt += 1

        print(other_caller, no_overlap_cnt)


def extract_paralogous_by_bedtools(input_path, input_vcf, output_path):

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
    paralogous_bed = '/mnt/d/workspace/svision-pro/real_trio_paralogous/2023NC.paralogous_list.txt'

    overlap_bed = os.path.join(output_path, input_vcf + ".simplified.paralogous.bed")

    os.system("bedtools intersect -wo -a {} -b {} > {}".format(simplifed_bed, paralogous_bed, overlap_bed))


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


if __name__ == '__main__':
    trios = ["hg002-3-4", "hg005-6-7", "na19240-38-39", "hg00514-2-3", "hg00733-1-2", "lcl5-6-7-8"]

    # for trio in trios:
        # extract_paralogous_by_bedtools("/mnt/d/workspace/svision-pro/real_trio/", "release_{}_256.svision_pro_v1.7.s10.vcf".format(trio), "/mnt/d/workspace/svision-pro/real_trio_paralogous/hifi/")

    # extract_paralogous_by_bedtools("/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/", "release_withbed.svision_pro_v1.7.s2.vcf", "/mnt/d/workspace/svision-pro/real_trio_paralogous/hifi/")


    trios = ["lcl5-6-7-8"]


    for trio in trios:
        overlap_svision_pro_with_others("/mnt/d/workspace/svision-pro/real_trio_paralogous/hifi/release_{}_256.svision_pro_v1.7.s10.vcf".format(trio), "/mnt/d/workspace/svision-pro/real_trio_paralogous/hifi/overlap_with_others", trio)