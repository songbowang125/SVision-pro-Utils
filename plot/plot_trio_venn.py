import venn5
from matplotlib import pyplot as plt

import os.path

import pysam
from pybedtools import BedTool

def overlap_high_confident(input_path, sample, tool):

    output_path = os.path.join(input_path, "{}_{}".format(sample, tool))

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    # # STEP: convert vcf to bed
    for vcf_file in os.listdir(input_path):
        if "vcf" not in vcf_file:
            continue
        if sample not in vcf_file:
            continue
        if tool not in vcf_file and "sniffles2" not in vcf_file and "svision_pro" not in vcf_file:
            continue

        bed_file = open(os.path.join(output_path, vcf_file.replace(".vcf", ".bed")), "w")
        vcf_file = pysam.VariantFile(os.path.join(input_path, vcf_file))

        for record in vcf_file:
            if "PASS" not in str(record):
                continue

            if "BND" in str(record):
                continue

            bed_file.write("{}\t{}\t{}\n".format(record.contig, record.start + 1, record.start + 1 + abs(record.info["SVLEN"]), record.info["SVTYPE"]))

        vcf_file.close()
        bed_file.close()

    # # STEP: using pybedtools to find overlap
    pro_bed = os.path.join(output_path, "release_{}_256.svision_pro_v1.7.s10.bed".format(sample))
    sniffles2_bed = os.path.join(output_path, "sniffles2.{}.s10.bed".format(sample))
    cutesv_bed = os.path.join(output_path, "cutesv_{}.{}.s10.bed".format(tool, sample))
    svision_bed = os.path.join(output_path, "svision_{}.{}.s10.bed".format(tool, sample))
    debreak_bed = os.path.join(output_path, "debreak_{}.{}.s10.bed".format(tool, sample))

    pro_bed_res = BedTool(pro_bed)
    sniffles2_bed_res = BedTool(sniffles2_bed)
    cutesv_bed_res = BedTool(cutesv_bed)
    svision_bed_res = BedTool(svision_bed)
    debreak_bed_res = BedTool(debreak_bed)

    merged_bed = os.path.join(output_path, "merged.{}.{}.bed".format(tool, sample))

    os.system("cat {} {} {} {} {} > {}".format(pro_bed, sniffles2_bed, cutesv_bed, svision_bed, debreak_bed, merged_bed))

    merged_bed2_res = BedTool(merged_bed).sort().merge()

    merged_bed2 = os.path.join(output_path, "merged2.{}.{}.bed".format(tool, sample))
    merged_bed2_res.saveas(merged_bed2)

    # # add id
    merged_bed3 = os.path.join(output_path, "merged3.{}.{}.bed".format(tool, sample))

    id = 1
    with open(merged_bed3, "w") as fout, open(merged_bed2) as fin:
        for line in fin:
            fout.write(line.strip() + "\t" + str(id) + "\n")

            id += 1

    merged_bed3_res = BedTool(merged_bed3)

    # # do overlap
    pro_include_ids = []
    for region in pro_bed_res.intersect(merged_bed3_res, wao=True):
        pro_include_ids.append(int(region[6]))

    sniffles2_include_ids = []
    for region in sniffles2_bed_res.intersect(merged_bed3_res, wao=True):
        sniffles2_include_ids.append(int(region[6]))

    cutesv_include_ids = []
    for region in cutesv_bed_res.intersect(merged_bed3_res, wao=True):
        cutesv_include_ids.append(int(region[6]))

    svision_include_ids = []
    for region in svision_bed_res.intersect(merged_bed3_res, wao=True):
        svision_include_ids.append(int(region[6]))

    debreak_include_ids = []
    for region in debreak_bed_res.intersect(merged_bed3_res, wao=True):
        debreak_include_ids.append(int(region[6]))

    print(len(pro_include_ids), len(sniffles2_include_ids), len(cutesv_include_ids), len(svision_include_ids), len(debreak_include_ids))

    return pro_include_ids, sniffles2_include_ids, cutesv_include_ids, svision_include_ids, debreak_include_ids


if __name__ == '__main__':


    # samples = ["hg002-3-4", "hg005-6-7", "hg00514-2-3", "hg00733-1-2", "na19240-38-39", "lcl5-6-7-8"]
    # tools = ['survivor', "jasmine"]

    samples = ["hg002-3-4"]
    tools = ["jasmine"]

    input_path = "/mnt/d/workspace/svision-pro/real_trio_highconf/supp10/hifi"

    for sample in samples:
        for tool in tools:
            pro_include_ids, sniffles2_include_ids, cutesv_include_ids, svision_include_ids, debreak_include_ids = overlap_high_confident(input_path, sample, tool)

            labels = venn5.get_labels([pro_include_ids, sniffles2_include_ids, cutesv_include_ids, svision_include_ids, debreak_include_ids], fill=['number', 'logic'])

            fig, ax = venn5.venn5(labels, names=['SVision-pro', 'Sniffles2', 'cutesv_{}'.format(tool), 'svision_{}'.format(tool), 'debreak_{}'.format(tool)])

            plt.savefig("test.png", dpi=1500, bbox_inches='tight')

    # labels = venn5.get_labels([range(10), range(5, 15), range(3, 8), range(8, 17), range(10, 20)], fill=['number', 'logic'])
    # fig, ax = venn5.venn5(labels, names=['list 1', 'list 2', 'list 3', 'list 4', 'list 5'])
    # fig.show()
    #
    # plt.savefig("test.png")