import os
import pysam


def filter_survivor_by_supp(vcf_path, min_supp=10):
    print(vcf_path)
    out_vcf_path = vcf_path.replace(".vcf", ".s{}.vcf".format(min_supp))

    vcf_file = pysam.VariantFile(vcf_path)

    with open(out_vcf_path, "w") as out_vcf_file:

        out_vcf_file.write(str(vcf_file.header))

        for record in vcf_file:

            if "survivor" in vcf_path:
                child_supp = int(str(record).strip().split("\t")[-3].split(":")[3].split(",")[1])
            else:
                if str(record).strip().split("\t")[-3].split(":")[-2] == "NA":
                    child_supp = 0
                else:
                    child_supp = int(str(record).strip().split("\t")[-3].split(":")[-2])


            if child_supp >= min_supp:
                out_vcf_file.write((str(record)))



if __name__ == '__main__':

    trios = ["hg002-3-4", "hg005-6-7", "na19240-38-39", "hg00514-2-3", "hg00733-1-2", "lcl5-6-7-8"]

    for trio in trios:
        filter_survivor_by_supp("/mnt/d/workspace/svision-pro/real_trio/supp5/cutesv_jasmine.{}.s5.vcf".format(trio))
    for trio in trios:
        filter_survivor_by_supp("/mnt/d/workspace/svision-pro/real_trio/supp5/cutesv_survivor.{}.s5.vcf".format(trio))

    # for trio in trios:
    #     filter_survivor_by_supp("/mnt/d/workspace/svision-pro/real_trio/supp5/svision_survivor.{}.s5.vcf".format(trio))
    #
    # for trio in trios:
    #     filter_survivor_by_supp("/mnt/d/workspace/svision-pro/real_trio/supp5/svision_jasmine.{}.s5.vcf".format(trio))


    for trio in trios:
        filter_survivor_by_supp("/mnt/d/workspace/svision-pro/real_trio/supp5/debreak_survivor.{}.s5.vcf".format(trio))

    for trio in trios:
        filter_survivor_by_supp("/mnt/d/workspace/svision-pro/real_trio/supp5/debreak_jasmine.{}.s5.vcf".format(trio))