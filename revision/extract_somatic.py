import os
import pysam


def extrack_somatic(vcf_path, output_folder, tool):

    min_supports = [2, 3, 5, 8, 10, 15, 20]
    # min_supports = [1, 2, 3, 5, 8, 10]
    # min_supports = [15, 20]
    # min_supports = [1]

    for min_support in min_supports:

        output_path = os.path.join(output_folder, vcf_path.split("/")[-1] + ".somatic.s{}.vcf".format(min_support))
        output_file = open(output_path, "w")

        vcf_file = pysam.VariantFile(vcf_path)

        output_file.write(str(vcf_file.header))

        cnt = 0
        for record in vcf_file:

            if "BND" in str(record):
                continue

            # if "PASS" not in str(record):
            #     continue

            if "POLY" in str(record):
                continue

            if tool == "svision-pro":
                if "NewCOMP" in str(record) and record.info["SUPPORT"] >= min_support:
                    cnt += 1

                    output_file.write(str(record))

            elif tool == "sniffles2":
                # if record.info["SUPPORT"] >= min_support:
                #     output_file.write(str(record))

                if "SUPP_VEC=10" in str(record) and record.info["SUPPORT"] >= min_support:
                    cnt += 1

                    output_file.write(str(record))

            elif tool == "sniffles2-non-germline":
                if record.info["SUPPORT"] >= min_support:
                    cnt += 1

                    output_file.write(str(record))
            elif tool == "nanomonsv":

                tumor_gt = str(record).strip().split("\t")[-2]
                tumor_sv_reads = int(tumor_gt.split(":")[-1])

                control_gt = str(record).strip().split("\t")[-1]

                control_sv_reads = int(control_gt.split(":")[-1])

                if control_sv_reads == 0 and tumor_sv_reads >= min_support:
                    cnt += 1

                    output_file.write(str(record))
        print(cnt)
        vcf_file.close()
        output_file.close()


def filter_by_support(vcf_path, output_folder, tool):

    min_supports = [3, 5, 8, 10, 15, 20]
    # min_supports = [15, 20]

    for min_support in min_supports:

        output_path = os.path.join(output_folder, vcf_path.split("/")[-1] + ".s{}.vcf".format(min_support))
        output_file = open(output_path, "w")

        vcf_file = pysam.VariantFile(vcf_path)

        output_file.write(str(vcf_file.header))
        for record in vcf_file:
            if tool == "svision-pro":

                if record.info["SUPPORT"] >= min_support:
                    output_file.write(str(record))

            elif tool == "sniffles2":
                if record.info["SUPPORT"] >= min_support:
                    output_file.write(str(record))

            elif tool == "nanomonsv":

                tumor_gt = str(record).strip().split("\t")[-2]
                tumor_sv_reads = int(tumor_gt.split(":")[-1])

                control_gt = str(record).strip().split("\t")[-1]

                control_sv_reads = int(control_gt.split(":")[-1])

                if tumor_sv_reads >= min_support:
                    output_file.write(str(record))


if __name__ == '__main__':
    # extrack_somatic("/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/release_withbed.svision_pro_v1.7.s1.vcf",
    #                 "/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/test",
    #                 tool="svision-pro")

    extrack_somatic("/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/sniffles2.merge.s2.vcf",
                    "/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/test",
                    tool="sniffles2")

    # extrack_somatic("/mnt/d/workspace/svision-pro/sim_somatic_ccs/complex/release_1024.svision_pro_v1.7.s1.vcf",
    #                 "/mnt/d/workspace/svision-pro/sim_somatic_ccs/complex",
    #                 tool="svision-pro")

    # extrack_somatic("/mnt/d/workspace/svision-pro/sim_somatic_ccs/simple/sniffles2.merge.s1.vcf",
    #                 "/mnt/d/workspace/svision-pro/sim_somatic_ccs/simple",
    #                 tool="sniffles2")

    # extrack_somatic("/mnt/d/workspace/svision-pro/sim_somatic_ont/simple/nanomonsv.s1.vcf",
    #                 "/mnt/d/workspace/svision-pro/sim_somatic_ont/simple",
    #                 tool="nanomonsv")




    # extrack_somatic("/mnt/d/workspace/svision-pro/sim_somatic_ont/simple_100/release_1024_m20.svision_pro_v1.7.s1.vcf",
    #                 "/mnt/d/workspace/svision-pro/sim_somatic_ont/simple_100",
    #                 tool="svision-pro")
    #
    # extrack_somatic("/mnt/d/workspace/svision-pro/sim_somatic_ont/simple_100/sniffles2.merge.s1.vcf",
    #                 "/mnt/d/workspace/svision-pro/sim_somatic_ont/simple_100",
    #                 tool="sniffles2")



    # extrack_somatic("/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/sniffles2.merge.s2.vcf",
    #                 "/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic",
    #                 tool="sniffles2")

    # filter_by_support("/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/release_withbed.svision_pro_v1.7.s2.vcf",
    #                 "/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi",
    #                 tool="svision-pro")

    # extrack_somatic("/mnt/d/workspace/svision-pro/paper_revision/real_nt_colo829/release_1024.svision_pro_v1.6.s2.vcf",
    #                 "/mnt/d/workspace/svision-pro/paper_revision/real_nt_colo829/somatic",
    #                 tool="svision-pro")

    # extrack_somatic("/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/release.svision_pro_v1.6.s2.vcf",
    #                 "/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic",
    #                 tool="svision-pro")

    # extrack_somatic("/mnt/d/workspace/svision-pro/real_nt_hcc1395/sniffles.non_germline.s2.vcf",
    #                 "/mnt/d/workspace/svision-pro/paper_revision/somatic",
    #                 tool="sniffles2")

    # extrack_somatic("/mnt/d/workspace/svision-pro/paper_revision/real_nt_colo829/sniffles2.non-germline.s1.m20.vcf",
    #                 "/mnt/d/workspace/svision-pro/paper_revision/real_nt_colo829/somatic",
    #                 tool="sniffles2")
    #

    # extrack_somatic("/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/sniffles2.non-germline.s2.vcf",
    #                 "/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic",
    #                 tool="sniffles2")

    # extrack_somatic("/mnt/d/workspace/svision-pro/real_nt_hcc1395/nanomonsv.s2.vcf",
    #                 "/mnt/d/workspace/svision-pro/paper_revision/somatic",
    #                 tool="nanomonsv")