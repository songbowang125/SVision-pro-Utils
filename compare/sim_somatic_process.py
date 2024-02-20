import pysam
import os


def add_af_to_somatic_benmark_vcf(benchmark_vcf, out_vcf, somatic_log):
    out_vcf = open(out_vcf, "w")

    benchmark_vcf = pysam.VariantFile(benchmark_vcf)

    out_vcf.write(str(benchmark_vcf.header))

    for record in benchmark_vcf:
        record_chrom = record.contig
        record_start = record.start + 1
        record_stop = record.stop
        record_svtype = record.info["SVTYPE"]

        search_info = "{}_{}".format(record_chrom, record_start)      # # for SSV
        # search_info = "{}_{}".format(record_chrom, record_stop)     # # for CSV

        found_flag = False
        # # search in log file
        for line in open(somatic_log):

            line_split = line.strip().split("\t")

            target_info = "{}_{}".format(line_split[1], line_split[2])    # # for SSV
            # target_info = "{}_{}".format(line_split[1], line_split[3])  # # for CSV

            if search_info == target_info:

                found_flag = True
                # somatic_af = line_split[3]  # f
                somatic_af = line_split[-1]

                break

        if not found_flag:
            print(search_info)
        else:
            record_split = str(record).split("\t")

            record_split[7] = "SOMATIC={};".format(somatic_af) + record_split[7]

            out_vcf.write("\t".join(record_split))

    out_vcf.close()


def evaluation_somatic_concordance_from_truvari_res_nanomonsv(truvari_res_folder):
    # af_set = ["all", 0.2, 0.15, 0.10, 0.08, 0.05, 0.04, 0.03, 0.02, 0.01]
    af_set = ["all",  0.10, 0.08, 0.05, 0.04, 0.03, 0.02, 0.01]

    total_events = {}
    notcalled_events = {}
    accurate_events = {}
    inaccurate_events = {}
    for af in af_set:
        af = str(af)
        notcalled_events[af] = 0
        accurate_events[af] = 0
        inaccurate_events[af] = 0
        total_events[af] = 0


    # # STEP: use FN file to count NotCalled events
    for record in pysam.VariantFile(os.path.join(truvari_res_folder, "fn.vcf")):

        record_af = record.info["SOMATIC"]

        notcalled_events[record_af] += 1
        notcalled_events["all"] += 1

        total_events[record_af] += 1
        total_events["all"] += 1

    # # STEP: use tpcall and tp base to count accurate and inaccurate events
    event_infos = []
    for record in pysam.VariantFile(os.path.join(truvari_res_folder, "tp-base.vcf")):
        record_str = "{}-{}".format(record.contig, record.start + 1)
        record_af = record.info["SOMATIC"]
        event_infos.append([record_str, record_af, "NA"])

    i = 0
    for record in pysam.VariantFile(os.path.join(truvari_res_folder, "tp-call.vcf")):
        called_af = int(str(record).split("\t")[-2].split(":")[-2]) /  int(str(record).split("\t")[-2].split(":")[-1])

        # if "NA" in str(record).split("\t")[-1].split(":")[-1]:
        #     normal_reads = 0
        # else:
        normal_reads = int(str(record).split("\t")[-1].split(":")[-1])

        if normal_reads == 0:
            event_infos[i][2] = called_af
        else:
            event_infos[i][2] = "germline"
        i += 1

    for record_str, af, called_af in event_infos:

        if "germline" == called_af:
            inaccurate_events[af] += 1
            inaccurate_events['all'] += 1
        else:
            accurate_events[af] += 1
            accurate_events['all'] += 1

        total_events[af] += 1
        total_events["all"] += 1

    for af in af_set:
        af = str(af)
        # print(total_events[af])
        # print(accurate_events[af])
        # print(inaccurate_events[af] + notcalled_events[af])
        print(round(accurate_events[af] / total_events[af], 3))


def evaluation_somatic_concordance_from_truvari_res_sniffles2(truvari_res_folder):


    # af_set = ["all", 0.01, 0.02, 0.03, 0.05, 0.08, 0.10, 0.15, 0.20, 0.25]
    # af_set = ["all", 0.2, 0.15, 0.10, 0.08, 0.05, 0.04, 0.03, 0.02, 0.01]
    af_set = ["all",  0.10, 0.08, 0.05, 0.04, 0.03, 0.02, 0.01]

    total_events = {}
    notcalled_events = {}
    accurate_events = {}
    inaccurate_events = {}

    for af in af_set:
        af = str(af)
        notcalled_events[af] = 0
        accurate_events[af] = 0
        inaccurate_events[af] = 0
        total_events[af] = 0

    # # STEP: use FN file to count NotCalled events
    for record in pysam.VariantFile(os.path.join(truvari_res_folder, "fn.vcf")):

        record_af = record.info["SOMATIC"]

        notcalled_events[record_af] += 1
        notcalled_events["all"] += 1

        total_events[record_af] += 1
        total_events["all"] += 1

    # # STEP: use tpcall and tp base to count accurate and inaccurate events
    event_infos = []
    for record in pysam.VariantFile(os.path.join(truvari_res_folder, "tp-base.vcf")):
        record_str = "{}-{}".format(record.contig, record.start + 1)
        record_af = record.info["SOMATIC"]
        event_infos.append([record_str, record_af, "NA"])

    i = 0
    for record in pysam.VariantFile(os.path.join(truvari_res_folder, "tp-call.vcf")):

        try:
            called_af = record.info["AF"]
        except:
            called_af = 0

        event_infos[i][2] = called_af

        i += 1

    for record_str, af, called_af in event_infos:

        accurate_events[af] += 1
        accurate_events['all'] += 1

        total_events[af] += 1
        total_events["all"] += 1

    for af in af_set:
        af = str(af)
        # print(total_events[af])
        # print(accurate_events[af])
        # print(inaccurate_events[af] + notcalled_events[af])
        print(round(accurate_events[af] / total_events[af], 3))


def evaluation_somatic_concordance_from_truvari_res_svision_pro(truvari_res_folder):


    # af_set = ["all", 0.01, 0.02, 0.03, 0.05, 0.08, 0.10, 0.15, 0.20, 0.25]
    # af_set = ["all", 0.2, 0.15, 0.10, 0.08, 0.05, 0.04, 0.03, 0.02, 0.01]
    af_set = ["all",  0.10, 0.08, 0.05, 0.04, 0.03, 0.02, 0.01]

    total_events = {}
    notcalled_events = {}
    accurate_events = {}
    inaccurate_events = {}

    for af in af_set:
        af = str(af)
        notcalled_events[af] = 0
        accurate_events[af] = 0
        inaccurate_events[af] = 0
        total_events[af] = 0

    # # STEP: use FN file to count NotCalled events
    for record in pysam.VariantFile(os.path.join(truvari_res_folder, "fn.vcf")):

        record_af = record.info["SOMATIC"]

        # if record_af == "0.01":
        #     print("not callded", str(record).strip())

        notcalled_events[record_af] += 1
        notcalled_events["all"] += 1

        total_events[record_af] += 1
        total_events["all"] += 1

    # # STEP: use tpcall and tp base to count accurate and inaccurate events
    event_infos = []
    for record in pysam.VariantFile(os.path.join(truvari_res_folder, "tp-base.vcf")):
        record_str = "{}-{}".format(record.contig, record.start + 1)
        record_af = record.info["SOMATIC"]
        event_infos.append([record_str, record_af, "NA", 0.0])

    i = 0
    for record in pysam.VariantFile(os.path.join(truvari_res_folder, "tp-call.vcf")):
        it_types = []
        for record_it in record.info["BKPSIT"]:
            it_types.append(record_it.split("_")[1])

        event_infos[i][2] = "-".join(it_types)
        event_infos[i][3] = record.info["VAF"]

        i += 1

    vaf_diff = []
    for record_str, af, called_it, called_vaf in event_infos:
        # if "9370132" in record_str:
        #     print(record_str, af, called_it, called_vaf)
        if "Germline" in called_it or "NewALLELE" in called_it:
        # if "NewCOMP" not in called_it:

        # if "Germline" in called_it:

            # if af == "0.04":
            #     print('inaccurate', record_str, called_it, called_vaf)

            # print(record_str[0: 100], af, called_it, called_vaf)
            inaccurate_events[af] += 1
            inaccurate_events['all'] += 1
        else:
            accurate_events[af] += 1
            accurate_events['all'] += 1

        total_events[af] += 1
        total_events["all"] += 1

        vaf_diff.append(float(af) - called_vaf if float(af) - called_vaf > 0.001 else 0)

    # for vaf in vaf_diff:
    #     print(vaf)
    for af in af_set:
        af = str(af)
        # print(total_events[af])
        # print(accurate_events[af])
        # print(inaccurate_events[af] + notcalled_events[af])
        print(round(accurate_events[af] / total_events[af], 3))
        # print(round(accurate_events[af] / (accurate_events[af] + inaccurate_events[af]), 3))


if __name__ == '__main__':
    # add_af_to_somatic_benmark_vcf("/mnt/d/workspace/svision-pro/sim_somatic_ccs/complex/benchmark_origin_clone0_300.csv.vcf", "/mnt/d/workspace/svision-pro/sim_somatic_ccs/complex_100/benchmark_origin_clone0_300.csv.addVAF.vcf", "/mnt/d/workspace/svision-pro/sim_somatic_ccs/complex_100/csv_list.txt")
    # add_af_to_somatic_benmark_vcf("/mnt/d/workspace/svision-pro/sim_somatic_ccs/simple/HG002_SVs_Tier1_v0.6_high_confident.vcf", "/mnt/d/workspace/svision-pro/sim_somatic_ccs/simple_100_30x/HG002_SVs_Tier1_v0.6_high_confident.addVAF.vcf", "/mnt/d/workspace/svision-pro/sim_somatic_ccs/simple_100_30x/csv_list.txt")

    print("svision-pro")
    evaluation_somatic_concordance_from_truvari_res_svision_pro("/mnt/d/workspace/svision-pro/sim_somatic_ccs/complex/release_1024.svision_pro_v1.7.s1.vcf.somatic.s1_bench_region")
    print("---")
    evaluation_somatic_concordance_from_truvari_res_svision_pro("/mnt/d/workspace/svision-pro/sim_somatic_ont/complex/release_1024.svision_pro_v1.7.s1.vcf.somatic.s1_bench_region")

    print("sniffles2")
    # evaluation_somatic_concordance_from_truvari_res_sniffles2("/mnt/d/workspace/svision-pro/sim_somatic_ont/simple/sniffles2.merge.s1.vcf.somatic.s1_bench_region")

    print("nanomonsv")
    # evaluation_somatic_concordance_from_truvari_res_nanomonsv("/mnt/d/workspace/svision-pro/sim_somatic_ont/simple/nanomonsv.s1.vcf.somatic.s1_bench_region")

