import os
import pysam
import multiprocessing
import collections

def convert_benchmark_vcf_into_vapor_bed(vcf_path):

    index = 1

    for record in pysam.VariantFile(vcf_path):

        sv_chrom = record.contig

        sv_start = record.start + 1

        sv_type = record.info["SVTYPE"]

        sv_length = record.info["SVLEN"]

        if sv_type == "INS":
            sv_type = "{}_{}".format(sv_type, sv_length)

            sv_end = sv_start
        else:
            sv_end = record.stop + 1


        print("{}\t{}\t{}\tSV_{}\t{}".format(sv_chrom, sv_start, sv_end, index, sv_type))


        index += 1


def vapor_process_one(sv_str, bed_path, sample, seq, tmp_out_path, vapor_out_path):

    tumor_ont = "/mnt/e/HCC1395.ONT.minimap2.srt.bam"
    tumor_hifi = "/mnt/h/data/hcc1395/HiFi_2023Pacbio/HCC1395.GRCh38.bam"
    tumor_clr = "/mnt/h/data/hcc1395/HCC1395.PB.minimap2.srt.bam"

    normal_ont = "/mnt/e/HCC1395BL.ONT.minimap2.srt.bam"
    normal_hifi = "/mnt/h/data/hcc1395/HiFi_2023Pacbio/HCC1395-BL.GRCh38.bam"
    normal_clr = "/mnt/h/data/hcc1395/HCC1395BL.PB.minimap2.srt.bam"

    os.system("echo {} > {}".format(sv_str, bed_path))

    input_path = eval("{}_{}".format(sample, seq))

    os.system("/home/ubuntu/anaconda3/envs/vapor-env/bin/vapor bed --sv-input {} --reference /mnt/f/ref/grch38/GRCh38.d1.vd1.fa --pacbio-input {} --output-path {} --output-file {}".format(bed_path, input_path, tmp_out_path, vapor_out_path))


def convert_svision_vcf_into_vapor_bed(vcf_path, output_path, tool):

    # if not os.path.exists(output_path):
    #     os.mkdir(output_path)

    process_pool = multiprocessing.Pool(20)

    index = 1

    unprocessed_types = {}

    for record in pysam.VariantFile(vcf_path):

        sv_chrom = record.contig

        sv_start = record.start + 1

        sv_type = record.info["SVTYPE"]

        sv_length = record.info["SVLEN"]

        support = record.info["SUPPORT"]

        if sv_type == "INS":
            sv_type = "{}_{}".format(sv_type, sv_length)

            sv_end = sv_start

        elif sv_type == "tDUP":
            sv_type = "DUP"

            sv_end = record.stop + 1

        elif sv_type == "DUP":
            sv_end = record.stop + 1

        elif sv_type == "DEL":
            sv_end = record.stop + 1

        elif sv_type == "INV":
            sv_end = record.stop + 1

        else:

            if sv_type not in unprocessed_types:
                unprocessed_types[sv_type] = 0

            unprocessed_types[sv_type] += 1

            continue

        sv_str = "{}\t{}\t{}\tSV_{}\t{}\t{}".format(sv_chrom, sv_start, sv_end, index, sv_type, support)
        # sv_str = "{}\t{}\t{}\tSV_{}\t{}".format(sv_chrom, sv_start, sv_end, index, sv_type)

        print(sv_str)

        # for sample in ["tumor", "normal"]:
        #
        #     for seq in ["hifi", "clr", "ont"]:
        #
        #         tmp_out_path = os.path.join(output_path, "vapor_{}_{}_{}".format(tool, sample, seq))
        #
        #         if not os.path.exists(tmp_out_path):
        #             os.mkdir(tmp_out_path)
        #
        #         bed_path = os.path.join(tmp_out_path, "sv{}_{}_{}.bed".format(index, sample, seq))
        #
        #         vapor_out_path = os.path.join(tmp_out_path, "sv{}_{}_{}.vapor.bed".format(index, sample, seq))
        #
        #         # vapor_process_one(sv_str, bed_path, sample, seq, tmp_out_path, vapor_out_path)
        #
        #         process_pool.apply_async(vapor_process_one, (sv_str, bed_path, sample, seq, tmp_out_path, vapor_out_path))

        index += 1

    process_pool.close()
    process_pool.join()

    print(unprocessed_types)

def parse_vapor_res_somatic_with_repeat_with_af(input_prefix, repeat_res, truvari_res, raw_res, supp_thresh=5, score_thresh=0.1):


    candidate_seqs = ["hifi", "ont", "clr"]

    # candidate_seqs = ["hifi", "ont"]

    # # STEP 1. load file
    vapor_res = {}

    for line in open(repeat_res):
        line_split = line.strip().split("\t")
        sv_id = ".".join(line_split[0: 3] + line_split[4: 5])

        sv_support = int(line_split[5])
        if sv_support < supp_thresh:
            continue

        if sv_id not in vapor_res:
            vapor_res[sv_id] = {"AF": None, "Repeat": None, "SUPPPORT": None, "tumor_hifi": None, "tumor_ont": None, "tumor_hifi": None, "normal_hifi": None, "normal_ont": None, "normal_hifi": None}

        repeat_ratio = float(line_split[-3])
        vapor_res[sv_id]["Repeat"] = repeat_ratio

        # # load af
        found_flag = False
        raw_res_file = pysam.VariantFile(raw_res)
        for record in raw_res_file:

            if record.contig == line_split[0] and int(line_split[1]) == record.start + 1:
                found_flag = True
                if "sniffles2" in input_prefix:
                    tumor_gt_field = str(record).split("\t")[-2].split(":")

                    dr = int(tumor_gt_field[2])
                    dv = int(tumor_gt_field[3])
                    vapor_res[sv_id]["AF"] = round(dv / (dr + dv), 4)

                elif "svision_pro" in input_prefix:
                    tumor_gt_field = str(record).split("\t")[-2].split(":")

                    dr = int(tumor_gt_field[1])
                    dv = int(tumor_gt_field[2])
                    vapor_res[sv_id]["AF"] = round(dv / (dr + dv), 4)

                break

        raw_res_file.close()

    for sample in ["tumor", "normal"]:
        for seq in candidate_seqs:
            #
            vapor_res_file = open(input_prefix + ".{}.{}.txt".format(sample, seq))
            #
            for line in vapor_res_file:
                if "SV_description" in line:
                    continue
                if "#CHR" in line:
                    continue
                line_split = line.strip().split("\t")

                # sv_id = ".".join(line_split[0: 5])
                sv_id = ".".join(line_split[0: 3] + line_split[4: 5])

                sv_support = int(line_split[5])
                if sv_support < supp_thresh:
                    continue

                if line_split[7] == "":
                    line_split[7] = "NA"
                if line_split[8] == "":
                    line_split[8] = "NA"

                sv_vapor_score = line_split[7]
                sv_vapor_vaf = line_split[8]

                # vapor_res[sv_id]["SUPPPORT"] = sv_support
                vapor_res[sv_id]["{}_{}".format(sample, seq)] = [sv_vapor_score, sv_vapor_vaf]

    # # # STEP: remove already reported events
    # for line in open(truvari_res):
    #     line_split = line.strip().split("\t")
    #     sv_id = ".".join(line_split[0: 3] + line_split[4: 5])
    #
    #     if sv_id in vapor_res:
    #         vapor_res.pop(sv_id)


    for af_thresh in [0.05, 0.1, 0.2, 0.5, 1, 2]:
        print(af_thresh)
        # # STEP: stat

        total_num = 0

        inconclusive = []
        supp_3 = []
        supp_2 = []
        supp_1 = []

        for sv_id in vapor_res:
            if vapor_res[sv_id]["AF"] > af_thresh:
                continue

            total_num += 1
            sv_id_somatic_res = {}

            # # STEP: determine repeat
            if vapor_res[sv_id]["Repeat"] >= 50:
                for seq in candidate_seqs:
                    sv_id_somatic_res[seq] = "Repeat"
            else:
                for seq in candidate_seqs:
                    tumor_score = vapor_res[sv_id]["tumor_{}".format(seq)][0]
                    tumor_vaf = vapor_res[sv_id]["tumor_{}".format(seq)][1]

                    normal_score = vapor_res[sv_id]["normal_{}".format(seq)][0]
                    normal_vaf = vapor_res[sv_id]["normal_{}".format(seq)][1]

                    # # STEP: determine NA output
                    if tumor_score == "NA" or normal_score == "NA":
                        sv_id_somatic_res[seq] = "NA"

                    else:
                        if float(tumor_score) >= score_thresh:
                            if float(normal_score) >= score_thresh:
                                sv_id_somatic_res[seq] = "Fail"

                            else:
                                sv_id_somatic_res[seq] = "Pass"
                        else:
                            sv_id_somatic_res[seq] = "Fail"

            # # STEP: re organize res
            sv_id_somatic_res_stats = {"Pass": 0, "Fail": 0, "NA": 0, "Repeat": 0}

            for seq in candidate_seqs:
                sv_id_somatic_res_stats[sv_id_somatic_res[seq]] += 1

            if sv_id_somatic_res_stats["Repeat"] >= 3:
                inconclusive.append(sv_id)

            elif sv_id_somatic_res_stats["NA"] >= 1:
                inconclusive.append(sv_id)
            else:
                if sv_id_somatic_res_stats["Pass"] == 3:
                    supp_3.append(sv_id)

                if sv_id_somatic_res_stats["Pass"] == 2:
                    supp_2.append(sv_id)

                if sv_id_somatic_res_stats["Pass"] == 1:
                    supp_1.append(sv_id)



        inconclusive_num = len(inconclusive)
        inconclusive_ratio = round(inconclusive_num / total_num, 4)

        total_num_no_inconclusive = total_num - inconclusive_num

        supp_1_num = len(supp_1)
        supp_1_ratio = round(supp_1_num / total_num_no_inconclusive, 4)

        supp_2_num = len(supp_2)
        supp_2_ratio = round(supp_2_num / total_num_no_inconclusive, 4)

        supp_3_num = len(supp_3)
        supp_3_ratio = round(supp_3_num / total_num_no_inconclusive, 4)

        supp_1_2_3_num =  supp_1_num + supp_2_num + supp_3_num
        supp_1_2_3_ratio = round(supp_1_2_3_num / total_num_no_inconclusive, 4)

        fail_num = total_num_no_inconclusive - supp_1_2_3_num
        fail_ratio = round(fail_num / total_num_no_inconclusive, 4)

        print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(total_num, inconclusive_num, supp_1_2_3_ratio, supp_1_ratio, supp_2_ratio, supp_3_ratio, fail_ratio))



def parse_vapor_res_somatic_with_repeat(input_prefix, repeat_res, truvari_res, supp_thresh=5, score_thresh=0.1):

    candidate_seqs = ["hifi", "ont", "clr"]

    # candidate_seqs = ["hifi", "ont"]

    # # STEP 1. load file
    vapor_res = {}

    for line in open(repeat_res):
        line_split = line.strip().split("\t")
        sv_id = ".".join(line_split[0: 3] + line_split[4: 5])

        sv_support = int(line_split[5])
        if sv_support < supp_thresh:
            continue

        if sv_id not in vapor_res:
            vapor_res[sv_id] = {"Repeat": None, "SUPPPORT": None, "tumor_hifi": None, "tumor_ont": None, "tumor_hifi": None, "normal_hifi": None, "normal_ont": None, "normal_hifi": None}

        repeat_ratio = float(line_split[-3])
        vapor_res[sv_id]["Repeat"] = repeat_ratio

    for sample in ["tumor", "normal"]:
        for seq in candidate_seqs:
            #
            vapor_res_file = open(input_prefix + ".{}.{}.txt".format(sample, seq))
            # vapor_res_file = open(input_prefix + ".{}.{}.bed".format(sample, seq))
            #
            for line in vapor_res_file:
                if "SV_description" in line:
                    continue
                if "#CHR" in line:
                    continue
                line_split = line.strip().split("\t")

                # sv_id = ".".join(line_split[0: 5])
                sv_id = ".".join(line_split[0: 3] + line_split[4: 5])

                sv_support = int(line_split[5])
                if sv_support < supp_thresh:
                    continue

                if line_split[7] == "":
                    line_split[7] = "NA"
                if line_split[8] == "":
                    line_split[8] = "NA"

                sv_vapor_score = line_split[7]
                sv_vapor_vaf = line_split[8]

                # vapor_res[sv_id]["SUPPPORT"] = sv_support
                vapor_res[sv_id]["{}_{}".format(sample, seq)] = [sv_vapor_score, sv_vapor_vaf]

    # # # STEP: remove already reported events
    # for line in open(truvari_res):
    #     line_split = line.strip().split("\t")
    #     sv_id = ".".join(line_split[0: 3] + line_split[4: 5])
    #
    #     if sv_id in vapor_res:
    #         vapor_res.pop(sv_id)

    # # STEP: stat
    inconclusive = []
    supp_3 = []
    supp_2 = []
    supp_1 = []

    for sv_id in vapor_res:
        sv_id_somatic_res = {}

        # # STEP: determine repeat
        if vapor_res[sv_id]["Repeat"] >= 50:
            for seq in candidate_seqs:
                sv_id_somatic_res[seq] = "Repeat"
        else:
            for seq in candidate_seqs:
                tumor_score = vapor_res[sv_id]["tumor_{}".format(seq)][0]
                tumor_vaf = vapor_res[sv_id]["tumor_{}".format(seq)][1]

                normal_score = vapor_res[sv_id]["normal_{}".format(seq)][0]
                normal_vaf = vapor_res[sv_id]["normal_{}".format(seq)][1]

                # # STEP: determine NA output
                if tumor_score == "NA" or normal_score == "NA":
                    sv_id_somatic_res[seq] = "NA"

                else:
                    if float(tumor_score) >= score_thresh:
                        if float(normal_score) >= score_thresh:
                            sv_id_somatic_res[seq] = "Fail"

                        else:
                            sv_id_somatic_res[seq] = "Pass"
                    else:
                        sv_id_somatic_res[seq] = "Fail"

        # # STEP: re organize res
        sv_id_somatic_res_stats = {"Pass": 0, "Fail": 0, "NA": 0, "Repeat": 0}

        for seq in candidate_seqs:
            sv_id_somatic_res_stats[sv_id_somatic_res[seq]] += 1

        if sv_id_somatic_res_stats["Repeat"] >= 3:
            inconclusive.append(sv_id)

        elif sv_id_somatic_res_stats["NA"] >= 1:
            inconclusive.append(sv_id)
        else:
            if sv_id_somatic_res_stats["Pass"] == 3:
                supp_3.append(sv_id)

            if sv_id_somatic_res_stats["Pass"] == 2:
                supp_2.append(sv_id)

            if sv_id_somatic_res_stats["Pass"] == 1:
                supp_1.append(sv_id)


    total_num = len(vapor_res)

    inconclusive_num = len(inconclusive)
    inconclusive_ratio = round(inconclusive_num / total_num, 4)

    total_num_no_inconclusive = total_num - inconclusive_num

    supp_1_num = len(supp_1)
    supp_1_ratio = round(supp_1_num / total_num_no_inconclusive, 4)

    supp_2_num = len(supp_2)
    supp_2_ratio = round(supp_2_num / total_num_no_inconclusive, 4)

    supp_3_num = len(supp_3)
    supp_3_ratio = round(supp_3_num / total_num_no_inconclusive, 4)

    supp_1_2_3_num =  supp_1_num + supp_2_num + supp_3_num
    supp_1_2_3_ratio = round(supp_1_2_3_num / total_num_no_inconclusive, 4)

    fail_num = total_num_no_inconclusive - supp_1_2_3_num
    fail_ratio = round(fail_num / total_num_no_inconclusive, 4)

    print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(total_num, inconclusive_num, supp_1_2_3_ratio, supp_1_ratio, supp_2_ratio, supp_3_ratio, fail_ratio))


def parse_vapor_res_somatic(input_prefix, supp_thresh=5, score_thresh=0.1):

    candidate_seqs = ["hifi", "ont", "clr"]

    # candidate_seqs = ["hifi", "ont"]

    # # STEP 1. load file
    vapor_res = {}

    for sample in ["tumor", "normal"]:
        for seq in candidate_seqs:
            #
            # vapor_res_file = open(input_prefix + ".{}.{}.txt".format(sample, seq))
            vapor_res_file = open(input_prefix + ".{}.{}.repeat_annot.RO50_large.bed".format(sample, seq))
            #
            for line in vapor_res_file:
                if "SV_description" in line:
                    continue
                if "#CHR" in line:
                    continue
                line_split = line.strip().split("\t")

                sv_id = ".".join(line_split[0: 5])

                # if "DEL" not in line:
                #     continue

                sv_support = int(line_split[5])
                if sv_support < supp_thresh:
                    continue

                if line_split[7] == "":
                    line_split[7] = "NA"
                if line_split[8] == "":
                    line_split[8] = "NA"

                sv_vapor_score = line_split[7]
                sv_vapor_vaf = line_split[8]
                # sv_vapor_score = float(line_split[7]) if line_split[7] != 'NA' else -1000
                # sv_vapor_vaf = float(line_split[8]) if line_split[8] != 'NA' else -1000
                #
                # if "NA" == sv_vapor_score:
                #     continue

                if sv_id not in vapor_res:
                    vapor_res[sv_id] = {"SUPPPORT": None, "tumor_hifi": None, "tumor_ont": None, "tumor_hifi": None, "normal_hifi": None, "normal_ont": None, "normal_hifi": None}

                # vapor_res[sv_id]["SUPPPORT"] = sv_support
                vapor_res[sv_id]["{}_{}".format(sample, seq)] = [sv_vapor_score, sv_vapor_vaf]

    # # STEP 2.
    seq_res = {}
    seq_res_somatic = {}

    three_seq_final = {"pass": 0, "fail": 0, "na": 0}
    three_seq_final_somatic = {"pass": 0, "fail": 0, "na": 0}

    poped_sv_id = []
    for sv_id in vapor_res:
        for seq in candidate_seqs:
            tumor_score = vapor_res[sv_id]["tumor_{}".format(seq)][0]
            normal_score = vapor_res[sv_id]["normal_{}".format(seq)][0]

            if tumor_score == "NA" or normal_score == "NA":

                poped_sv_id.append(sv_id)
                break

    for sv_id in poped_sv_id:
        three_seq_final['na'] += 1
        three_seq_final_somatic['na'] += 1

        vapor_res.pop(sv_id)

    for sv_id in vapor_res:
        # seq_flags = {"pass": 0, "fail": 0}
        sv_id_res = {"pass": 0, "fail": 0}
        sv_id_res_somatic = {"pass": 0, "fail": 0}

        for seq in candidate_seqs:

            if seq not in seq_res:
                seq_res[seq] = {"pass": 0, "fail": 0, "NA": 0}
                seq_res_somatic[seq] = {"pass": 0, "fail": 0, "NA": 0}

            tumor_score = vapor_res[sv_id]["tumor_{}".format(seq)][0]
            tumor_vaf = vapor_res[sv_id]["tumor_{}".format(seq)][1]

            normal_score = vapor_res[sv_id]["normal_{}".format(seq)][0]
            normal_vaf = vapor_res[sv_id]["normal_{}".format(seq)][1]

            if float(tumor_score) >= score_thresh:
                seq_res[seq]["pass"] += 1
                sv_id_res["pass"] += 1

                if float(normal_score) >= score_thresh:
                    seq_res_somatic[seq]["fail"] += 1
                    sv_id_res_somatic["fail"] += 1

                else:
                    seq_res_somatic[seq]["pass"] += 1
                    sv_id_res_somatic["pass"] += 1

            else:
                seq_res[seq]["fail"] += 1
                sv_id_res["fail"] += 1

        if sv_id_res["pass"] >= 1:
            three_seq_final["pass"] += 1
        else:
            three_seq_final["fail"] += 1

        if sv_id_res_somatic["pass"] >= 1:
            three_seq_final_somatic["pass"] += 1
        else:
            three_seq_final_somatic["fail"] += 1

    for seq in seq_res:

        seq_res[seq]["ratio"] = round(seq_res[seq]["pass"] / (seq_res[seq]["pass"] + seq_res[seq]["fail"]), 3)
        seq_res_somatic[seq]["ratio"] = round(seq_res_somatic[seq]["pass"] / (seq_res_somatic[seq]["pass"] + seq_res_somatic[seq]["fail"]), 3)

    # print(seq_res)
    # print(seq_res_somatic)
    #
    # print(three_seq_final, three_seq_final["pass"] / (three_seq_final["pass"] + three_seq_final["fail"]))
    # print(three_seq_final_somatic, three_seq_final_somatic["pass"] / (three_seq_final_somatic["pass"] + three_seq_final_somatic["fail"]))
    print(three_seq_final_somatic["pass"] + three_seq_final_somatic["fail"],  round(three_seq_final_somatic["pass"] / (three_seq_final_somatic["pass"] + three_seq_final_somatic["fail"]), 4))

    # print(round(three_seq_final_somatic["pass"] / (three_seq_final_somatic["pass"] + three_seq_final_somatic["fail"]), 3))

    # print(len(vapor_res))


def parse_vapor_res_somatic2(input_prefix, supp_thresh=5, score_thresh=0.1):
    # # STEP 1. load file
    vapor_res = {}

    for sample in ["tumor", "normal"]:
        for seq in ["hifi", "ont", "clr"]:

            # vapor_res_file = open(input_prefix + ".{}.{}.txt".format(sample, seq))
            vapor_res_file = open(input_prefix + ".{}.{}.repeat_annot.RO50.bed".format(sample, seq))

            for line in vapor_res_file:

                if "SV_description" in line:
                    continue
                if "#CHR" in line:
                    continue

                line_split = line.strip().split("\t")

                sv_id = ".".join(line_split[0: 5])

                sv_support = int(line_split[5])
                sv_length = int(line_split[6])

                if sv_support < supp_thresh:
                    continue

                if line_split[7] == "":
                    line_split[7] = "NA"
                if line_split[8] == "":
                    line_split[8] = "NA"
                sv_vapor_score = float(line_split[7]) if line_split[7] != 'NA' else -1000
                sv_vapor_vaf = float(line_split[8]) if line_split[8] != 'NA' else -1000

                # sv_vapor_score = line_split[6]
                # sv_vapor_vaf = line_split[7]

                # if "NA" == sv_vapor_score:
                #     continue

                if sv_id not in vapor_res:
                    vapor_res[sv_id] = {"SUPPPORT": None, "tumor_hifi": None, "tumor_ont": None, "tumor_hifi": None, "normal_hifi": None, "normal_ont": None, "normal_hifi": None}

                vapor_res[sv_id]["SUPPPORT"] = sv_support
                vapor_res[sv_id]["{}_{}".format(sample, seq)] = [sv_vapor_score, sv_vapor_vaf]

    somatic_tot_num = len(vapor_res.keys())
    somatic_pass_num = 0
    for sv_id in vapor_res:

        tumor_scores = [vapor_res[sv_id]["tumor_hifi"][0], vapor_res[sv_id]["tumor_ont"][0], vapor_res[sv_id]["tumor_clr"][0]]

        if tumor_scores.count(-1000) >= 1:
            somatic_tot_num -= 1
            continue

        normal_scores = [vapor_res[sv_id]["normal_hifi"][0], vapor_res[sv_id]["normal_ont"][0], vapor_res[sv_id]["normal_clr"][0]]

        if normal_scores.count(-1000) >= 1:
            somatic_tot_num -= 1
            continue

        tumor_flags = ["pass" if score >= score_thresh else "fail" for score in tumor_scores]
        normal_flags = ["pass" if score >= score_thresh else "fail" for score in normal_scores]

        tumor_flags = collections.Counter(tumor_flags)
        normal_flags = collections.Counter(normal_flags)

        if "pass" not in tumor_flags:
            tumor_flags["pass"] = 0
        if "pass" not in normal_flags:
            normal_flags["pass"] = 0

        # print(sv_id, tumor_flags, normal_flags)

        # # at least support by two techs
        # # that is for tumor, at least two techs; for normal, less than two techs

        if tumor_flags['pass'] >= 3 and normal_flags["pass"] < 3:
            somatic_pass_num += 1

        else:
            pass
            # print(sv_id)
    # print(somatic_pass_num, somatic_tot_num, round(somatic_pass_num / somatic_tot_num, 3))
    print(round(somatic_pass_num / somatic_tot_num, 3))


def parse_vapor_res_fn(input_prefix, score_thresh=0.1):

    # # STEP 1. load file
    vapor_res = {}

    for sample in ["tumor", "normal"]:
        for seq in ["hifi", "ont", "clr"]:

            vapor_res_file = open(input_prefix + ".{}.{}.txt".format(sample, seq))
            # vapor_res_file = open(input_prefix + ".{}.{}.repeat_annot.RO50_large.bed".format(sample, seq))

            for line in vapor_res_file:
                if "SV_description" in line:
                    continue
                if "#CHR" in line:
                    continue
                line_split = line.strip().split("\t")

                sv_id = ".".join(line_split[0: 5])

                if line_split[6] == "":
                    line_split[6] = "NA"
                if line_split[7] == "":
                    line_split[7] = "NA"

                sv_vapor_score = line_split[6]
                sv_vapor_vaf = line_split[7]

                if sv_id not in vapor_res:
                    vapor_res[sv_id] = {"SUPPPORT": None, "tumor_hifi": None, "tumor_ont": None, "tumor_hifi": None, "normal_hifi": None, "normal_ont": None, "normal_hifi": None}

                # vapor_res[sv_id]["SUPPPORT"] = sv_support
                vapor_res[sv_id]["{}_{}".format(sample, seq)] = [sv_vapor_score, sv_vapor_vaf]

    # # STEP 2.
    seq_res = {}
    seq_res_somatic = {}

    three_seq_final = {"pass": 0, "fail": 0}
    three_seq_final_somatic = {"pass": 0, "fail": 0}

    poped_sv_id = []
    for sv_id in vapor_res:
        for seq in ["hifi", "ont", "clr"]:
            tumor_score = vapor_res[sv_id]["tumor_{}".format(seq)][0]
            normal_score = vapor_res[sv_id]["normal_{}".format(seq)][0]

            if tumor_score == "NA" or normal_score == "NA":

                poped_sv_id.append(sv_id)
                break

    for sv_id in poped_sv_id:
        vapor_res.pop(sv_id)

    for sv_id in vapor_res:
        # seq_flags = {"pass": 0, "fail": 0}
        sv_id_res = {"pass": 0, "fail": 0}
        sv_id_res_somatic = {"pass": 0, "fail": 0}

        for seq in ["hifi", "ont", "clr"]:

            if seq not in seq_res:
                seq_res[seq] = {"pass": 0, "fail": 0, "NA": 0}
                seq_res_somatic[seq] = {"pass": 0, "fail": 0, "NA": 0}

            tumor_score = vapor_res[sv_id]["tumor_{}".format(seq)][0]
            tumor_vaf = vapor_res[sv_id]["tumor_{}".format(seq)][1]

            normal_score = vapor_res[sv_id]["normal_{}".format(seq)][0]
            normal_vaf = vapor_res[sv_id]["normal_{}".format(seq)][1]

            if float(tumor_score) >= score_thresh:
                seq_res[seq]["pass"] += 1
                sv_id_res["pass"] += 1

                if float(normal_score) >= score_thresh:
                    seq_res_somatic[seq]["fail"] += 1
                    sv_id_res_somatic["fail"] += 1

                else:
                    seq_res_somatic[seq]["pass"] += 1
                    sv_id_res_somatic["pass"] += 1

            else:
                seq_res[seq]["fail"] += 1
                sv_id_res["fail"] += 1

        if sv_id_res["pass"] >= 2:
            three_seq_final["pass"] += 1
        else:
            three_seq_final["fail"] += 1

        if sv_id_res_somatic["pass"] >= 2:
            three_seq_final_somatic["pass"] += 1
            print(sv_id)
        else:
            three_seq_final_somatic["fail"] += 1

    for seq in seq_res:

        seq_res[seq]["ratio"] = round(seq_res[seq]["pass"] / (seq_res[seq]["pass"] + seq_res[seq]["fail"]), 3)
        seq_res_somatic[seq]["ratio"] = round(seq_res_somatic[seq]["pass"] / (seq_res_somatic[seq]["pass"] + seq_res_somatic[seq]["fail"]), 3)

    # print(seq_res)
    # print(seq_res_somatic)
    #
    # print(three_seq_final, three_seq_final["pass"] / (three_seq_final["pass"] + three_seq_final["fail"]))
    print(three_seq_final_somatic, three_seq_final_somatic["pass"] / (three_seq_final_somatic["pass"] + three_seq_final_somatic["fail"]))
    # print(round(three_seq_final_somatic["pass"] / (three_seq_final_somatic["pass"] + three_seq_final_somatic["fail"]), 3))

    # print(len(vapor_res))


def add_length(bed_path):

    with open(bed_path) as fin, open(bed_path.replace(".bed", ".length.bed"), "w") as fout:

        for line in fin:

            line_split = line.strip().split()

            if "INS" in line:
                line_split.append(line_split[4].split("_")[1])
                line_split[2] = str(int(line_split[1]) + int(line_split[4].split("_")[1]))
            else:

                line_split.append(str(int(line_split[2]) - int(line_split[1])))


            fout.write("\t".join(line_split) + "\n")


if __name__ == '__main__':

    # add_length("/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic/sniffles2.merge.s2.vcf.somatic.s2.vapor.support.bed")


    # convert_benchmark_vcf_into_vapor_bed("/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/release_withbed.svision_pro_v1.7.s2_bench_region/fn.vcf")

    #
    # convert_svision_vcf_into_vapor_bed("/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic/sniffles2.merge.s2.vcf.somatic.s2.vcf", "/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic", "sniffles2")
    # convert_svision_vcf_into_vapor_bed("/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic/sniffles2.merge.s2.vcf.somatic.s2_bench_region/tp-call.vcf", "/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic", "sniffles2")
    # exit()
    # convert_svision_vcf_into_vapor_bed("/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic/release.svision_pro_v1.7.s2.vcf.somatic.s2.vcf", "/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic", "sniffles2")

    # convert_svision_vcf_into_vapor_bed("/mnt/d/workspace/svision-pro/paper_revision/denovo_vapor/denovo.sv.hg00733-2-3.vcf", "/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic", "sniffles2")

    # for supp in [2, 3, 5, 8, 10, 15, 20]:
    #
    #     print("----------", supp)
    #     parse_vapor_res_somatic("/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic/release.svision_pro_v1.7.s2.vcf.somatic.s2", supp_thresh=supp)
    #     parse_vapor_res_somatic("/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic/sniffles2.merge.s2.vcf.somatic.s2", supp_thresh=supp)


    # # regroup
    # for sample in ["tumor", "normal"]:
    #     for seq in ["hifi", "ont", "clr"]:
    #
    #         cmd_str = "paste <(cut -f 1-7 ../release_withbed.svision_pro_v1.7.s2.vcf.somatic.s2.vapor.support.length.bed) <(cut -f 6- raw/release_withbed.svision_pro_v1.7.s2.vcf.somatic.s2.{}.{}.txt) > release_withbed.svision_pro_v1.7.s2.vcf.somatic.s2.{}.{}.txt".format(sample, seq, sample, seq)
    #         print(cmd_str)
    #
    #         # cmd_str = "paste <(cut -f 1-7 release.svision_pro_v1.7.s2.vcf.somatic.s2.vapor.bed) <(cut -f 6- raw/release.svision_pro_v1.7.s2.vcf.somatic.s2.{}.{}.txt) > release.svision_pro_v1.7.s2.vcf.somatic.s2.{}.{}.txt".format(sample, seq, sample, seq)
    #         # print(cmd_str)
    #
    #         cmd_str = "paste <(cut -f 1-7 ../sniffles2.merge.s2.vcf.somatic.s2.vapor.support.length.bed) <(cut -f 6- raw/sniffles2.merge.s2.vcf.somatic.s2.{}.{}.txt) > sniffles2.merge.s2.vcf.somatic.s2.{}.{}.txt".format(sample, seq, sample, seq)
    #         print(cmd_str)


    # print("svision-pro")
    # print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format("total_num", "inconclusive_num", "supp_1_2_3_ratio", "supp_1_ratio", "supp_2_ratio", "supp_3_ratio", "fail_ratio"))
    #
    # for supp in [2, 3, 5, 8, 10, 15, 20]:
    #     parse_vapor_res_somatic_with_repeat("/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic/modified/release_withbed.svision_pro_v1.7.s2.vcf.somatic.s2",
    #                                         "/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic/modified/repeat_annot/release_withbed.svision_pro_v1.7.s2.vcf.somatic.s2.tumor.hifi.repeat_annot.bed",
    #                                         "/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic/release_withbed.svision_pro_v1.7.s2.vcf.somatic.s2_bench_region/tp-call.bed",
    #                                         supp_thresh=supp, score_thresh=0.1)
    #
    #
    # print("sniffles")
    # print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format("total_num", "inconclusive_num", "supp_1_2_3_ratio", "supp_1_ratio", "supp_2_ratio", "supp_3_ratio", "fail_ratio"))
    #
    # for supp in [2, 3, 5, 8, 10, 15, 20]:
    #     parse_vapor_res_somatic_with_repeat("/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic/modified/sniffles2.merge.s2.vcf.somatic.s2",
    #                                         "/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic/modified/repeat_annot/sniffles2.merge.s2.vcf.somatic.s2.tumor.hifi.repeat_annot.bed",
    #                                         "/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic/sniffles2.merge.s2.vcf.somatic.s2_bench_region/tp-call.bed",
    #                                         supp_thresh=supp, score_thresh=0.1)

    print("svision-pro")

    parse_vapor_res_somatic_with_repeat_with_af("/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic/modified/release_withbed.svision_pro_v1.7.s2.vcf.somatic.s2",
                                            "/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic/modified/repeat_annot/release_withbed.svision_pro_v1.7.s2.vcf.somatic.s2.tumor.hifi.repeat_annot.bed",
                                            "/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic/release_withbed.svision_pro_v1.7.s2.vcf.somatic.s2_bench_region/tp-call.bed",
                                            "/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic/release_withbed.svision_pro_v1.7.s2.vcf.somatic.s2.vcf",

                                            supp_thresh=2, score_thresh=0.1)
    print("sniffles")

    parse_vapor_res_somatic_with_repeat_with_af("/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic/modified/sniffles2.merge.s2.vcf.somatic.s2",
                                                "/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic/modified/repeat_annot/sniffles2.merge.s2.vcf.somatic.s2.tumor.hifi.repeat_annot.bed",
                                                "/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic/sniffles2.merge.s2.vcf.somatic.s2_bench_region/tp-call.bed",
                                                "/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic/sniffles2.merge.s2.vcf.somatic.s2.vcf",

                                                supp_thresh=2, score_thresh=0.1)