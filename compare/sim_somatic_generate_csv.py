import os
import random

import pysam


class CSV:
    def __init__(self, id, included_SSVs):

        self.id = id

        self.ref_chrom = included_SSVs[0].ref_chrom
        self.ref_start = min([ssv.ref_start for ssv in included_SSVs])
        self.ref_end = max([ssv.ref_end for ssv in included_SSVs])

        self.included_SSVs = included_SSVs

        self.random_set_vaf()

    def random_set_vaf(self):
        self.vaf = random.choice([0.01, 0.02, 0.03, 0.04, 0.05, 0.08, 0.10])

    def to_string(self):
        return "{}\t{}\t{}\t{}\t{}\t{}".format(self.id, self.ref_chrom, self.ref_start, self.ref_end, ";".join([ssv.to_string() for ssv in self.included_SSVs]), self.vaf)

class SSV:
    def __init__(self, ref_chrom, ref_start, ref_end, ssv_type, ssv_info):

        self.ref_chrom = ref_chrom
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.ssv_type = ssv_type
        self.ssv_info = ssv_info

    def to_string(self):

        return "{}_{}_{}_{}".format(self.ref_chrom, self.ref_start, self.ref_end, self.ssv_type)


class READ:
    def __init__(self, ref_chrom, ref_start, ref_end, read_name, read_seq):

        self.ref_chrom = ref_chrom
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.read_seq = read_seq
        self.read_name = read_name


def split_sim_benchmark_to_visor_bed(input_file, ref_file):

    ref_file = pysam.FastaFile(ref_file)

    output_file = input_file.replace(".bed", ".split.bed")

    id = 0
    with open(input_file) as fin, open(output_file, "w") as fout:

        for line in fin:
            id += 1
            line_split = line.strip().split("\t")

            chrom = line_split[0]
            start = int(line_split[1])
            end = int(line_split[2])
            sv_type = line_split[3]
            sv_type_split = sv_type.split(";")

            # # is a ssv
            if len(sv_type_split) == 1:
                if len(line.strip().split("\t")) == 5:
                    fout.write("{}\t".format(id) + line.strip() + "\t2\n")

            else:
                sv_bkps_str = line_split[4].split("><")
                for i in range(len(sv_bkps_str)):
                    sv_bkps_str[i] = sv_bkps_str[i].replace("<", "")
                    sv_bkps_str[i] = sv_bkps_str[i].replace(">", "")
                    sv_bkps_str[i] = sv_bkps_str[i].split(",")
                    sv_bkps_str[i][1] = int(sv_bkps_str[i][1])
                    sv_bkps_str[i][2] = int(sv_bkps_str[i][2])

                if sv_type in ["tandem duplication;insertion"]:

                    dup_source_start = sv_bkps_str[0][1]
                    dup_source_end = sv_bkps_str[0][2]
                    more = "2"

                    fout.write("{}\t{}\t{}\t{}\t{}\t{}\t2\n".format(id, chrom, dup_source_start, dup_source_end, "tandem duplication", more))

                    more = ref_file.fetch(chrom, dup_source_start, dup_source_end + 1)

                    fout.write("{}\t{}\t{}\t{}\t{}\t{}\t2\n".format(id, chrom, sv_bkps_str[1][1], sv_bkps_str[1][2], "insertion", more))

                elif sv_type in ["dispersed duplication;deletion", "dispersed inverted duplication;deletion"]:
                    dup_ins = sv_bkps_str[1][1] - 1
                    dup_strand = "forward" if "inverted" not in sv_type else "reverse"
                    dup_source_start = sv_bkps_str[0][1]
                    dup_source_end = sv_bkps_str[0][2]

                    more = "h1:{}:{}:{}".format(chrom, dup_ins, dup_strand)
                    fout.write("{}\t{}\t{}\t{}\t{}\t{}\t2\n".format(id, chrom, dup_source_start, dup_source_end, sv_bkps_str[0][0], more))

                    more = "None"
                    fout.write("{}\t{}\t{}\t{}\t{}\t{}\t2\n".format(id, chrom, sv_bkps_str[1][1], sv_bkps_str[1][2], sv_bkps_str[1][0], more))

                else:
                    for i in range(len(sv_bkps_str)):
                        if "duplication" in sv_bkps_str[i][0]:
                            more = "2"
                        else:
                            more = "None"

                        fout.write("{}\t{}\t{}\t{}\t{}\t{}\t2\n".format(id, chrom, sv_bkps_str[i][1], sv_bkps_str[i][2], sv_bkps_str[i][0], more))


def load_from_split_bed(split_bed_file, out_path):

    csv_list = []

    previous_csv_id = -1
    included_SSVs = []

    for line in open(split_bed_file):
        line_split = line.strip().split("\t")

        ssv_ref_chrom = line_split[1]
        ssv_ref_start = int(line_split[2])
        ssv_ref_end = int(line_split[3])
        ssv_type = line_split[4]

        current_csv_id = line_split[0]

        if current_csv_id != previous_csv_id:

            if previous_csv_id != -1:
                csv_list.append(CSV(previous_csv_id, included_SSVs))

            previous_csv_id = current_csv_id
            included_SSVs = []

            if "dispersed" in ssv_type:
                dup_info = line_split[5].split(":")

                ins_chrom = dup_info[1]
                ins_start = int(dup_info[2])
                ins_end = ins_start

                included_SSVs.append(SSV(ins_chrom, ins_start, ins_end, ssv_type, [ssv_ref_chrom, ssv_ref_start, ssv_ref_end]))

            elif "tandem" in ssv_type:
                included_SSVs.append(SSV(ssv_ref_chrom, ssv_ref_end, ssv_ref_end, ssv_type, [ssv_ref_chrom, ssv_ref_start, ssv_ref_end]))

            else:
                included_SSVs.append(SSV(ssv_ref_chrom, ssv_ref_start, ssv_ref_end, ssv_type, None))

        else:
            if "dispersed" in ssv_type:
                dup_info = line_split[5].split(":")

                ins_chrom = dup_info[1]
                ins_start = int(dup_info[2])
                ins_end = ins_start + 1

                included_SSVs.append(SSV(ins_chrom, ins_start, ins_end, ssv_type, [ssv_ref_chrom, ssv_ref_start, ssv_ref_end]))
            else:
                included_SSVs.append(SSV(ssv_ref_chrom, ssv_ref_start, ssv_ref_end, ssv_type, None))

    with open(os.path.join(output_path, "csv_list.txt"), "w") as fout:
        for csv in csv_list:
            fout.write(csv.to_string() + "\n")

    return csv_list


def reverse_seq(seq):

    reversed_seq = ""

    for i in range(len(seq) - 1, -1, -1):
        if seq[i] == "A" or seq[i] == "a":
            reversed_seq += "T"
        elif seq[i] == "T" or seq[i] == "t":
            reversed_seq += "A"
        elif seq[i] == "C" or seq[i] == "c":
            reversed_seq += "G"
        elif seq[i] == "G" or seq[i] == "g":
            reversed_seq += "C"
        else:
            reversed_seq += "N"

    return reversed_seq


def alter_read(csv, read, ref_file):


    alter_seq = read.read_seq

    for ssv_index in range(len(csv.included_SSVs) - 1, -1, -1):
        ssv = csv.included_SSVs[ssv_index]

        ssv_start_index = ssv.ref_start - read.ref_start
        ssv_end_index = ssv.ref_end - read.ref_start

        if ssv.ssv_type == "deletion":
            alter_seq = alter_seq[: ssv_start_index] + alter_seq[ssv_end_index: ]

        elif ssv.ssv_type == "inversion":
            alter_seq = alter_seq[: ssv_start_index] + reverse_seq(alter_seq[ssv_start_index: ssv_end_index]) + alter_seq[ssv_end_index: ]

        elif ssv.ssv_type in ["dispersed inverted duplication", "dispersed duplication", "inverted tandem duplication", "tandem duplication"]:
            source_seq = ref_file.fetch(ssv.ssv_info[0], ssv.ssv_info[1], ssv.ssv_info[2] + 1)

            if "inverted" in ssv.ssv_type:
                source_seq = reverse_seq(source_seq)

            alter_seq = alter_seq[: ssv_start_index] + source_seq + alter_seq[ssv_end_index: ]
        else:
            print("no this sv type", ssv.ssv_type)
            exit()


    return alter_seq


def surgone_raw_bam(csv_list, raw_bam_path, output_path, ref_file, raw_bam_coverage=100):

    raw_bam_file = pysam.AlignmentFile(raw_bam_path)

    out_bam_path = os.path.join(output_path, "unchanged_bam.bam")
    out_fastq_path = os.path.join(output_path, "changed_reads.fastq")

    with pysam.AlignmentFile(out_bam_path, "wb", header=raw_bam_file.header) as out_bam_file, open(out_fastq_path, "w") as out_fastq_file:

        all_altered_read_names = {}

        # # traverse each csv
        for csv in csv_list:

            # if csv.ref_chrom != "chr20":
            #     continue

            csv_total_alter_read_num = int(raw_bam_coverage * csv.vaf)

            csv_altered_reads = []
            for align in raw_bam_file.fetch(csv.ref_chrom, csv.ref_start - 50, csv.ref_end + 50):

                if len(csv_altered_reads) >= csv_total_alter_read_num:
                    break

                # # filter aligns
                if align.is_supplementary or align.cigarstring is None or align.is_unmapped or align.is_secondary or align.mapq < 20:
                    continue

                try:
                    sa_tags = align.get_tag("SA").split(";")
                except KeyError:
                    sa_tags = []

                # # must contain no SA
                if len(sa_tags) != 0:
                    continue

                if align.is_reverse:
                    continue

                if align.reference_start < csv.ref_start - 1000 and align.reference_end > csv.ref_end + 1000:
                    read = READ(align.reference_name, align.reference_start, align.reference_end, align.qname, align.query_sequence)
                    csv_altered_reads.append(read)

                    # # record this read name
                    if align.reference_name not in all_altered_read_names:
                        all_altered_read_names[align.reference_name] = []
                    all_altered_read_names[align.reference_name].append(align.qname)

                    # # output altered read
                    out_fastq_file.write(">{}\n".format(align.qname))
                    out_fastq_file.write("{}\n".format(alter_read(csv, read, ref_file)))

            print(csv.id, csv.vaf, len(csv_altered_reads))


        for chrom in all_altered_read_names:
            print(chrom)
            for align in raw_bam_file.fetch(chrom):
                if align.qname not in all_altered_read_names[chrom]:
                    out_bam_file.write(align)

if __name__ == '__main__':

    # ref_file = "/mnt/d/data/ref/grch38/GRCh38.d1.vd1.fa"
    #
    # bed_file = "/mnt/c/workspace/test/sim_somatic_ccs/complex_bamsurgeon/data_benchmark_origin_clone0_300.csv.bed"
    #
    # split_bed_file = "/mnt/c/workspace/test/sim_somatic_ccs/complex_bamsurgeon/benchmark_origin_clone0_300.csv.split.bed"


    ref_file = pysam.FastaFile("/data/home/songbo/workspace/ref/chr1-X.fa")
    split_bed_file = "/data/home/songbo/workspace/svision-pro/sim_somatic_ccs/complex_bamsurgeon/benchmark_origin_clone0_300.csv.split.bed"

    raw_bam_path = "/data/home/songbo/workspace/svision-pro/sim_somatic_ccs/complex/simulated.100x.srt.bam"

    output_path = "/data/home/songbo/workspace/svision-pro/sim_somatic_ccs/complex_bamsurgeon"


    csv_list = load_from_split_bed(split_bed_file, output_path)

    print(len(csv_list))

    surgone_raw_bam(csv_list, raw_bam_path, output_path, ref_file)