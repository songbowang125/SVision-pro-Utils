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


def load_from_split_bed(split_bed_file, out_path, candidate_ids=None, candidate_ids_vafs=None):
    print(candidate_ids)
    csv_list = []
    id = 0
    for line in open(split_bed_file):
        id += 1

        if candidate_ids is not None and id not in candidate_ids:
            continue

        line_split = line.strip().split("\t")

        ssv_ref_chrom = line_split[0]
        ssv_ref_start = int(line_split[1])
        ssv_ref_end = int(line_split[2])
        ssv_type = line_split[3]

        ssv_info = line_split[4]

        csv = CSV(id, [SSV(ssv_ref_chrom, ssv_ref_start, ssv_ref_end, ssv_type, ssv_info)])

        if candidate_ids_vafs is not None:
            csv.vaf = candidate_ids_vafs[candidate_ids.index(id)]

        csv_list.append(csv)

        print(candidate_ids_vafs[candidate_ids.index(id)], csv.vaf)

    with open(os.path.join(out_path, "csv_list.txt"), "w") as fout:
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
        if ssv.ssv_type == "insertion":
            alter_seq = alter_seq[: ssv_start_index] + ssv.ssv_info + alter_seq[ssv_end_index: ]

        elif ssv.ssv_type == "deletion":
            alter_seq = alter_seq[: ssv_start_index] + alter_seq[ssv_end_index: ]

        elif ssv.ssv_type == "inversion":
            alter_seq = alter_seq[: ssv_start_index] + reverse_seq(alter_seq[ssv_start_index: ssv_end_index]) + alter_seq[ssv_end_index: ]

        elif ssv.ssv_type in ["dispersed inverted duplication", "dispersed duplication"]:
            source_seq = ref_file.fetch(ssv.ssv_info[0], ssv.ssv_info[1], ssv.ssv_info[2] + 1)

            if "inverted" in ssv.ssv_type:
                source_seq = reverse_seq(source_seq)

            alter_seq = alter_seq[: ssv_start_index] + source_seq + alter_seq[ssv_end_index: ]

        elif ssv.ssv_type in ["inverted tandem duplication", "tandem duplication"]:
            source_start_index = ssv.ssv_info[1] - read.ref_start
            source_end_index = ssv.ssv_info[2] + 1 - read.ref_start

            # source_seq = ref_file.fetch(ssv.ssv_info[0], ssv.ssv_info[1], ssv.ssv_info[2] + 1)
            source_seq = alter_seq[source_start_index: source_end_index]

            if "inverted" in ssv.ssv_type:
                source_seq = reverse_seq(source_seq)

            alter_seq = alter_seq[: ssv_start_index] + source_seq + alter_seq[ssv_end_index:]

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

            if len(csv_altered_reads) != csv_total_alter_read_num:
                print(csv.id, csv.vaf, len(csv_altered_reads))


        for chrom in all_altered_read_names:
            print(chrom)
            for align in raw_bam_file.fetch(chrom):
                if align.qname not in all_altered_read_names[chrom]:
                    out_bam_file.write(align)


def down_sample_to_100(all_sv_list_path, output_path):

    # # STEP: load all sv list file
    line_num = 0

    all_sv_collect = {}

    all_sv_list_file = open(all_sv_list_path)

    for line in all_sv_list_file:

        line_split = line.strip().split("\t")

        line_vaf = line_split[-1]

        if line_vaf not in all_sv_collect:
            all_sv_collect[line_vaf] = []

        all_sv_collect[line_vaf].append(line_num)

        line_num += 1

    all_sv_list_file.close()

    # # down sample
    import random

    left_lines = []
    left_ids = []
    left_ids_vafs = []

    for vaf in all_sv_collect:
        left_lines.extend(random.sample(all_sv_collect[vaf], 100))

    # # output after down sample
    output_file =open(os.path.join(output_path, "sv_list_100.txt"), "w")

    line_num = 0

    all_sv_list_file = open(all_sv_list_path)

    for line in all_sv_list_file:
        if line_num in left_lines:
            left_ids.append(int(line.split("\t")[0]))
            left_ids_vafs.append(float(line.split("\t")[-1]))
            output_file.write(line)

        line_num += 1

    all_sv_list_file.close()
    output_file.close()


    return left_ids, left_ids_vafs



def load_from_down_100(sv_list_100_path):
    candidate_ids = []
    candidate_ids_vafs = []

    for line in open(sv_list_100_path):
        line_split = line.strip().split("\t")

        candidate_ids.append(int(line_split[0]))

        candidate_ids_vafs.append(float(line_split[-1]))

    return candidate_ids, candidate_ids_vafs

if __name__ == '__main__':

    #
    # ref_file = pysam.FastaFile("/data/home/songbo/workspace/ref/hs37d5.1-Y.fa")
    # split_bed_file = "/data/home/songbo/workspace/svision-pro/sim_somatic_ont/simple_bamsurgeon/data_HG002_SVs_Tier1_v0.6.final.bed"
    #
    # raw_bam_path = "/data/home/songbo/workspace/svision-pro/sim_somatic_ont/simple/simulated.100x.srt.bam"
    #
    # output_path = "/data/home/songbo/workspace/svision-pro/sim_somatic_ont/simple_bamsurgeon"
    #
    #
    # csv_list = load_from_split_bed(split_bed_file, output_path)
    #
    # print(len(csv_list))
    #
    # surgone_raw_bam(csv_list, raw_bam_path, output_path, ref_file)

    # # down to 100
    # raw_bam_path = "/data/home/songbo/workspace/svision-pro/sim_somatic_ccs/simple_100_80x/simulated.80x.srt.bam"
    #
    # ref_file = pysam.FastaFile("/data/home/songbo/workspace/ref/hs37d5.1-Y.fa")
    #
    # split_bed_file = "/data/home/songbo/workspace/svision-pro/sim_somatic_ccs/simple/data_HG002_SVs_Tier1_v0.6.final.bed"
    #
    # output_path = "/data/home/songbo/workspace/svision-pro/sim_somatic_ccs/simple_100_80x"
    #
    # all_sv_list_path = "/data/home/songbo/workspace/svision-pro/sim_somatic_ccs/simple/sv_list.txt"
    #
    # candidate_ids, candidate_ids_vafs = down_sample_to_100(all_sv_list_path, output_path)
    #
    #
    # csv_list = load_from_split_bed(split_bed_file, output_path, candidate_ids, candidate_ids_vafs)
    #
    # print(len(csv_list))
    #
    # surgone_raw_bam(csv_list, raw_bam_path, output_path, ref_file, raw_bam_coverage=80)


    # # down to low coverage
    raw_bam_coverage = 80

    raw_bam_path = "/data/home/songbo/workspace/svision-pro/sim_somatic_ccs/simple_100_{}x/simulated.{}x.srt.bam".format(raw_bam_coverage, raw_bam_coverage)

    output_path = "/data/home/songbo/workspace/svision-pro/sim_somatic_ccs/simple_100_{}x".format(raw_bam_coverage)

    ref_file = pysam.FastaFile("/data/home/songbo/workspace/ref/hs37d5.1-Y.fa")

    split_bed_file = "/data/home/songbo/workspace/svision-pro/sim_somatic_ccs/simple/data_HG002_SVs_Tier1_v0.6.final.bed"

    sv_list_100_path = "/data/home/songbo/workspace/svision-pro/sim_somatic_ccs/simple_100/sv_list_100.txt"

    candidate_ids, candidate_ids_vafs = load_from_down_100(sv_list_100_path)

    csv_list = load_from_split_bed(split_bed_file, output_path, candidate_ids, candidate_ids_vafs)

    surgone_raw_bam(csv_list, raw_bam_path, output_path, ref_file, raw_bam_coverage=raw_bam_coverage)

    os.system("cd {}".format(output_path))
    os.system("minimap2 -a -t 20 --MD -Y -R \'@RG\\tID:unchanged\\tSM:unchanged\\tLB:unchanged\' -x map-pb /data/home/songbo/workspace/ref/hs37d5.1-Y.fa changed_reads.fastq > changed_reads.sam")
    os.system("samtools sort -@ 20 -o changed_reads.srt.bam changed_reads.sam")
    os.system("samtools index -@ 20  changed_reads.srt.bam")
    os.system("samtools merge -@ 20 -o merged.srt.bam changed_reads.srt.bam unchanged_bam.bam")
    os.system("samtools index -@ 20 merged.srt.bam")
