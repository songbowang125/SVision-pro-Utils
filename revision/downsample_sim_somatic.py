import os
import pysam
import random

def collect_all_reads(origin_bam_path):

    read_list = []

    origin_bam = pysam.AlignmentFile(origin_bam_path)

    for align in origin_bam:

        if align.qname not in read_list:
            read_list.append(align.qname)

    origin_bam.close()

    return read_list


def downsample_bam(origin_bam_path, candidate_read_list, candidate_coverage):

    print("downsampling ", candidate_coverage)
    origin_bam = pysam.AlignmentFile(origin_bam_path)

    output_bam_path = origin_bam_path.replace(".bam", ".tmp.bam".format(candidate_coverage))

    with pysam.AlignmentFile(output_bam_path, "wb", header=origin_bam.header) as output_bam:

        for align in origin_bam:

            if align.qname in candidate_read_list:
                output_bam.write(align)

    output_srt_bam_path = origin_bam_path.replace(".bam", ".{}.bam".format(candidate_coverage))
    os.system("samtools sort -@ 8 -o {} {}".format(output_srt_bam_path, output_bam_path))

    os.remove(output_bam_path)

if __name__ == '__main__':
    origin_bam_path = "/mnt/d/workspace/svision-pro/sim_somatic_ccs/simple_100/changed_reads.srt.bam"

    read_list = collect_all_reads(origin_bam_path)

    read_list_80x = random.sample(read_list, int(len(read_list) * 0.8))

    read_list_50x = random.sample(read_list_80x, int(len(read_list_80x) * 0.625))

    read_list_30x = random.sample(read_list_50x, int(len(read_list_50x) * 0.6))

    read_list_20x = random.sample(read_list_30x, int(len(read_list_30x) * 0.667))

    read_list_10x = random.sample(read_list_20x, int(len(read_list_20x) * 0.5))

    downsample_bam(origin_bam_path, read_list_80x, "80x")
    downsample_bam(origin_bam_path, read_list_50x, "50x")
    downsample_bam(origin_bam_path, read_list_30x, "30x")
    downsample_bam(origin_bam_path, read_list_20x, "20x")
    downsample_bam(origin_bam_path, read_list_10x, "10x")

    print(len(read_list), len(read_list_80x), len(read_list_50x), len(read_list_30x), len(read_list_20x), len(read_list_10x))