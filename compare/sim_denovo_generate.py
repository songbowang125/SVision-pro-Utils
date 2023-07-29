import pysam
import random
import os

def allele_alter(sv_info):

    type = random.choice(["denovo", "inherit"])

    denovo_num = 0
    if type == "denovo":
        father_allele = 0.0
        mother_allele = 0.0
    else:
        father_allele = 0.0
        mother_allele = 0.0

        while father_allele + mother_allele == 0:
            father_allele = random.choice([0.0, 0.5, 1.0])
            mother_allele = random.choice([0.0, 0.5, 1.0])

    if father_allele == 0.0 and mother_allele == 0.0:
        # child_allele = random.choice([0.0, 0.5, 1.0])
        child_allele = random.choice([0.5])
        if child_allele == 0.5:
            denovo_num += 1
            child_out_bed0.write("\t".join(sv_info) + "\n")
        if child_allele == 1.0:
            denovo_num += 1
            child_out_bed0.write("\t".join(sv_info) + "\n")
            child_out_bed1.write("\t".join(sv_info) + "\n")

    elif father_allele == 0.5 and mother_allele == 0.5:
        father_out_bed0.write("\t".join(sv_info) + "\n")
        mother_out_bed0.write("\t".join(sv_info) + "\n")

        child_allele = random.choice([0.0, 0.5, 0.5, 1.0])

        if child_allele == 0.5:
            child_out_bed0.write("\t".join(sv_info) + "\n")
        if child_allele == 1.0:
            child_out_bed0.write("\t".join(sv_info) + "\n")
            child_out_bed1.write("\t".join(sv_info) + "\n")

    elif father_allele == 1.0 and mother_allele == 1.0:
        father_out_bed0.write("\t".join(sv_info) + "\n")
        father_out_bed1.write("\t".join(sv_info) + "\n")
        mother_out_bed0.write("\t".join(sv_info) + "\n")
        mother_out_bed1.write("\t".join(sv_info) + "\n")

        child_allele = 1.0

        child_out_bed0.write("\t".join(sv_info) + "\n")
        child_out_bed1.write("\t".join(sv_info) + "\n")

    elif 0.0 in [father_allele, mother_allele] and 1.0 in [father_allele, mother_allele]:

        if father_allele == 0.0:
            mother_out_bed0.write("\t".join(sv_info) + "\n")
            mother_out_bed1.write("\t".join(sv_info) + "\n")
        else:
            father_out_bed0.write("\t".join(sv_info) + "\n")
            father_out_bed1.write("\t".join(sv_info) + "\n")

        child_allele = 0.5

        child_out_bed0.write("\t".join(sv_info) + "\n")

    elif 0.5 in [father_allele, mother_allele] and 1.0 in [father_allele, mother_allele]:
        if father_allele == 0.5:
            mother_out_bed0.write("\t".join(sv_info) + "\n")
            mother_out_bed1.write("\t".join(sv_info) + "\n")

            father_out_bed0.write("\t".join(sv_info) + "\n")
        else:
            father_out_bed0.write("\t".join(sv_info) + "\n")
            father_out_bed1.write("\t".join(sv_info) + "\n")

            mother_out_bed0.write("\t".join(sv_info) + "\n")

        child_allele = random.choice([0.5, 1.0])

        if child_allele == 0.5:
            child_out_bed0.write("\t".join(sv_info) + "\n")

        if child_allele == 1.0:
            child_out_bed0.write("\t".join(sv_info) + "\n")
            child_out_bed1.write("\t".join(sv_info) + "\n")

    elif 0.0 in [father_allele, mother_allele] and 0.5 in [father_allele, mother_allele]:
        if father_allele == 0.0:
            mother_out_bed0.write("\t".join(sv_info) + "\n")
        else:
            father_out_bed0.write("\t".join(sv_info) + "\n")

        child_allele = random.choice([0.0, 0.5])

        if child_allele == 0.5:
            child_out_bed0.write("\t".join(sv_info) + "\n")

    else:
        child_allele = 0.0
        print("no this allele: {}-{}-{}".format(father_allele, mother_allele, child_allele))

    return "{}:{}-{}-{}\n".format(type, father_allele, mother_allele, child_allele)


def generate_somatic_benchmark(input_file):

    # # STEP: traverse
    with open(input_file) as fin, open(input_file.replace(".bed", ".denovo_allele.log"), "w") as log_out:

        for line in fin:
            sv_info = line.strip().split("\t")

            # # STEP: deal with events
            new_sv_info = allele_alter(sv_info)

            log_out.write("{}\t{}".format("\t".join(sv_info), new_sv_info))


def split_sim_benchmark_to_visor_bed(input_file, ref_file):
    ref_file = pysam.FastaFile(ref_file)
    output_file = input_file.replace(".bed", ".split.bed")

    with open(input_file) as fin, open(output_file, "w") as fout:

        for line in fin:
            line_split = line.strip().split("\t")

            chrom = line_split[0]
            start = int(line_split[1])

            end = int(line_split[2])
            sv_type = line_split[3]
            sv_type_split = sv_type.split(";")

            if len(sv_type_split) == 1:
                if len(line.strip().split("\t")) == 5:
                    fout.write(line.strip() + "\t2\n")
                else:
                    fout.write(line)
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

                    fout.write(
                        "{}\t{}\t{}\t{}\t{}\t2\n".format(chrom, dup_source_start, dup_source_end, "tandem duplication", more))

                    more = ref_file.fetch(chrom, dup_source_start, dup_source_end + 1)

                    fout.write(
                        "{}\t{}\t{}\t{}\t{}\t2\n".format(chrom, sv_bkps_str[1][1], sv_bkps_str[1][2], "insertion", more))

                elif sv_type in ["dispersed duplication;deletion", "dispersed inverted duplication;deletion"]:
                    dup_ins = sv_bkps_str[1][1] - 1
                    dup_strand = "forward" if "inverted" not in sv_type else "reverse"
                    dup_source_start = sv_bkps_str[0][1]
                    dup_source_end = sv_bkps_str[0][2]

                    more = "h1:{}:{}:{}".format(chrom, dup_ins, dup_strand)
                    fout.write(
                        "{}\t{}\t{}\t{}\t{}\t2\n".format(chrom, dup_source_start, dup_source_end, sv_bkps_str[0][0], more))

                    more = "None"
                    fout.write(
                        "{}\t{}\t{}\t{}\t{}\t2\n".format(chrom, sv_bkps_str[1][1], sv_bkps_str[1][2], sv_bkps_str[1][0], more))

                else:
                    for i in range(len(sv_bkps_str)):
                        if "duplication" in sv_bkps_str[i][0]:
                            more = "2"
                        else:
                            more = "None"

                        fout.write("{}\t{}\t{}\t{}\t{}\t2\n".format(chrom, sv_bkps_str[i][1], sv_bkps_str[i][2], sv_bkps_str[i][0], more))


if __name__ == '__main__':

    benchmark_bed = "/mnt/c/workspace/test/sim_denovo_ccs/ground_truth2/data_HG002_SVs_Tier1_v0.6.final.bed"
    child_bed0 = benchmark_bed.replace(".bed", ".child_allele0.bed")
    child_bed1 = benchmark_bed.replace(".bed", ".child_allele1.bed")

    father_bed0 = benchmark_bed.replace(".bed", ".father_allele0.bed")
    father_bed1 = benchmark_bed.replace(".bed", ".father_allele1.bed")

    mother_bed0 = benchmark_bed.replace(".bed", ".mother_allele0.bed")
    mother_bed1 = benchmark_bed.replace(".bed", ".mother_allele1.bed")

    child_out_bed0 = open(child_bed0, "w")
    child_out_bed1 = open(child_bed1, "w")
    father_out_bed0 = open(father_bed0, "w")
    father_out_bed1 = open(father_bed1, "w")
    mother_out_bed0 = open(mother_bed0, "w")
    mother_out_bed1 = open(mother_bed1, "w")

    generate_somatic_benchmark(benchmark_bed)

    child_out_bed0.close()
    child_out_bed1.close()
    father_out_bed0.close()
    father_out_bed1.close()
    mother_out_bed0.close()
    mother_out_bed1.close()


    work_path = "/data/home/songbo/workspace/svision-pro/sim_denovo_ccs/simple"
    benchmark_bed = "/data/home/songbo/workspace/svision-pro/sim_denovo_ccs/simple/data_HG002_SVs_Tier1_v0.6.final.bed"
    ref_file = "/data/home/songbo/workspace/svision-pro/sim_denovo_ccs/simple/1-Y.fa"


    # work_path = "/data/home/songbo/workspace/svision-pro/sim_denovo_ccs/complex"
    # benchmark_bed = "/data/home/songbo/workspace/svision-pro/sim_denovo_ccs/complex/data_benchmark_origin_clone0_300.csv.bed"
    # ref_file = "/data/home/songbo/workspace/svision-pro/sim_denovo_ccs/complex/chr1-X.fa"

    # benchmark_bed = "/mnt/c/workspace/test/sim_denovo/ground_truth/data_benchmark_origin_clone0_300.csv.bed"

    child_bed0 = benchmark_bed.replace(".bed", ".child_allele0.bed")
    child_bed1 = benchmark_bed.replace(".bed", ".child_allele1.bed")

    father_bed0 = benchmark_bed.replace(".bed", ".father_allele0.bed")
    father_bed1 = benchmark_bed.replace(".bed", ".father_allele1.bed")

    mother_bed0 = benchmark_bed.replace(".bed", ".mother_allele0.bed")
    mother_bed1 = benchmark_bed.replace(".bed", ".mother_allele1.bed")

    child_out_bed0 = open(child_bed0, "w")
    child_out_bed1 = open(child_bed1, "w")
    father_out_bed0 = open(father_bed0, "w")
    father_out_bed1 = open(father_bed1, "w")
    mother_out_bed0 = open(mother_bed0, "w")
    mother_out_bed1 = open(mother_bed1, "w")

    generate_somatic_benchmark(benchmark_bed)

    child_out_bed0.close()
    child_out_bed1.close()
    father_out_bed0.close()
    father_out_bed1.close()
    mother_out_bed0.close()
    mother_out_bed1.close()

    # # # STEP: split
    for i in range(2):
        child_split_bed = eval("child_bed{}".format(i)).replace(".bed", ".split.bed")
        father_split_bed = eval("father_bed{}".format(i)).replace(".bed", ".split.bed")
        mother_split_bed = eval("mother_bed{}".format(i)).replace(".bed", ".split.bed")

        child_fa = os.path.join(work_path, "child_allele{}.h1.fa".format(i))
        father_fa = os.path.join(work_path, "father_allele{}.h1.fa".format(i))
        mother_fa = os.path.join(work_path, "mother_allele{}.h1.fa".format(i))

        split_sim_benchmark_to_visor_bed(eval("child_bed{}".format(i)), ref_file)
        split_sim_benchmark_to_visor_bed(eval("father_bed{}".format(i)), ref_file)
        split_sim_benchmark_to_visor_bed(eval("mother_bed{}".format(i)), ref_file)

        cmd_str = "python /data/home/songbo/workspace/svision-pro/v0.15.4_more_somatic_types/src/post_process/sim/simulate.py sim_clone -g {} -c {} -bed {} -o {}".format(ref_file, "child_allele{}".format(i), child_split_bed, work_path)
        os.system(cmd_str)

        cmd_str = "python /data/home/songbo/workspace/svision-pro/v0.15.4_more_somatic_types/src/post_process/sim/simulate.py sim_clone -g {} -c {} -bed {} -o {}".format(ref_file, "father_allele{}".format(i), father_split_bed, work_path)
        os.system(cmd_str)

        cmd_str = "python /data/home/songbo/workspace/svision-pro/v0.15.4_more_somatic_types/src/post_process/sim/simulate.py sim_clone -g {} -c {} -bed {} -o {}".format(ref_file, "mother_allele{}".format(i), mother_split_bed, work_path)
        os.system(cmd_str)

        # # generate reads
        cmd_str = "/data/home/songbo/workspace/tools/pbsim2-master/src/pbsim --seed 1234 --depth 30 --sample-fastq /data/DATA/GIAB/RAWDATA/HiFi/HG002/PacBio_CCS_15kb/reads/m54238_180901_011437.Q20.fastq {}".format(child_fa)
        # cmd_str = "/data/home/songbo/workspace/tools/pbsim2-master/src/pbsim --seed 1234 --depth 30 --hmm_model /data/home/songbo/workspace/tools/pbsim2-master/data/R95.model  {}".format(child_fa)
        os.system(cmd_str)
        os.system("cat sd*fastq > {}".format(os.path.join(work_path, "child_allele{}.fastq".format(i))))
        os.system("rm sd*ref")
        os.system("rm sd*maf")
        os.system("rm sd*fastq")

        cmd_str = "/data/home/songbo/workspace/tools/pbsim2-master/src/pbsim --seed 1234 --depth 30 --sample-fastq /data/DATA/GIAB/RAWDATA/HiFi/HG002/PacBio_CCS_15kb/reads/m54238_180901_011437.Q20.fastq {}".format(father_fa)
        # cmd_str = "/data/home/songbo/workspace/tools/pbsim2-master/src/pbsim --seed 1234 --depth 30 --hmm_model /data/home/songbo/workspace/tools/pbsim2-master/data/R95.model  {}".format(father_fa)
        os.system(cmd_str)
        os.system("cat sd*fastq > {}".format(os.path.join(work_path, "father_allele{}.fastq".format(i))))
        os.system("rm sd*ref")
        os.system("rm sd*maf")
        os.system("rm sd*fastq")

        cmd_str = "/data/home/songbo/workspace/tools/pbsim2-master/src/pbsim --seed 1234 --depth 30 --sample-fastq /data/DATA/GIAB/RAWDATA/HiFi/HG002/PacBio_CCS_15kb/reads/m54238_180901_011437.Q20.fastq {}".format(mother_fa)
        # cmd_str = "/data/home/songbo/workspace/tools/pbsim2-master/src/pbsim --seed 1234 --depth 30 --hmm_model /data/home/songbo/workspace/tools/pbsim2-master/data/R95.model  {}".format(mother_fa)
        os.system(cmd_str)
        os.system("cat sd*fastq > {}".format(os.path.join(work_path, "mother_allele{}.fastq".format(i))))
        os.system("rm sd*ref")
        os.system("rm sd*maf")
        os.system("rm sd*fastq")

    child_final_fastq = open(os.path.join(work_path, "child.fastq"), "w")
    father_final_fastq = open(os.path.join(work_path, "father.fastq"), "w")
    mother_final_fastq = open(os.path.join(work_path, "mother.fastq"), "w")

    for i in range(2):
        for line in open(os.path.join(work_path, "child_allele{}.fastq".format(i))):
            if len(line) < 100 and (line.startswith("@S") or line.startswith("+S")):
                child_final_fastq.write(line.strip() + "_allele{}\n".format(i))
            else:
                child_final_fastq.write(line)

        for line in open(os.path.join(work_path, "father_allele{}.fastq".format(i))):
            if len(line) < 100 and (line.startswith("@S") or line.startswith("+S")):
                father_final_fastq.write(line.strip() + "_allele{}\n".format(i))
            else:
                father_final_fastq.write(line)

        for line in open(os.path.join(work_path, "mother_allele{}.fastq".format(i))):
            if len(line) < 100 and (line.startswith("@S") or line.startswith("+S")):
                mother_final_fastq.write(line.strip() + "_allele{}\n".format(i))
            else:
                mother_final_fastq.write(line)


    child_final_fastq.close()
    father_final_fastq.close()
    mother_final_fastq.close()