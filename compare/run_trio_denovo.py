import os
import pysam


def fetch_denovo_sv(input_vcf, tool):

    output_path = "/".join(input_vcf.split("/")[: -1]) + "/trio_denovo"

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    print(os.path.join(output_path, input_vcf.split("/")[-1].replace(".vcf", ".denovo.vcf")))
    output_vcf = open(os.path.join(output_path, input_vcf.split("/")[-1].replace(".vcf", ".denovo.vcf")), "w")

    input_vcf = pysam.VariantFile(input_vcf)

    output_vcf.write(str(input_vcf.header))

    for record in input_vcf:

        if record.contig not in ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", ]:
            continue

        if "BND" in str(record):
            continue

        # if "PASS" not in str(record):
        #     continue

        if abs(record.info["SVLEN"]) > 50000:
            continue

        if tool == "svision-pro" and "_NewCOMP_NewCOMP;" in str(record):
            output_vcf.write(str(record))

        if tool == "merge" and ("SUPP_VEC=100;" in str(record) or "SUPP_VEC=1000;" in str(record) or "SUPP_VEC=1100;" in str(record)):
            output_vcf.write(str(record))

        if tool == "sniffles2" and ("SUPP_VEC=100\t" in str(record) or "SUPP_VEC=1000\t" in str(record) or "SUPP_VEC=1100\t" in str(record)):
            output_vcf.write(str(record))

    input_vcf.close()
    output_vcf.close()


def collect_and_summary_denovo(input_path):

    svision_pro_summary = open(os.path.join(input_path, "collection.svision_pro.denovo.bed"), "w")

    sniffles2_summary = open(os.path.join(input_path, "collection.sniffles2.denovo.bed"), "w")
    merge_summary = open(os.path.join(input_path, "collection.merge.denovo.bed"), "w")

    for vcf_file in os.listdir(input_path):
        if ".denovo.vcf" not in vcf_file:
            continue

        sample_name = vcf_file.split(".")[1]
        tool_name = vcf_file.split(".")[0]

        vcf_file = pysam.VariantFile(os.path.join(input_path, vcf_file))
        for record in vcf_file:

            # if "IMPRECISE" in str(record):
            #     continue

            record_contig = record.contig
            record_start = record.start + 1
            record_end = record.stop + 1
            record_svtype = record.info["SVTYPE"]
            record_length = abs(int(record.info["SVLEN"]))

            if record_length > 50000 or record_length < 50:
                continue
            if "release_" in tool_name:
                svision_pro_summary.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(record_contig, record_start, record_end, record_length, record_svtype, tool_name.split("_")[1], "svision_pro", str(record).split("\t")[6]))
            elif "sniffles2" in tool_name:
                sniffles2_summary.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(record_contig, record_start, record_end, record_length, record_svtype, sample_name, tool_name, str(record).split("\t")[7].split(";")[0]))
            else:
                merge_summary.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(record_contig, record_start, record_end, record_length, record_svtype, sample_name, tool_name, str(record).split("\t")[6]))

        vcf_file.close()
    svision_pro_summary.close()
    merge_summary.close()
    sniffles2_summary.close()


def fetch_bam(input_file, output_path):

    bam_dict = {"LCL5": "/data/DATA/Human_genomes/reference_based_analysis/aligned_reads/ChineseQuartet/LCL5/ChineseQuartet.LCL5.GRCh38.HiFi.minimap2.bam",
                "LCL6": "/data/DATA/Human_genomes/reference_based_analysis/aligned_reads/ChineseQuartet/LCL6/ChineseQuartet.LCL6.GRCh38.HiFi.minimap2.bam",
                "LCL7": "/data/DATA/Human_genomes/reference_based_analysis/aligned_reads/ChineseQuartet/LCL7/ChineseQuartet.LCL7.GRCh38.HiFi.minimap2.bam",
                "LCL8": "/data/DATA/Human_genomes/reference_based_analysis/aligned_reads/ChineseQuartet/LCL8/ChineseQuartet.LCL8.GRCh38.HiFi.minimap2.bam",
                "NA19238": "/data/DATA/Human_genomes/reference_based_analysis/aligned_reads/HGSVC/NA19238/HGSVC.NA19238.GRCh38.HiFi.minimap2.bam",
                "NA19239": "/data/DATA/Human_genomes/reference_based_analysis/aligned_reads/HGSVC/NA19239/HGSVC.NA19239.GRCh38.HiFi.minimap2.bam",
                "NA19240": "/data/DATA/Human_genomes/reference_based_analysis/aligned_reads/HGSVC/NA19240/HGSVC.NA19240.GRCh38.HiFi.minimap2.bam",
                "HG00512": "/data/DATA/Human_genomes/reference_based_analysis/aligned_reads/HGSVC/HG00512/HGSVC.HG00512.GRCh38.HiFi.minimap2.bam",
                "HG00513": "/data/DATA/Human_genomes/reference_based_analysis/aligned_reads/HGSVC/HG00513/HGSVC.HG00513.GRCh38.HiFi.minimap2.bam",
                "HG00514": "/data/DATA/Human_genomes/reference_based_analysis/aligned_reads/HGSVC/HG00514/HGSVC.HG00514.GRCh38.HiFi.minimap2.bam",
                "HG00731": "/data/DATA/Human_genomes/reference_based_analysis/aligned_reads/HGSVC/HG00731/HGSVC.HG00731.GRCh38.HiFi.minimap2.bam",
                "HG00732": "/data/DATA/Human_genomes/reference_based_analysis/aligned_reads/HGSVC/HG00732/HGSVC.HG00732.GRCh38.HiFi.minimap2.bam",
                "HG00733": "/data/DATA/Human_genomes/reference_based_analysis/aligned_reads/HGSVC/HG00733/HGSVC.HG00733.GRCh38.HiFi.minimap2.bam",
                "HG002": "/data/DATA/Human_genomes/reference_based_analysis/aligned_reads/GIAB/HG002/GIAB.HG002.GRCh38.HiFi.minimap2.bam",
                "HG003": "/data/DATA/Human_genomes/reference_based_analysis/aligned_reads/GIAB/HG003/GIAB.HG003.GRCh38.HiFi.minimap2.bam",
                "HG004": "/data/DATA/Human_genomes/reference_based_analysis/aligned_reads/GIAB/HG004/GIAB.HG004.GRCh38.HiFi.minimap2.bam",
                "HG005": "/data/DATA/Human_genomes/reference_based_analysis/aligned_reads/GIAB/HG005/GIAB.HG005.GRCh38.HiFi.minimap2.bam",
                "HG006": "/data/DATA/Human_genomes/reference_based_analysis/aligned_reads/GIAB/HG006/GIAB.HG006.GRCh38.HiFi.minimap2.bam",
                "HG007": "/data/DATA/Human_genomes/reference_based_analysis/aligned_reads/GIAB/HG007/GIAB.HG007.GRCh38.HiFi.minimap2.bam", }

    for line in open(input_file):
        line_split = line.strip().split("\t")
        print(line.strip())
        chrom = line_split[0]
        start = int(line_split[1])
        end = int(line_split[2])

        sample = line_split[5]

        if "hg002-" in sample:
            bam_list = [bam_dict["HG002"], bam_dict["HG003"], bam_dict["HG004"]]
        elif "hg005-" in sample:
            bam_list = [bam_dict["HG005"], bam_dict["HG006"], bam_dict["HG007"]]
        elif "hg00514-" in sample:
            bam_list = [bam_dict["HG00514"], bam_dict["HG00512"], bam_dict["HG00513"]]
        elif "hg00733-" in sample:
            bam_list = [bam_dict["HG00731"], bam_dict["HG00732"], bam_dict["HG00733"]]
        elif "na19240-" in sample:
            bam_list = [bam_dict["NA19238"], bam_dict["NA19239"], bam_dict["NA19240"]]
        elif "LCL5-6-7-8" in sample:
            bam_list = [bam_dict["LCL5"], bam_dict["LCL6"], bam_dict["LCL7"], bam_dict["LCL8"]]
        else:
            print(sample)
            exit()


        for bam_path in bam_list:
            out_bam_path = os.path.join(output_path, "{}-{}-{}.".format(chrom, start, end) + bam_path.split("/")[-1])

            cmd_str = "samtools view -@ 20 -b -h {} {}:{}-{} > {}".format(bam_path, chrom, start - 10000, end + 10000, out_bam_path)
            os.system(cmd_str)

            cmd_str = "samtools index -@ 20 {}".format(out_bam_path)
            os.system(cmd_str)


def convert_svision_pro_to_merged_bed(input_path):
    trios = ["hg002-3-4", "hg005-6-7", "na19240-38-39", "hg00514-2-3", "hg00733-1-2", "lcl5-6-7-8"]

    with open(os.path.join(input_path, "release_six-trios_256.svision_pro_v1.6.s10.bed"), "w") as fout:
        for trio in trios:
            vcf_path = os.path.join(input_path, "release_{}_256.svision_pro_v1.6.s10.vcf".format(trio))

            for record in pysam.VariantFile(vcf_path):
                chrom = record.contig
                start = record.pos + 1
                end = record.stop + 2
                svlength = record.info["SVLEN"]
                svtype = record.info['SVTYPE']

                bkp = record.info["BKPSIT"]

                bed_str = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(chrom, start, end, svlength, svtype, trio, ",".join(list(bkp)))
                fout.write(bed_str)



if __name__ == '__main__':
    # convert_svision_pro_to_merged_bed()

    #
    trios = ["hg002-3-4", "hg005-6-7", "na19240-38-39", "hg00514-2-3", "hg00733-1-2", "lcl5-6-7-8"]

    # trios = ["hg002-3-4_ONT", "hg005-6-7_ONT", "lcl5-6-7-8_ONT"]
    # #
    print("------------SVision-pro")
    # # convert_svision_pro_to_merged_bed("/mnt/c/workspace/test_new/real_trio/")
    #
    for trio in trios:
        fetch_denovo_sv("/mnt/d/workspace/svision-pro/real_trio/release_{}_256.svision_pro_v1.7.s10.vcf".format(trio), tool="svision-pro")

    # print("------------sniffles2")
    # for trio in trios:
    #     fetch_denovo_sv("/mnt/d/workspace/svision-pro/real_trio/supp10_supp2/sniffles2.{}.s10.vcf".format(trio), tool="sniffles2")

    # print("------------svision_jasmine")
    # for trio in trios:
    #     fetch_denovo_sv("/mnt/d/workspace/svision-pro/real_trio/supp10_supp2/svision_jasmine.{}.s10.s2.vcf".format(trio), tool="merge")
    # #
    # print("------------svision_survivor")
    # for trio in trios:
    #     fetch_denovo_sv("/mnt/d/workspace/svision-pro/real_trio/supp10_supp2/svision_survivor.{}.s10.s2.vcf".format(trio), tool="merge")
    #
    # print("------------cutesv_jasmine")
    # for trio in trios:
    #     fetch_denovo_sv("/mnt/d/workspace/svision-pro/real_trio/supp10_supp2/cutesv_jasmine.{}.s10.s2.vcf".format(trio), tool="merge")
    #
    # print("------------cutesv_survivor")
    # for trio in trios:
    #     fetch_denovo_sv("/mnt/d/workspace/svision-pro/real_trio/supp10_supp2/cutesv_survivor.{}.s10.s2.vcf".format(trio), tool="merge")
    #
    # print("------------debreak_jasmine")
    # for trio in trios:
    #     fetch_denovo_sv("/mnt/d/workspace/svision-pro/real_trio/supp10_supp2/debreak_jasmine.{}.s10.s2.vcf".format(trio), tool="merge")
    #
    # print("------------debreak_survivor")
    # for trio in trios:
    #     fetch_denovo_sv("/mnt/d/workspace/svision-pro/real_trio/supp10_supp2/debreak_survivor.{}.s10.s2.vcf".format(trio), tool="merge")
    # # #
    # collect_and_summary_denovo("/mnt/d/workspace/svision-pro/real_trio/supp10_supp2/trio_denovo/")


    # fetch_bam("/data/home/songbo/workspace/svision-pro/real_trio/sv_call_release/denovo_bams/collection.svision_pro.sniffles2.denovo.bed", "/data/home/songbo/workspace/svision-pro/real_trio/sv_call_release/denovo_bams")