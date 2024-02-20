import os
import pysam


def convert_vcf_to_bed(vcf_path):

    bed_path = vcf_path.replace(".vcf", ".bed")

    with open(bed_path, "w") as bed_out:
        for record in pysam.VariantFile(vcf_path):

            if "BND" in str(record):
                continue

            record_chrom = record.contig

            if "_" in record_chrom:
                continue

            if "jasmine" in vcf_path and "SUPP_VEC_EXT=1" not in str(record):
                continue

            if "survivor" in vcf_path and "SUPP_VEC=1" not in str(record):
                continue

            if "sniffles" in vcf_path and "SUPP_VEC=1" not in str(record):
                continue
            record_start = record.start

            record_length = record.info["SVLEN"]

            record_end = record_start + abs(record_length)

            bed_out.write("{}\t{}\t{}\t{}\t{}\n".format(record_chrom, record_start, record_end, record_length, record.info["SVTYPE"]))



if __name__ == '__main__':


    convert_vcf_to_bed("/mnt/d/workspace/svision-pro/paper_revision/pcr_quartet_ssv/release_lcl5-6-7-8_256.svision_pro_v1.6.s10.vcf")
    convert_vcf_to_bed("/mnt/d/workspace/svision-pro/paper_revision/pcr_quartet_ssv/sniffles2.lcl5-6-7-8.s10.vcf")
    convert_vcf_to_bed("/mnt/d/workspace/svision-pro/paper_revision/pcr_quartet_ssv/cutesv_jasmine.lcl5-6-7-8.s10.vcf")
    convert_vcf_to_bed("/mnt/d/workspace/svision-pro/paper_revision/pcr_quartet_ssv/cutesv_survivor.lcl5-6-7-8.s10.vcf")

    convert_vcf_to_bed("/mnt/d/workspace/svision-pro/paper_revision/pcr_quartet_ssv/debreak_jasmine.lcl5-6-7-8.s10.vcf")
    convert_vcf_to_bed("/mnt/d/workspace/svision-pro/paper_revision/pcr_quartet_ssv/debreak_survivor.lcl5-6-7-8.s10.vcf")

    convert_vcf_to_bed("/mnt/d/workspace/svision-pro/paper_revision/pcr_quartet_ssv/svision_jasmine.lcl5-6-7-8.s10.vcf")
    convert_vcf_to_bed("/mnt/d/workspace/svision-pro/paper_revision/pcr_quartet_ssv/svision_survivor.lcl5-6-7-8.s10.vcf")