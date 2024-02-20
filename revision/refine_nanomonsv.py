import os
import pysam

# # NanomonSV's VCF format does not fit the truvari comparison, so we need refine its format


# # STEP 1: bcftools sort

# # STEP 2: SVINSLEN --> SVLEN

# # STEP 3: reformat GT


raw_vcf_path = "/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/ont/nanomonsv.s2.vcf"

output_vcf_path = raw_vcf_path.replace(".vcf", ".tmp.vcf")


with open(output_vcf_path, "w") as fout, pysam.VariantFile(raw_vcf_path) as fin:

    fout.write("\n".join(str(fin.header).split("\n")[0: -2]))

    fout.write("\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n")
    fout.write(str(fin.header).split("\n")[-2] + "\n")

    for record in fin:

        start = record.start + 1

        end = str(record).split("\t")[7].split(";")[0]

        if "END=" in end and int(end[4: ]) < start:

            continue

        record_split = str(record).strip().split("\t")

        record_split[7] = record_split[7].replace("SVINSLEN=", "SVLEN=")

        record_split[-1] = "./.:" + record_split[-1]
        record_split[-2] = "./.:" + record_split[-2]
        record_split[-3] = "GT:" + record_split[-3]

        fout.write("\t".join(record_split) + "\n")


os.remove(raw_vcf_path)
os.system("mv {} {}".format(output_vcf_path, raw_vcf_path))
