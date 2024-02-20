import os
import pysam


def parse_ttmars_res(work_dir):

    # samples = ["HG00514", "HG00512", "HG00513", "HG00733", "HG00731", "HG00732", "NA19240", "NA19239", "NA19238"]
    samples = ["HG00514", "HG00512", "HG00513"]

    ttmars_res = {}

    vcf_file = pysam.VariantFile(os.path.join(work_dir, "denovo.sv.vcf"))

    for record in vcf_file:

        record_chrom = record.contig
        record_start = record.start + 1
        record_end = record.stop

        record_sv = record.info["SVTYPE"]
        record_length = record.info["SVLEN"]

        record_sample = record.id

        sv_id = "{}_{}_{}_{}".format(record_chrom, record_start, record_end, record_sv)

        ttmars_res[sv_id] = {"HG00514": "NA", "HG00512": "NA", "HG00513": "NA", "HG00733": "NA", "HG00731": "NA", "HG00732": "NA", "NA19240": "NA", "NA19238": "NA", "NA19239": "NA"}

    vcf_file.close()

    # # load each file

    for sample in samples:

        sample_file = open(os.path.join(work_dir, "ttmars_combined_res.{}.txt".format(sample)))

        for line in sample_file:
            line_split = line.strip().split("\t")

            sv_id = "_".join(line_split[-4: ])

            ttmars_length_score = float(line_split[1])

            ttmars_valid = line_split[3]

            if abs(ttmars_length_score - 1) > 0.1:
                ttmars_valid = "False"

            ttmars_res[sv_id][sample] = ttmars_valid

    # # output
    print("{}\t{}\t{}\t{}\t{}\tHG00514 Trio\tHG00733 Trio\tNA19240 Trio".format("#Chrom", "Start", "End", "Type", "\t".join(samples)))

    for sv_id in ttmars_res:

        sample_valids = [ttmars_res[sv_id][sample] for sample in samples]

        hg00514_trio_res = "_".join([ttmars_res[sv_id][sample] for sample in ["HG00514", "HG00512", "HG00513"]])
        hg00733_trio_res = "_".join([ttmars_res[sv_id][sample] for sample in ["HG00733", "HG00731", "HG00732"]])
        na19240_trio_res = "_".join([ttmars_res[sv_id][sample] for sample in ["NA19240", "NA19239", "NA19238"]])

        if "NA" in [ttmars_res[sv_id][sample] for sample in ["HG00514", "HG00512", "HG00513"]]:
            hg00514_trio_flag = "NA"
        elif ttmars_res[sv_id]["HG00514"] == "True" and ttmars_res[sv_id]["HG00512"] == "False" and ttmars_res[sv_id]["HG00513"] == "False":
            hg00514_trio_flag = "De novo"
        elif ttmars_res[sv_id]["HG00514"] == "True" and "True" in [ttmars_res[sv_id][sample] for sample in ["HG00512", "HG00513"]]:
            hg00514_trio_flag = "Inherited"
        else:
            hg00514_trio_flag = "Bad"

        if "NA" in [ttmars_res[sv_id][sample] for sample in ["HG00733", "HG00731", "HG00732"]]:
            hg00733_trio_flag = "NA"
        elif ttmars_res[sv_id]["HG00733"] == "True" and ttmars_res[sv_id]["HG00731"] == "False" and ttmars_res[sv_id]["HG00732"] == "False":
            hg00733_trio_flag = "De novo"
        elif ttmars_res[sv_id]["HG00733"] == "True" and "True" in [ttmars_res[sv_id][sample] for sample in ["HG00731", "HG00732"]]:
            hg00733_trio_flag = "Inherited"
        else:
            hg00733_trio_flag = "Bad"

        if "NA" in [ttmars_res[sv_id][sample] for sample in ["NA19240", "NA19239", "NA19238"]]:
            na19240_trio_flag = "NA"
        elif ttmars_res[sv_id]["NA19240"] == "True" and ttmars_res[sv_id]["NA19239"] == "False" and ttmars_res[sv_id]["NA19238"] == "False":
            na19240_trio_flag = "De novo"
        elif ttmars_res[sv_id]["NA19240"] == "True" and "True" in [ttmars_res[sv_id][sample] for sample in ["NA19239", "NA19238"]]:
            na19240_trio_flag = "Inherited"
        else:
            na19240_trio_flag = "Bad"

        sv_id_split = sv_id.split("_")

        print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(sv_id_split[0], sv_id_split[1], sv_id_split[2], sv_id_split[3], "\t".join(sample_valids), hg00514_trio_flag, hg00733_trio_flag, na19240_trio_flag))


if __name__ == '__main__':
    parse_ttmars_res("/mnt/d/workspace/svision-pro/paper_revision/denovo_ttmars/hg00514")