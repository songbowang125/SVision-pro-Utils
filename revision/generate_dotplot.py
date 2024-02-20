import os
import pysam


def generate_plot():
    cmd = "java -cp ./utils/gepard/Gepard_altered.jar org.gepard.client.cmdline.CommandLine -seq {} {} -matrix ./utils/gepard/edna.mat -outfile {} -silent -word {} -zoom {}"


def fetch_ref_seq(ref_file, chrom, start, end):

    return ref_file.fetch(chrom, start, end + 1)


def plot_from_bed(bed_path, bam_path, ref_path, out_path, max_read_num=10):

    with open(bed_path) as bed_input, pysam.AlignmentFile(bam_path) as bam_input, pysam.FastaFile(ref_path) as ref_input:

        for line in bed_input:

            line_split = line.strip().split("\t")

            chrom = line_split[0]
            start = int(line_split[1])
            end = int(line_split[2])

            ref_seq = fetch_ref_seq(ref_input, chrom, start - 5000, end + 5000)
            ref_seq_out_path = os.path.join(out_path, "{}_{}_{}.ref.fa".format(chrom, start, end))
            with open(ref_seq_out_path, "w") as fout:
                fout.write(">ref\n")
                fout.write(ref_seq + "\n")

            read_num = 0
            for primary_align in bam_input.fetch(chrom, start - 100, end + 100):

                if primary_align.is_supplementary or primary_align.is_secondary:
                    continue

                primary_align_seq = primary_align.query_sequence

                read_seq_out_path = os.path.join(out_path, "{}_{}_{}_{}.read.fa".format(chrom, start, end, read_num))

                with open(read_seq_out_path, "w") as fout:
                    fout.write(">read\n")
                    fout.write(primary_align_seq + "\n")

                plot_out_path = os.path.join(out_path, "{}_{}_{}_{}.png".format(chrom, start, end, read_num))
                os.system("java -cp /mnt/d/tools/Gepard/Gepard-2.1.jar org.gepard.client.cmdline.CommandLine -seq {} {} -matrix /mnt/d/tools/Gepard/edna.mat -outfile {} -silent".format(ref_seq_out_path, read_seq_out_path, plot_out_path))

                os.remove(read_seq_out_path)
                read_num += 1

                if read_num > max_read_num:
                    break

            os.remove(ref_seq_out_path)


if __name__ == '__main__':

    bed_path = "/mnt/d/workspace/svision-pro/paper_revision/pcr_quartet_ssv2/collection.all.denovo.lcl5-6-7-8.highconf.repeat_annot.RO50.bed"
    bam_path = "/mnt/h/data/chinese_quartet/ChineseQuartet.LCL5.GRCh38.HiFi.minimap2.bam"
    # bam_path = "/mnt/f/hg002/GIAB.HG003.GRCh38.HiFi.minimap2.bam"


    ref_path = "/mnt/f/ref/grch38/GRCh38.d1.vd1.fa"

    out_path = "/mnt/d/workspace/svision-pro/paper_revision/pcr_quartet_ssv2/dotplot/"

    plot_from_bed(bed_path, bam_path, ref_path, out_path, max_read_num=10)