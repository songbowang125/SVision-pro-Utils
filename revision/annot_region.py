import pandas as pd
import os
import pysam


def convert_vcf2bed_by_tools_raw_svtype(vcf_file, bed_out, tool):
    bed_writer = open(bed_out, "w")

    svtypes = {}

    vcf_base = pysam.VariantFile(vcf_file, 'r')

    for record in vcf_base.fetch():
        chrom = record.contig
        start = record.pos + 1
        end = record.stop + 2
        svtype = record.info['SVTYPE']

        if svtype not in svtypes.keys():
            svtypes[svtype] = 1
        else:
            svtypes[svtype] += 1

        if record.info['SVTYPE'] == 'INS':
            if tool == 'pbsv':
                size = record.info['SVLEN'][0]
            else:
                size = record.info['SVLEN']

            if size is None:
                continue

            # end = start + 1
            end = start + size

        else:
            size = abs(end - start) + 1

        bed_str = "{0}\t{1}\t{2}\t{3}\t{4}\n".format(chrom, start, end, size, svtype)
        bed_writer.write(bed_str)

    print(svtypes)
    bed_writer.close()


def convert_vcf2bed_by_tools(vcf_file, bed_out, tool):
    bed_writer = open(bed_out, "w")

    svtypes = {}

    vcf_base = pysam.VariantFile(vcf_file, 'r')

    for record in vcf_base.fetch():
        chrom = record.contig
        start = record.pos + 1
        end = record.stop + 2
        svtype = record.info['SVTYPE']

        if svtype not in svtypes.keys():
            svtypes[svtype] = 1
        else:
            svtypes[svtype] += 1

        if tool in ["survivor"]:
            svtype = "DEL"
        elif tool in ["sniffles", "svim", "pbsv", "cutesv", ]:
            if svtype == "BND":
                continue
            elif svtype == "DUP":
                svtype = 'tDUP'
            elif svtype == 'INV/INVDUP':
                svtype = 'itDUP'
            elif svtype == 'DUP/INS':
                svtype = 'dDUP'
            elif svtype == 'INVDUP' or svtype == 'DUP:INT':
                svtype = 'itDUP'
            elif svtype == 'DEL/INV':
                svtype = 'DEL+INV'
            elif svtype == "DUP:TANDEM":
                svtype = 'tDUP'
            elif svtype == 'DUP_INT':
                svtype = 'itDUP'
            elif svtype == 'DEL/INV/INVDUP':
                svtype = 'DEL+DUP+INV'
            else:
                if svtype not in ['INS', 'DEL', "INV", 'cnv']:
                    print("Not converted ", svtype)
        elif tool in ["svision"]:
            if svtype == "INS+tDUP":
                svtype = "tDUP"
            elif svtype == "INS+INV+tDUP":
                svtype = "itDUP"
            elif svtype == "INS+DUP":
                svtype = "dDUP"
            elif svtype == "INS+INV+DUP":
                svtype = "idDUP"
            elif svtype == "DEL+DUP":
                svtype = "DEL+dDUP"
            elif svtype == "DEL+INV+DUP":
                svtype = "DEL+idDUP"
            elif svtype == "DEL+tDUP":
                svtype = "DEL+tDUP"
            elif svtype == "DEL+INV+tDUP":
                svtype = "DEL+itDUP"


        elif tool in ['svision-pro']:
            pass
        elif tool in ["benchmark", "hg002_benchmark"]:
            pass
        else:
            print("error: no this tool")
            exit()

        if record.info['SVTYPE'] == 'INS':
            if tool == 'pbsv':
                size = record.info['SVLEN'][0]
            elif tool == "hg002_benchmark":
                size = abs(record.info['SVLEN'][0])

            elif tool == "benchmark":
                size = len(record.alts[0])
            else:
                size = record.info['SVLEN']

            if size is None:
                continue

            # end = start + 1
            end = start + size

        else:
            size = abs(end - start) + 1

        # # STEP: parse detail bkps
        if tool in ["svision"]:
            bkps_str = []
            for bkp in record.info["BKPS"]:
                bkp = bkp.replace(":", "_")
                bkp = bkp.replace("-", "_")

                bkps_str.append(bkp)

            bkps_str = ",".join(bkps_str)
            it_str = "."

        elif tool in ["svision-pro"]:
            bkps_str = []

            for bkp in record.info["BKPS"]:
                bkp_split = bkp.split("_")
                bkp = "{}_{}_{}".format(bkp_split[0], bkp_split[2], bkp_split[3])

                bkps_str.append(bkp)
            bkps_str = ",".join(bkps_str)

            try:
                it_str = record.info["BKPSIT"]
            except:
                it_str = "."
        elif tool in ['hg002_benchmark']:
            bkps_str = "{}_{}_{}".format(svtype, start, end)

            it_str = "."

        elif tool in ["benchmark"]:
            bkps_str = []

            for bkp in record.info["BKPS"]:
                bkp_split = bkp.replace(":", "_").split("_")
                bkp = "{}_{}_{}".format(bkp_split[0], bkp_split[2], bkp_split[3])

                bkps_str.append(bkp)
            bkps_str = ",".join(bkps_str)
            it_str = "."

        else:
            bkps_str = "NA"
            it_str = "."

        if tool in ["survivor"]:
            supp_vec = record.info["SUPP_VEC"]
        else:
            supp_vec = "."
        bed_str = "{0}\t{1}\t{2}\t{3}\t{4}\tsupp={5}\n".format(chrom, start, end, size, svtype, supp_vec)
        bed_writer.write(bed_str)

    print(svtypes)
    bed_writer.close()


def overlap(a,b,c,d):
    r = 0 if a == c and b == d else min(b, d)-max(a, c)
    if r >= 0:
        return r


def overlap_repeat_elements(in_bed, out_dir, rmsk_path, bedtools='bedtools'):

    bed_prefix = ".".join(os.path.basename(in_bed).split(".")[0:-1])

    overlap_out_file = "{}/{}.repeat_overlap.bed".format(out_dir, bed_prefix)

    rmsk_cmd = '{0} intersect -wa -wb -a {1} -b {2} > {3}'.format(bedtools, in_bed, rmsk_path, overlap_out_file)
    os.system(rmsk_cmd)

    return overlap_out_file


def assign_region_rmsk(in_bed, rmsk_overlaps, annot_out_path, all_possible_chrs):
    bed_prefix = ".".join(os.path.basename(in_bed).split(".")[0:-1])

    annot_out_origin = os.path.join(annot_out_path, "{}.repeat_annot.bed".format(bed_prefix))

    # df_rmsk = pd.read_csv(rmsk_overlaps, sep="\t", names=['chrom', 'start', 'end', 'svlen', 'type', 'rpchrom', 'rpstart', 'rpend', 'rpsubtype', 'rptype'])
    # df_rmsk = pd.read_csv(rmsk_overlaps, sep="\t", names=['chrom', 'start', 'end', 'svlen', 'type', 'rpchrom', 'rpstart', 'rpend', 'rpsubtype', 'rptype'])

    # df_rmsk = pd.read_csv(rmsk_overlaps, sep="\t", names=['chrom', 'start', 'end', 'id', 'type', 'support', 'length', "VaPoR_gs",	"VaPoR_vaf", "VaPoR_GT", "VaPoR_GQ", "VaPoR_Rec", 'rpchrom', 'rpstart', 'rpend', 'rpsubtype', 'rptype'])
    df_rmsk = pd.read_csv(rmsk_overlaps, sep="\t", names=['chrom', 'start', 'end', 'length', 'type', 'sample', 'tool', "flag", 'rpchrom', 'rpstart', 'rpend', 'rpsubtype', 'rptype'])

    annot_by_sv = {}

    for idx, row in df_rmsk.iterrows():
        sv_id = "{0}-{1}-{2}".format(row['chrom'], row['start'], row['end'])
        overlap_size = overlap(int(row['start']), int(row['end']), int(row['rpstart']), int(row['rpend']))
        if sv_id in annot_by_sv:
            annot_by_sv[sv_id].append((row['rpsubtype'], row['rptype'], overlap_size))
        else:
            annot_by_sv[sv_id] = [(row['rpsubtype'], row['rptype'], overlap_size)]

    df_svs = pd.read_csv(in_bed, sep="\t", names=['chrom', 'start', 'end', 'length', 'type', 'sample', 'tool', "flag"])

    sv_annots = list()
    count_by_rptype = {}
    for idx, row in df_svs.iterrows():
        sv_id = "{0}-{1}-{2}".format(row['chrom'], row['start'], row['end'])

        # this_chrom_exclude_.regions = exclude_region_dict[row['chrom']]
        # if is_filtered_sv_by_region(this_chrom_exclude_regions, int(row['start']), int(row['end'])):
        #     continue
        rptype = 'None'
        rpsubtype = 'None'
        rep_overlaps = 0
        if sv_id in annot_by_sv:
            annots = annot_by_sv[sv_id]
            if len(annots) == 1:
                # print(sv_id)

                rpsubtype = annots[0][0]
                rptype = annots[0][1]
                if rptype == 'Simple_repeat':
                    motif = annots[0][0][1:-2]
                    if len(motif) <= 7:
                        rpsubtype = 'STR'
                    else:
                        rpsubtype = 'VNTR'
                rep_overlaps = annots[0][2]
            else:
                sorted_annots_by_size = sorted(annots, key=lambda x: x[2], reverse=True)
                rpsubtype = sorted_annots_by_size[0][0]
                rptype = sorted_annots_by_size[0][1]

                if rptype == 'Simple_repeat':
                    motif = sorted_annots_by_size[0][0][1:-2]
                    if len(motif) <= 7:
                        rpsubtype = 'STR'
                    else:
                        rpsubtype = 'VNTR'
                rep_overlaps = sorted_annots_by_size[0][2]
        msk_pcrt = 100 * rep_overlaps / (int(row['end']) - int(row['start']) + 1)

        this_sv = row.tolist()
        this_sv.extend([round(msk_pcrt, 2), rpsubtype, rptype])

        sv_annots.append(this_sv)

        if rptype in count_by_rptype:
            count_by_rptype[rptype] += 1
        else:
            count_by_rptype[rptype] = 1

    df_sv_annots = pd.DataFrame(sv_annots, columns=['chrom', 'start', 'end', 'length', 'type', 'sample', 'tool', "flag", 'pcrt', 'rpsubtype', 'rptype'])

    sorter_index = dict(zip(all_possible_chrs, range(len(all_possible_chrs))))
    df_sv_annots['chrom_rank'] = df_sv_annots['chrom'].map(sorter_index)
    df_sv_annots.drop('chrom_rank', 1, inplace=True)

    # output origin annot
    df_sv_annots_origin = df_sv_annots.copy()
    # del df_sv_annots_origin['support_reads']
    df_sv_annots_origin.to_csv(annot_out_origin, index=False, header=False, sep="\t")

    # df_sv_annots['rpsubtype'][df_sv_annots.pcrt < min_overlap_ratio] = 'None'
    # df_sv_annots['rptype'][df_sv_annots.pcrt < min_overlap_ratio] = 'None'
    #
    # # rm by min_overlap_ratio
    # if rm_TE == True:
    #     df_sv_annots = df_sv_annots.drop(df_sv_annots[df_sv_annots.pcrt >= min_overlap_ratio].index)
    # else:
    #     df_sv_annots = df_sv_annots.drop(df_sv_annots[(df_sv_annots.pcrt >= min_overlap_ratio) & (df_sv_annots.rptype == "Simple_repeat")].index)
    #     df_sv_annots = df_sv_annots.drop(df_sv_annots[(df_sv_annots.pcrt >= min_overlap_ratio) & (df_sv_annots.rptype == "Satellite")].index)
    #
    # df_sv_annots.to_csv(annot_out, index=False, header=False, sep="\t")

    annot_out_ro50 = os.path.join(annot_out_path, "{}.repeat_annot.RO50.bed".format(bed_prefix))
    df_sv_annots_ro50 = df_sv_annots_origin[df_sv_annots_origin.pcrt < 50]
    df_sv_annots_ro50.to_csv(annot_out_ro50, index=False, header=False, sep="\t")


    annot_out_ro50_large = os.path.join(annot_out_path, "{}.repeat_annot.RO50_large.bed".format(bed_prefix))
    df_sv_annots_ro50_large = df_sv_annots_origin[df_sv_annots_origin.pcrt >= 50]
    df_sv_annots_ro50_large.to_csv(annot_out_ro50_large, index=False, header=False, sep="\t")

    return annot_out_origin, count_by_rptype



def annot_pipe(bed_file, out_path, rp_file="/mnt/f/ref/grch38/repeat/hg38.trf.convert.bed"):

    # all_possible_chrs = ["1", "2", "3", "4", "5", "6","7", "8", "9", "10","11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X",]
    all_possible_chrs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6","chr7", "chr8", "chr9", "chr10","chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX",]

    overlap_out_file = overlap_repeat_elements(bed_file, out_path, rp_file)
    #
    assign_region_rmsk(bed_file, overlap_out_file, out_path, all_possible_chrs)


if __name__ == '__main__':
    # overlap_repeat_elements("/mnt/c/workspace/test/real_ccs/svision-pro_bench_region/tp-base.bed", "/mnt/c/workspace/test/real_ccs/svision-pro_bench_region", "/mnt/d/data/ref/grch37/repeat/hg19.rmsk_trf.bed")


    # all_possible_chrs = ["1", "2", "3", "4", "5", "6","7", "8", "9", "10","11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X",]
    #
    # assign_region_rmsk("/mnt/c/workspace/test/real_ccs/svision-pro_bench_region/tp-base.bed", "/mnt/c/workspace/test/real_ccs/svision-pro_bench_region/tp-base.repeat_overlap.bed" , "/mnt/c/workspace/test/real_ccs/svision-pro_bench_region", all_possible_chrs)

    # out_path = "/mnt/c/workspace/test_new/real_trio/twin_discordant"
    #
    # vcf_file = "/mnt/c/workspace/test_new/real_trio/twin_discordant/release_LCL5-6-7-8_256.svision_pro_v1.1.s10.denovo_pass.vcf"
    # bed_file = "/mnt/c/workspace/test_new/real_trio/twin_discordant/release_LCL5-6-7-8_256.svision_pro_v1.1.s10.denovo_pass.bed"
    # convert_vcf2bed_by_tools_raw_svtype(vcf_file, bed_file, "svision-pro")

    # vcf_file = "/mnt/c/workspace/test_new/real_quartet/sniffles2.LCL5-6-7-8.s10.vcf"
    # bed_file = "/mnt/c/workspace/test_new/real_quartet/sniffles2.LCL5-6-7-8.s10.bed"
    # convert_vcf2bed_by_tools_raw_svtype(vcf_file, bed_file, "sniffles")


    # vcf_file = "/mnt/c/workspace/test_new/real_quartet/cutesv_jasmine.LCL5-6-7-8.s10.vcf"
    # bed_file = "/mnt/c/workspace/test_new/real_quartet/cutesv_jasmine.LCL5-6-7-8.s10.bed"
    # convert_vcf2bed_by_tools_raw_svtype(vcf_file, bed_file, "jasmine")


    # vcf_file = "/mnt/c/workspace/test_new/real_quartet/twin_discordant/svision_jasmine.LCL5-6-7-8.s10.twin_discordant.vcf"
    # bed_file = "/mnt/c/workspace/test_new/real_quartet/twin_discordant/svision_jasmine.LCL5-6-7-8.s10.twin_discordant.bed"
    # convert_vcf2bed_by_tools_raw_svtype(vcf_file, bed_file, "jasmine")
    #
    # vcf_file = "/mnt/c/workspace/test_new/real_quartet/twin_discordant/debreak_jasmine.LCL5-6-7-8.s10.twin_discordant.vcf"
    # bed_file = "/mnt/c/workspace/test_new/real_quartet/twin_discordant/debreak_jasmine.LCL5-6-7-8.s10.twin_discordant.bed"
    # convert_vcf2bed_by_tools_raw_svtype(vcf_file, bed_file, "jasmine")

    annot_pipe( "/mnt/d/workspace/svision-pro/paper_revision/pcr_quartet_ssv2/collection.all.denovo.lcl5-6-7-8.highconf.bed", "/mnt/d/workspace/svision-pro/paper_revision/pcr_quartet_ssv2")

    # out_path = "/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic/modified/repeat_annot_80"
    #
    # for sample in ["tumor", "normal"]:
    #     for seq in ["hifi", "ont", "clr"]:
    #         bed_file = "/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic/modified/release_withbed.svision_pro_v1.7.s2.vcf.somatic.s2.{}.{}.txt".format(sample, seq)
    #         annot_pipe(bed_file, out_path)
    #
    #         bed_file = "/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic/modified/sniffles2.merge.s2.vcf.somatic.s2.{}.{}.txt".format(sample, seq)
    #         annot_pipe(bed_file, out_path)

    # vcf_file = "/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic/release.svision_pro_v1.7.s2.vcf.somatic.s2.vcf"
    # bed_file = "/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic/release.svision_pro_v1.7.s2.vcf.somatic.s2.bed"
    # convert_vcf2bed_by_tools_raw_svtype(vcf_file, bed_file, "svision-pro")
    #
    # out_path = "/mnt/d/workspace/svision-pro/paper_revision/real_nt_hcc1395/hifi/somatic"
    #
    # annot_pipe(bed_file, out_path)