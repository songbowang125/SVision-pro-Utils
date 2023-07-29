
import os
import shutil


def compare(vcf_file, tool_out_path, compare_type):

    print(vcf_file)

    gz_vcf = vcf_file + '.gz'

    cmd_str = 'bgzip -c {0} > {1}'.format(vcf_file, gz_vcf)
    os.system(cmd_str)

    cmd_str = 'tabix -p vcf {0}'.format(gz_vcf)
    os.system(cmd_str)

    if compare_type == 'region':

        if base_bed is None:
            if not pass_only:
                cmd_str = 'truvari bench -b {} -c {} --typeignore -r 1000 -p 0 -f {} -o {}'.format(base_gz, gz_vcf, ref, tool_out_path)
            else:
                cmd_str = 'truvari bench -b {} -c {} --typeignore -r 1000 -p 0 -f {} -o {} --passonly'.format(base_gz, gz_vcf, ref, tool_out_path)
        else:
            if not pass_only:
                cmd_str = 'truvari bench -b {} -c {} --includebed {} --typeignore  -r 1000 -p 0 -f {} -o {}'.format(base_gz, gz_vcf, base_bed, ref, tool_out_path)
            else:
                cmd_str = 'truvari bench -b {} -c {} --includebed {} --typeignore  -r 1000 -p 0 -f {} -o {} --passonly'.format(base_gz, gz_vcf, base_bed, ref, tool_out_path)

        os.system(cmd_str)

    elif compare_type == 'type':
        if base_bed is None:
            if not pass_only:
                cmd_str = 'truvari bench -b {} -c {}  -r 1000 -p 0 -f {} -o {}'.format(base_gz, gz_vcf, ref, tool_out_path)
            else:
                cmd_str = 'truvari bench -b {} -c {}  -r 1000 -p 0 -f {} -o {} --passonly'.format(base_gz, gz_vcf, ref, tool_out_path)
        else:
            if not pass_only:
                cmd_str = 'truvari bench -b {} -c {} --includebed {} -r 1000 -p 0 -f {} -o {}'.format(base_gz, gz_vcf, base_bed, ref, tool_out_path)
            else:
                cmd_str = 'truvari bench -b {} -c {} --includebed {} -r 1000 -p 0 -f {} -o {} --passonly'.format(base_gz, gz_vcf, base_bed, ref, tool_out_path)

        os.system(cmd_str)
    else:
        pass


class real_ssv_params:
    def __init__(self):
        self.work_path = '/mnt/c/workspace/test_new/real_ccs'
        self.ref = '/mnt/d/data/ref/grch37/hs37d5.fa'
        self.base_gz = '/mnt/c/workspace/test_new/real_ccs/ground_truth/HG002_SVs_Tier1_v0.6.vcf.gz'
        # self.base_gz = '/mnt/c/workspace/test_new/sim_csv_ccs/ground_truth/extended/benchmark_origin_clone0_300.csv.vcf.gz'

        self.base_bed = '/mnt/c/workspace/test_new/real_ccs/ground_truth/HG002_SVs_Tier1_v0.6.bed'
        self.pass_only = True
        self.type = "region"

        self.tools = ["release.svision_pro_v1.6.s10"]
        # self.tools = ["HG002.cutesv.s10.m50", "HG002.sniffles2.s10.m50", "HG002.svision.s10", "HG002.pbsv.s10.m50", "HG002.debreak.s10.m50"]
        # self.tools = ["HG002.debreak.s10.m50"]
        # self.tools = ["HG002-ONT.cutesv.s10.m50", "HG002-ONT.sniffles2.s10.m50", "HG002-ONT.svision.s10", "HG002-ONT.pbsv.s10.m50", "HG002-ONT.debreak.s10.m50"]


class sim_csv_params:
    def __init__(self):
        self.work_path = '/mnt/c/workspace/test_new/sim_csv_ccs/'
        self.ref = "/mnt/d/data/ref/grch38/GRCh38.d1.vd1.fa"
        self.base_gz = '/mnt/c/workspace/test_new/sim_csv_ccs/ground_truth/benchmark_origin_clone0_300.csv.vcf.gz'
        # self.base_gz = '/mnt/c/workspace/test_new/sim_csv_ccs/ground_truth/extended/benchmark_origin_clone0_300.csv.vcf.gz'

        self.base_bed = None
        self.pass_only = False
        self.type = "region"

        # self.tools = [ "sim_csv.cutesv.s5.m50", "sim_csv.debreak.s5.m50", "sim_csv.sniffles2.s5.m50", "sim_csv.svision.s5"]
        self.tools = ["release.svision_pro_v1.6.s5", "release.svision_pro_v1.6.s10"]

        # self.tools = [ "sim_csv.cutesv.s10.m50", "sim_csv.debreak.s10.m50", "sim_csv.sniffles2.s10.m50", "sim_csv.svision.s10"]
        # self.tools = ["sim_csv.svision.s5", "sim_csv.svision.s10"]


class sim_denovo_csv_params:
    def __init__(self):
        self.work_path = '/mnt/c/workspace/test_new/sim_denovo_ccs/complex/'

        self.ref = "/mnt/d/data/ref/grch38/GRCh38.d1.vd1.fa"
        # self.base_gz = '/mnt/c/workspace/test_new/sim_denovo_ccs/ground_truth/extended/benchmark_origin_clone0_300.csv.addTrio.vcf.gz'

        self.base_gz = '/mnt/c/workspace/test_new/sim_denovo_ccs/ground_truth/benchmark_origin_clone0_300.csv.addTrio.vcf.gz'

        self.base_bed = None
        self.pass_only = False
        self.type = "region"
        # self.tools = ["cutesv_jasmine.child-father-mother.s10", "cutesv_survivor.child-father-mother.s10", "debreak_jasmine.child-father-mother.s10", "debreak_survivor.child-father-mother.s10", "svision_jasmine.child-father-mother.s10", "svision_survivor.child-father-mother.s10", "sniffles2.child-father-mother.s10"]
        #
        self.tools = ["release_256.svision_pro_v1.6.s10"]


class sim_denovo_ssv_params:
    def __init__(self):
        self.work_path = '/mnt/c/workspace/test_new/sim_denovo_ccs/simple/'
        self.ref = '/mnt/d/data/ref/grch37/hs37d5.fa'

        self.base_gz = '/mnt/c/workspace/test_new/sim_denovo_ccs/ground_truth/HG002_SVs_Tier1_v0.6.final.addTrio.vcf.gz'

        self.base_bed = None
        self.pass_only = False
        self.type = "region"
        self.tools = ["cutesv.s10.jasmine", "cutesv.s10.survivor", "debreak.s10.jasmine", "debreak.s10.survivor", "svision.s10.jasmine", "svision.s10.survivor", "sniffles2.s10"]


class real_colo829:
    def __init__(self):
        self.work_path = '/mnt/c/workspace/test_new/real_nt_colo829/'
        self.ref = '/mnt/d/data/ref/grch37/hs37d5.fa'
        self.base_gz = '/mnt/c/workspace/test_new/real_nt_colo829/truthset_somaticSVs_COLO829.vcf.gz'

        self.base_bed = None
        self.pass_only = False

        self.type = "region"

        # self.tools = ["release_rxivdata.svision_pro_v1.0.s3"]

        self.tools = ["sniffles2.non-germline.s3"]


class sim_somatic_complex:
    def __init__(self):
        self.work_path = '/mnt/c/workspace/test_new/sim_somatic_ccs/complex'
        self.ref = "/mnt/d/data/ref/grch38/GRCh38.d1.vd1.fa"
        self.base_gz = '/mnt/c/workspace/test_new/sim_somatic_ccs/complex/benchmark_origin_clone0_300.csv.addVAF.vcf.gz'

        self.base_bed = None
        self.pass_only = False

        self.type = "region"

        self.tools = [ "release_1024.svision_pro_v1.6.s1",  "release_512.svision_pro_v1.6.s1",  "release_256.svision_pro_v1.6.s1"]
        # self.tools = [ "release_1024.svision_pro_v1.6.s1"]


class sim_somatic_simple:
    def __init__(self):
        self.work_path = '/mnt/c/workspace/test_new/sim_somatic_ccs/simple/'
        self.ref = '/mnt/d/data/ref/grch37/hs37d5.fa'

        self.base_gz = '/mnt/c/workspace/test_new/sim_somatic_ccs/simple/HG002_SVs_Tier1_v0.6_high_confident.addVAF.vcf.gz'

        self.base_bed = None
        self.pass_only = False

        self.type = "region"

        self.tools = [ "release_1024.svision_pro_v1.6.s1",  "release_512.svision_pro_v1.6.s1",  "release_256.svision_pro_v1.6.s1"]
        # self.tools = [ "release_1024.svision_pro_v1.6.s1"]

        # self.tools = ["sniffles2.non-germline.s1", "nanomonsv.s1"]
        # self.tools = ["nanomonsv.s1"]


class real_hcc1395:
    def __init__(self):
        self.work_path = '/mnt/c/workspace/test_new/real_nt_hcc1395'
        self.ref = "/mnt/d/data/ref/grch38/GRCh38.d1.vd1.fa"
        self.base_gz = '/mnt/c/workspace/test_new/real_nt_hcc1395/GB_1788_SVs.vcf.gz'

        self.base_bed = None
        self.pass_only = False

        self.type = "region"

        # self.tools = ["release.svision_pro_v1.6.s2"]
        # self.tools = ["sniffles.non_germline.s2"]
        self.tools = ["nanomonsv.s2"]

if __name__ == '__main__':

    params = real_hcc1395()

    work_path = params.work_path
    ref = params.ref
    base_gz = params.base_gz
    base_bed = params.base_bed
    pass_only = params.pass_only
    compare_type = params.type

    tools = params.tools


    for tool in tools:
        # # STEP: run
        tool_vcf = os.path.join(work_path, '{0}.vcf'.format(tool))

        sorted_vcf = os.path.join(work_path, '{0}.sorted.vcf'.format(tool))

        cmd_str = 'bcftools sort -o {0} {1}'.format(sorted_vcf, tool_vcf)
        os.system(cmd_str)

        tool_out_path = os.path.join(work_path, tool + '_bench_{0}'.format(compare_type))
        if os.path.exists(tool_out_path):
            shutil.rmtree(tool_out_path)

        compare(sorted_vcf, tool_out_path, compare_type)

        os.system("rm {}*".format(sorted_vcf))


