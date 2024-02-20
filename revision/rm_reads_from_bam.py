import os
import pysam

input_bam = "/data/DATA/CellLine/HCC1395/HiFi_2023Pacbio/HCC1395.GRCh38.bam"

# input_bam = "/data/DATA/CellLine/HCC1395/HiFi_2023Pacbio/HCC1395-BL.GRCh38.bam"


output_bam = "/data/DATA/CellLine/HCC1395/HiFi_2023Pacbio/HCC1395.GRCh38.nanomonsv.bam"

# output_bam = "/data/DATA/CellLine/HCC1395/HiFi_2023Pacbio/HCC1395-BL.GRCh38.nanomonsv.bam"



with pysam.AlignmentFile(input_bam) as fin:

    with pysam.AlignmentFile(output_bam, "wb", header=fin.header) as fout:


        for align in fin:

            if align.qname not in ["m84039_230414_235240_s2/249105983/ccs", "m84039_230415_005427_s4/264705331/ccs"]:

                fout.write(align)

            else:
                print(align.qname)
