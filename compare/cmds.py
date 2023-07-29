

"""

# # # # # # STEP: germline
python SVision-pro.py --target_path /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/GIAB/HG002/GIAB.HG002.GRCh37.HiFi.minimap2.bam --access_path ./src/pre_process/hg19.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_256_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/real_ccs/sv_call_release/ --genome_path /data/DATA/Human_genomes/reference_based_analysis/ref/GRCh37.fasta --sample_name release --min_supp 10 --preset hifi --process_num 48
python SVision-pro.py --target_path /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/GIAB/HG002/GIAB.HG002.GRCh37.ONTUL.minimap2.bam --access_path ./src/pre_process/hg19.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_256_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/real_ont/sv_call_release/ --genome_path /data/DATA/Human_genomes/reference_based_analysis/ref/GRCh37.fasta --sample_name release --min_supp 10 --preset error-prone --process_num 48

python SVision-pro.py --target_path /data/home/songbo/workspace/svision-pro/sim_csv_ccs/sim_csv.srt.bam --access_path ./src/pre_process/hg38.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_256_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/sim_csv_ccs/sv_call_release/ --genome_path ~/workspace/ref/chr1-X.fa --sample_name release --min_supp 10 --preset hifi --process_num 48 --force_cluster
python SVision-pro.py --target_path /data/home/songbo/workspace/svision-pro/sim_csv_ont/sim_csv.srt.bam --access_path ./src/pre_process/hg38.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_256_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/sim_csv_ont/sv_call_release/ --genome_path ~/workspace/ref/chr1-X.fa --sample_name release --min_supp 10 --preset error-prone --process_num 48 --force_cluster


# # # # # # STEP: simulate denovo csv
python SVision-pro.py --target_path /data/home/songbo/workspace/svision-pro/sim_denovo_ccs/complex/child.srt.bam --base_path  /data/home/songbo/workspace/svision-pro/sim_denovo_ccs/complex/father.srt.bam /data/home/songbo/workspace/svision-pro/sim_denovo_ccs/complex/mother.srt.bam --access_path ./src/pre_process/hg38.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_256_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/sim_denovo_ccs/complex/sv_call_release --genome_path ~/workspace/ref/chr1-X.fa --sample_name release_256 --min_supp 10 --preset hifi --process_num 48 --detect_mode denovo --force_cluster
python SVision-pro.py --target_path /data/home/songbo/workspace/svision-pro/sim_denovo_ont/complex/child.srt.bam --base_path  /data/home/songbo/workspace/svision-pro/sim_denovo_ont/complex/father.srt.bam /data/home/songbo/workspace/svision-pro/sim_denovo_ont/complex/mother.srt.bam --access_path ./src/pre_process/hg38.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_256_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/sim_denovo_ont/complex/sv_call_release --genome_path ~/workspace/ref/chr1-X.fa --sample_name release_256 --min_supp 10 --preset error-prone --process_num 48 --detect_mode denovo --force_cluster


# # # # # # STEP: simulate somatic csv

python SVision-pro.py --target_path ~/workspace/svision-pro/sim_somatic_ccs/complex/merged.srt.bam --base_path ~/workspace/svision-pro/sim_somatic_ccs/complex/simulated.100x.srt.bam --access_path ./src/pre_process/hg38.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_256_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/sim_somatic_ccs/complex/sv_call_release/ --genome_path ~/workspace/ref/chr1-X.fa --sample_name release_256 --min_supp 1 --preset hifi --process_num 48 --detect_mode somatic --img_size 256 --min_mapq 0 --force_cluster
python SVision-pro.py --target_path ~/workspace/svision-pro/sim_somatic_ccs/complex/merged.srt.bam --base_path ~/workspace/svision-pro/sim_somatic_ccs/complex/simulated.100x.srt.bam --access_path ./src/pre_process/hg38.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_512_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/sim_somatic_ccs/complex/sv_call_release/ --genome_path ~/workspace/ref/chr1-X.fa --sample_name release_512 --min_supp 1 --preset hifi --process_num 48 --detect_mode somatic --img_size 512 --min_mapq 0 --force_cluster
python SVision-pro.py --target_path ~/workspace/svision-pro/sim_somatic_ccs/complex/merged.srt.bam --base_path ~/workspace/svision-pro/sim_somatic_ccs/complex/simulated.100x.srt.bam --access_path ./src/pre_process/hg38.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_1024_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/sim_somatic_ccs/complex/sv_call_release/ --genome_path ~/workspace/ref/chr1-X.fa --sample_name release_1024 --min_supp 1 --preset hifi --process_num 48 --detect_mode somatic --img_size 1024 --min_mapq 0 --force_cluster

python SVision-pro.py --target_path ~/workspace/svision-pro/sim_somatic_ont/complex/merged.srt.bam --base_path ~/workspace/svision-pro/sim_somatic_ont/complex/simulated.100x.srt.bam --access_path ./src/pre_process/hg38.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_256_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/sim_somatic_ont/complex/sv_call_release/ --genome_path ~/workspace/ref/chr1-X.fa --sample_name release_256 --min_supp 1 --preset error-prone --process_num 48 --detect_mode somatic --img_size 256 --min_mapq 0 --force_cluster
python SVision-pro.py --target_path ~/workspace/svision-pro/sim_somatic_ont/complex/merged.srt.bam --base_path ~/workspace/svision-pro/sim_somatic_ont/complex/simulated.100x.srt.bam --access_path ./src/pre_process/hg38.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_512_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/sim_somatic_ont/complex/sv_call_release/ --genome_path ~/workspace/ref/chr1-X.fa --sample_name release_512 --min_supp 1 --preset error-prone --process_num 48 --detect_mode somatic --img_size 512 --min_mapq 0 --force_cluster
python SVision-pro.py --target_path ~/workspace/svision-pro/sim_somatic_ont/complex/merged.srt.bam --base_path ~/workspace/svision-pro/sim_somatic_ont/complex/simulated.100x.srt.bam --access_path ./src/pre_process/hg38.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_1024_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/sim_somatic_ont/complex/sv_call_release/ --genome_path ~/workspace/ref/chr1-X.fa --sample_name release_1024 --min_supp 1 --preset error-prone --process_num 48 --detect_mode somatic --img_size 1024 --min_mapq 0 --force_cluster


# # # # # # STEP: simulate somatic ssv

python SVision-pro.py --target_path ~/workspace/svision-pro/sim_somatic_ccs/simple/merged.srt.bam --base_path ~/workspace/svision-pro/sim_somatic_ccs/simple/simulated.100x.srt.bam --access_path ./src/pre_process/hg19.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_256_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/sim_somatic_ccs/simple/sv_call_release/ --genome_path ~/workspace/ref/hs37d5.1-Y.fa --sample_name release_256 --min_supp 1 --preset hifi --process_num 48 --detect_mode somatic --img_size 256 --min_mapq 0
python SVision-pro.py --target_path ~/workspace/svision-pro/sim_somatic_ccs/simple/merged.srt.bam --base_path ~/workspace/svision-pro/sim_somatic_ccs/simple/simulated.100x.srt.bam --access_path ./src/pre_process/hg19.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_512_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/sim_somatic_ccs/simple/sv_call_release/ --genome_path ~/workspace/ref/hs37d5.1-Y.fa --sample_name release_512 --min_supp 1 --preset hifi --process_num 48 --detect_mode somatic --img_size 512 --min_mapq 0
python SVision-pro.py --target_path ~/workspace/svision-pro/sim_somatic_ccs/simple/merged.srt.bam --base_path ~/workspace/svision-pro/sim_somatic_ccs/simple/simulated.100x.srt.bam --access_path ./src/pre_process/hg19.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_1024_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/sim_somatic_ccs/simple/sv_call_release/ --genome_path ~/workspace/ref/hs37d5.1-Y.fa --sample_name release_1024 --min_supp 1 --preset hifi --process_num 48 --detect_mode somatic --img_size 1024 --min_mapq 0

python SVision-pro.py --target_path ~/workspace/svision-pro/sim_somatic_ont/simple/merged.srt.bam --base_path ~/workspace/svision-pro/sim_somatic_ont/simple/simulated.100x.srt.bam --access_path ./src/pre_process/hg19.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_256_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/sim_somatic_ont/simple/sv_call_release/ --genome_path ~/workspace/ref/hs37d5.1-Y.fa --sample_name release_256 --min_supp 1 --preset error-prone --process_num 48 --detect_mode somatic --img_size 256 --min_mapq 0
python SVision-pro.py --target_path ~/workspace/svision-pro/sim_somatic_ont/simple/merged.srt.bam --base_path ~/workspace/svision-pro/sim_somatic_ont/simple/simulated.100x.srt.bam --access_path ./src/pre_process/hg19.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_512_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/sim_somatic_ont/simple/sv_call_release/ --genome_path ~/workspace/ref/hs37d5.1-Y.fa --sample_name release_512 --min_supp 1 --preset error-prone --process_num 48 --detect_mode somatic --img_size 512 --min_mapq 0
python SVision-pro.py --target_path ~/workspace/svision-pro/sim_somatic_ont/simple/merged.srt.bam --base_path ~/workspace/svision-pro/sim_somatic_ont/simple/simulated.100x.srt.bam --access_path ./src/pre_process/hg19.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_1024_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/sim_somatic_ont/simple/sv_call_release/ --genome_path ~/workspace/ref/hs37d5.1-Y.fa --sample_name release_1024 --min_supp 1 --preset error-prone --process_num 48 --detect_mode somatic --img_size 1024 --min_mapq 0



sniffles --input simulated.100x.srt.bamsurgeon.round1.srt.bam --vcf sv_call_round1_100x/sniffles2.non-germline.s1.vcf --minsupport 1 --minsupport-auto-mult 0.01 --minsvlen 50 --mapq 20 --non-germline -t 48

nanomonsv parse $tumor $tumor_prefix
nanomonsv parse $normal $normal_prefix
nanomonsv get --min_tumor_variant_read_num 1 --min_tumor_VAF 0.01 --processes $thread --control_prefix $normal_prefix --control_bam $normal $tumor_prefix $tumor $REF


# # # # # # STEP: normal-tumor paired cell line

python SVision-pro.py --target_path /data/DATA/CellLine/HCC1395/HCC1395.PB.minimap2.srt.bam --base_path /data/DATA/CellLine/HCC1395/HCC1395BL.PB.minimap2.srt.bam --access_path ./src/pre_process/hg38.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_1024_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/real_hcc1395/sv_call_release/ --genome_path ~/workspace/ref/chr1-X.fa --sample_name release --min_supp 2 --preset error-prone --process_num 48 --detect_mode somatic --img_size 1024 --min_mapq 0 --force_cluster

# # # # # # STEP: SVision-pro multisample genotype

python SVision-pro.py --target_path /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/ChineseQuartet/LCL7/ChineseQuartet.LCL7.GRCh38.HiFi.minimap2.bam --base_path /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/ChineseQuartet/LCL6/ChineseQuartet.LCL6.GRCh38.HiFi.minimap2.bam /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/ChineseQuartet/LCL5/ChineseQuartet.LCL5.GRCh38.HiFi.minimap2.bam /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/ChineseQuartet/LCL8/ChineseQuartet.LCL8.GRCh38.HiFi.minimap2.bam /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/HGSVC/NA19238/HGSVC.NA19238.GRCh38.HiFi.minimap2.bam /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/HGSVC/NA19239/HGSVC.NA19239.GRCh38.HiFi.minimap2.bam /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/HGSVC/NA19240/HGSVC.NA19240.GRCh38.HiFi.minimap2.bam /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/HGSVC/HG00512/HGSVC.HG00512.GRCh38.HiFi.minimap2.bam /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/HGSVC/HG00513/HGSVC.HG00513.GRCh38.HiFi.minimap2.bam /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/HGSVC/HG00514/HGSVC.HG00514.GRCh38.HiFi.minimap2.bam /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/HGSVC/HG00731/HGSVC.HG00731.GRCh38.HiFi.minimap2.bam /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/HGSVC/HG00732/HGSVC.HG00732.GRCh38.HiFi.minimap2.bam /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/HGSVC/HG00733/HGSVC.HG00733.GRCh38.HiFi.minimap2.bam /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/GIAB/HG002/GIAB.HG002.GRCh38.HiFi.minimap2.bam /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/GIAB/HG003/GIAB.HG003.GRCh38.HiFi.minimap2.bam /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/GIAB/HG004/GIAB.HG004.GRCh38.HiFi.minimap2.bam /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/GIAB/HG005/GIAB.HG005.GRCh38.HiFi.minimap2.bam /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/GIAB/HG006/GIAB.HG006.GRCh38.HiFi.minimap2.bam /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/GIAB/HG007/GIAB.HG007.GRCh38.HiFi.minimap2.bam --access_path ./src/pre_process/hg38.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_256_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/real_trio/ --genome_path ~/workspace/ref/chr1-X.fa --sample_name test_gt --min_supp 10 --preset hifi --process_num 48 --detect_mode genotype --region chr11:34680676-34689676



# # # # # # STEP: five trios HiFi
python SVision-pro.py --target_path /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/GIAB/HG002/GIAB.HG002.GRCh38.HiFi.minimap2.bam --base_path /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/GIAB/HG003/GIAB.HG003.GRCh38.HiFi.minimap2.bam /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/GIAB/HG004/GIAB.HG004.GRCh38.HiFi.minimap2.bam --access_path ./src/pre_process/hg38.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_256_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/real_trio/sv_call_release/ --genome_path /data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta --sample_name release_hg002-3-4_256 --min_supp 10 --preset hifi --process_num 48 --detect_mode denovo
python SVision-pro.py --target_path /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/GIAB/HG005/GIAB.HG005.GRCh38.HiFi.minimap2.bam --base_path /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/GIAB/HG006/GIAB.HG006.GRCh38.HiFi.minimap2.bam /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/GIAB/HG007/GIAB.HG007.GRCh38.HiFi.minimap2.bam --access_path ./src/pre_process/hg38.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_256_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/real_trio/sv_call_release/ --genome_path /data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta --sample_name release_hg005-6-7_256 --min_supp 10 --preset hifi --process_num 48 --detect_mode denovo
python SVision-pro.py --target_path /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/HGSVC/HG00733/HGSVC.HG00733.GRCh38.HiFi.minimap2.bam --base_path /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/HGSVC/HG00731/HGSVC.HG00731.GRCh38.HiFi.minimap2.bam /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/HGSVC/HG00732/HGSVC.HG00732.GRCh38.HiFi.minimap2.bam --access_path ./src/pre_process/hg38.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_256_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/real_trio/sv_call_release/ --genome_path /data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta --sample_name release_hg00733-1-2_256 --min_supp 10 --preset hifi --process_num 48 --detect_mode denovo
python SVision-pro.py --target_path /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/HGSVC/HG00514/HGSVC.HG00514.GRCh38.HiFi.minimap2.bam --base_path /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/HGSVC/HG00512/HGSVC.HG00512.GRCh38.HiFi.minimap2.bam /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/HGSVC/HG00513/HGSVC.HG00513.GRCh38.HiFi.minimap2.bam --access_path ./src/pre_process/hg38.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_256_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/real_trio/sv_call_release/ --genome_path /data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta --sample_name release_hg00514-2-3_256 --min_supp 10 --preset hifi --process_num 48 --detect_mode denovo
python SVision-pro.py --target_path /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/HGSVC/NA19240/HGSVC.NA19240.GRCh38.HiFi.minimap2.bam --base_path /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/HGSVC/NA19238/HGSVC.NA19238.GRCh38.HiFi.minimap2.bam /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/HGSVC/NA19239/HGSVC.NA19239.GRCh38.HiFi.minimap2.bam --access_path ./src/pre_process/hg38.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_256_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/real_trio/sv_call_release/ --genome_path /data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta --sample_name release_na19240-38-39_256 --min_supp 10 --preset hifi --process_num 48 --detect_mode denovo
python SVision-pro.py --target_path /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/ChineseQuartet/LCL5/ChineseQuartet.LCL5.GRCh38.HiFi.minimap2.bam --base_path /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/ChineseQuartet/LCL6/ChineseQuartet.LCL6.GRCh38.HiFi.minimap2.bam /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/ChineseQuartet/LCL7/ChineseQuartet.LCL7.GRCh38.HiFi.minimap2.bam /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/ChineseQuartet/LCL8/ChineseQuartet.LCL8.GRCh38.HiFi.minimap2.bam --access_path ./src/pre_process/hg38.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_256_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/real_trio/sv_call_release/ --genome_path /data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta --sample_name release_lcl5-6-7-8_256 --min_supp 10 --preset hifi --process_num 48 --detect_mode denovo


sniffles --input HG002.sniffles2.s10.m50.snf HG003.sniffles2.s10.m50.snf HG004.sniffles2.s10.m50.snf --vcf sniffles2.hg002-3-4.s10.vcf --minsupport 10 --minsvlen 50 -t 48 --mapq 20
sniffles --input HG005.sniffles2.s10.m50.snf HG006.sniffles2.s10.m50.snf HG007.sniffles2.s10.m50.snf --vcf sniffles2.hg005-6-7.s10.vcf --minsupport 10 --minsvlen 50 -t 48 --mapq 20
sniffles --input NA19240.sniffles2.s10.m50.snf NA19238.sniffles2.s10.m50.snf NA19239.sniffles2.s10.m50.snf --vcf sniffles2.na19240-38-39.s10.vcf --minsupport 10 --minsvlen 50 -t 48 --mapq 20
sniffles --input HG00514.sniffles2.s10.m50.snf HG00512.sniffles2.s10.m50.snf HG00513.sniffles2.s10.m50.snf --vcf sniffles2.hg00514-2-3.s10.vcf --minsupport 10 --minsvlen 50 -t 48 --mapq 20
sniffles --input HG00733.sniffles2.s10.m50.snf HG00731.sniffles2.s10.m50.snf HG00732.sniffles2.s10.m50.snf --vcf sniffles2.hg00733-1-2.s10.vcf --minsupport 10 --minsvlen 50 -t 48 --mapq 20
sniffles --input LCL5.sniffles2.s10.m50.snf LCL6.sniffles2.s10.m50.snf LCL7.sniffles2.s10.m50.snf LCL8.sniffles2.s10.m50.snf --vcf sniffles2.lcl5-6-7-8.s10.vcf --minsupport 10 --minsvlen 50 -t 48 --mapq 20


========================svision========================

ls HG002.svision.s10.vcf HG003.svision.s10.vcf HG004.svision.s10.vcf > sample_files
jasmine --output_genotypes file_list=sample_files out_file=svision_jasmine.hg002-3-4.s10.vcf genome_file=/data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta
~/workspace/tools/SURVIVOR-master/Debug/SURVIVOR merge sample_files 1000 1 1 1 0 50 svision_survivor.hg002-3-4.s10.vcf

ls HG005.svision.s10.vcf HG006.svision.s10.vcf HG007.svision.s10.vcf > sample_files
jasmine --output_genotypes file_list=sample_files out_file=svision_jasmine.hg005-6-7.s10.vcf genome_file=/data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta
~/workspace/tools/SURVIVOR-master/Debug/SURVIVOR merge sample_files 1000 1 1 1 0 50 svision_survivor.hg005-6-7.s10.vcf

ls NA19240.svision.s10.vcf > sample_files
ls NA19238.svision.s10.vcf >> sample_files
ls NA19239.svision.s10.vcf >> sample_files
jasmine --output_genotypes file_list=sample_files out_file=svision_jasmine.na19240-38-39.s10.vcf genome_file=/data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta
~/workspace/tools/SURVIVOR-master/Debug/SURVIVOR merge sample_files 1000 1 1 1 0 50 svision_survivor.na19240-38-39.s10.vcf

ls HG00514.svision.s10.vcf > sample_files
ls HG00512.svision.s10.vcf >> sample_files
ls HG00513.svision.s10.vcf >> sample_files
jasmine --output_genotypes file_list=sample_files out_file=svision_jasmine.hg00514-2-3.s10.vcf genome_file=/data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta
~/workspace/tools/SURVIVOR-master/Debug/SURVIVOR merge sample_files 1000 1 1 1 0 50 svision_survivor.hg00514-2-3.s10.vcf

ls HG00733.svision.s10.vcf > sample_files
ls HG00731.svision.s10.vcf >> sample_files
ls HG00732.svision.s10.vcf >> sample_files
jasmine --output_genotypes file_list=sample_files out_file=svision_jasmine.hg00733-1-2.s10.vcf genome_file=/data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta
~/workspace/tools/SURVIVOR-master/Debug/SURVIVOR merge sample_files 1000 1 1 1 0 50 svision_survivor.hg00733-1-2.s10.vcf

ls LCL5.svision.s10.vcf LCL6.svision.s10.vcf LCL7.svision.s10.vcf LCL8.svision.s10.vcf > sample_files
jasmine --output_genotypes file_list=sample_files out_file=svision_jasmine.lcl5-6-7-8.s10.vcf genome_file=/data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta
~/workspace/tools/SURVIVOR-master/Debug/SURVIVOR merge sample_files 1000 1 1 1 0 50 svision_survivor.lcl5-6-7-8.s10.vcf


========================cutesv========================


ls HG002.cutesv.s10.m50.vcf HG003.cutesv.s10.m50.vcf HG004.cutesv.s10.m50.vcf > sample_files
jasmine --output_genotypes file_list=sample_files out_file=cutesv_jasmine.hg002-3-4.s10.vcf genome_file=/data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta
~/workspace/tools/SURVIVOR-master/Debug/SURVIVOR merge sample_files 1000 1 1 1 0 50 cutesv_survivor.hg002-3-4.s10.vcf

ls HG005.cutesv.s10.m50.vcf HG006.cutesv.s10.m50.vcf HG007.cutesv.s10.m50.vcf > sample_files
jasmine --output_genotypes file_list=sample_files out_file=cutesv_jasmine.hg005-6-7.s10.vcf genome_file=/data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta
~/workspace/tools/SURVIVOR-master/Debug/SURVIVOR merge sample_files 1000 1 1 1 0 50 cutesv_survivor.hg005-6-7.s10.vcf

ls NA19240.cutesv.s10.m50.vcf > sample_files
ls NA19238.cutesv.s10.m50.vcf >> sample_files
ls NA19239.cutesv.s10.m50.vcf >> sample_files
jasmine --output_genotypes file_list=sample_files out_file=cutesv_jasmine.na19240-38-39.s10.vcf genome_file=/data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta
~/workspace/tools/SURVIVOR-master/Debug/SURVIVOR merge sample_files 1000 1 1 1 0 50 cutesv_survivor.na19240-38-39.s10.vcf

ls HG00514.cutesv.s10.m50.vcf > sample_files
ls HG00512.cutesv.s10.m50.vcf >> sample_files
ls HG00513.cutesv.s10.m50.vcf >> sample_files
jasmine --output_genotypes file_list=sample_files out_file=cutesv_jasmine.hg00514-2-3.s10.vcf genome_file=/data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta
~/workspace/tools/SURVIVOR-master/Debug/SURVIVOR merge sample_files 1000 1 1 1 0 50 cutesv_survivor.hg00514-2-3.s10.vcf

ls HG00733.cutesv.s10.m50.vcf > sample_files
ls HG00731.cutesv.s10.m50.vcf >> sample_files
ls HG00732.cutesv.s10.m50.vcf >> sample_files
jasmine --output_genotypes file_list=sample_files out_file=cutesv_jasmine.hg00733-1-2.s10.vcf genome_file=/data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta
~/workspace/tools/SURVIVOR-master/Debug/SURVIVOR merge sample_files 1000 1 1 1 0 50 cutesv_survivor.hg00733-1-2.s10.vcf

ls LCL5.cutesv.s10.m50.vcf LCL6.cutesv.s10.m50.vcf LCL7.cutesv.s10.m50.vcf LCL8.cutesv.s10.m50.vcf > sample_files
jasmine --output_genotypes file_list=sample_files out_file=cutesv_jasmine.lcl5-6-7-8.s10.vcf genome_file=/data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta
~/workspace/tools/SURVIVOR-master/Debug/SURVIVOR merge sample_files 1000 1 1 1 0 50 cutesv_survivor.lcl5-6-7-8.s10.vcf


========================debreak========================

ls HG002.debreak.s10.m50.vcf HG003.debreak.s10.m50.vcf HG004.debreak.s10.m50.vcf > sample_files
jasmine --output_genotypes file_list=sample_files out_file=debreak_jasmine.hg002-3-4.s10.vcf genome_file=/data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta
~/workspace/tools/SURVIVOR-master/Debug/SURVIVOR merge sample_files 1000 1 1 1 0 50 debreak_survivor.hg002-3-4.s10.vcf

ls HG005.debreak.s10.m50.vcf HG006.debreak.s10.m50.vcf HG007.debreak.s10.m50.vcf > sample_files
jasmine --output_genotypes file_list=sample_files out_file=debreak_jasmine.hg005-6-7.s10.vcf genome_file=/data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta
~/workspace/tools/SURVIVOR-master/Debug/SURVIVOR merge sample_files 1000 1 1 1 0 50 debreak_survivor.hg005-6-7.s10.vcf

ls NA19240.debreak.s10.m50.vcf > sample_files
ls NA19238.debreak.s10.m50.vcf >> sample_files
ls NA19239.debreak.s10.m50.vcf >> sample_files
jasmine --output_genotypes file_list=sample_files out_file=debreak_jasmine.na19240-38-39.s10.vcf genome_file=/data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta
~/workspace/tools/SURVIVOR-master/Debug/SURVIVOR merge sample_files 1000 1 1 1 0 50 debreak_survivor.na19240-38-39.s10.vcf

ls HG00514.debreak.s10.m50.vcf > sample_files
ls HG00512.debreak.s10.m50.vcf >> sample_files
ls HG00513.debreak.s10.m50.vcf >> sample_files
jasmine --output_genotypes file_list=sample_files out_file=debreak_jasmine.hg00514-2-3.s10.vcf genome_file=/data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta
~/workspace/tools/SURVIVOR-master/Debug/SURVIVOR merge sample_files 1000 1 1 1 0 50 debreak_survivor.hg00514-2-3.s10.vcf

ls HG00733.debreak.s10.m50.vcf > sample_files
ls HG00731.debreak.s10.m50.vcf >> sample_files
ls HG00732.debreak.s10.m50.vcf >> sample_files
jasmine --output_genotypes file_list=sample_files out_file=debreak_jasmine.hg00733-1-2.s10.vcf genome_file=/data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta
~/workspace/tools/SURVIVOR-master/Debug/SURVIVOR merge sample_files 1000 1 1 1 0 50 debreak_survivor.hg00733-1-2.s10.vcf

ls LCL5.debreak.s10.m50.vcf LCL6.debreak.s10.m50.vcf LCL7.debreak.s10.m50.vcf LCL8.debreak.s10.m50.vcf > sample_files
jasmine --output_genotypes file_list=sample_files out_file=debreak_jasmine.lcl5-6-7-8.s10.vcf genome_file=/data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta
~/workspace/tools/SURVIVOR-master/Debug/SURVIVOR merge sample_files 1000 1 1 1 0 50 debreak_survivor.lcl5-6-7-8.s10.vcf



# # # # # # five trios ONT

python SVision-pro.py --target_path /data/DATA/ChineseQuartet/ref_based_analysis/aligned_reads/ChineseQuartet/LCL5/ChineseQuartet.LCL5.GRCh38.ONT.minimap2.bam --base_path /data/DATA/ChineseQuartet/ref_based_analysis/aligned_reads/ChineseQuartet/LCL6/ChineseQuartet.LCL6.GRCh38.ONT.minimap2.bam /data/DATA/ChineseQuartet/ref_based_analysis/aligned_reads/ChineseQuartet/LCL7/ChineseQuartet.LCL7.GRCh38.ONT.minimap2.bam /data/DATA/ChineseQuartet/ref_based_analysis/aligned_reads/ChineseQuartet/LCL8/ChineseQuartet.LCL8.GRCh38.ONT.minimap2.bam --access_path ./src/pre_process/hg38.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_256_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/real_trio/sv_call_release/ont --genome_path /data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta --sample_name release_lcl5-6-7-8_ONT_256 --min_supp 10 --preset error-prone --process_num 48 --detect_mode denovo
python SVision-pro.py --target_path /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/GIAB/HG002/GIAB.HG002.GRCh38.ONTUL.minimap2.bam --base_path /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/GIAB/HG003/GIAB.HG003.GRCh38.ONTUL.minimap2.bam /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/GIAB/HG004/GIAB.HG004.GRCh38.ONTUL.minimap2.bam  --access_path ./src/pre_process/hg38.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_256_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/real_trio/sv_call_release/ont/ --genome_path /data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta --sample_name release_hg002-3-4_ONT_256 --min_supp 10 --preset error-prone --process_num 48 --detect_mode denovo
python SVision-pro.py --target_path /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/GIAB/HG005/GIAB.HG005.GRCh38.ONTUL.minimap2.bam --base_path /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/GIAB/HG006/GIAB.HG006.GRCh38.ONTUL.minimap2.bam /data/DATA/Human_genomes/reference_based_analysis/aligned_reads/GIAB/HG007/GIAB.HG007.GRCh38.ONTUL.minimap2.bam  --access_path ./src/pre_process/hg38.access.10M.bed --model_path /data/home/songbo/workspace/svision-pro/sim_train_ccs/models/model_liteunet_256_8_16_32_32_32.pth --out_path /data/home/songbo/workspace/svision-pro/real_trio/sv_call_release/ont/ --genome_path /data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta --sample_name release_hg005-6-7_ONT_256 --min_supp 10 --preset error-prone --process_num 48 --detect_mode denovo


sniffles --input HG002-ONT.sniffles2.s10.m50.snf HG003-ONT.sniffles2.s10.m50.snf HG004-ONT.sniffles2.s10.m50.snf --vcf sniffles2.hg002-3-4_ONT.s10.vcf --minsupport 10 --minsvlen 50 -t 48 --mapq 20
sniffles --input HG005-ONT.sniffles2.s10.m50.snf HG006-ONT.sniffles2.s10.m50.snf HG007-ONT.sniffles2.s10.m50.snf --vcf sniffles2.hg005-6-7_ONT.s10.vcf --minsupport 10 --minsvlen 50 -t 48 --mapq 20
sniffles --input LCL5-ONT.sniffles2.s10.m50.snf LCL6-ONT.sniffles2.s10.m50.snf LCL7-ONT.sniffles2.s10.m50.snf LCL8-ONT.sniffles2.s10.m50.snf --vcf sniffles2.lcl5-6-7-8_ONT.s10.vcf --minsupport 10 --minsvlen 50 -t 48 --mapq 20


========================svision========================

ls HG002-ONT.svision.s10.vcf HG003-ONT.svision.s10.vcf HG004-ONT.svision.s10.vcf > sample_files
jasmine --output_genotypes file_list=sample_files out_file=svision_jasmine.hg002-3-4_ONT.s10.vcf genome_file=/data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta
~/workspace/tools/SURVIVOR-master/Debug/SURVIVOR merge sample_files 1000 1 1 1 0 50 svision_survivor.hg002-3-4_ONT.s10.vcf

ls HG005-ONT.svision.s10.vcf HG006-ONT.svision.s10.vcf HG007-ONT.svision.s10.vcf > sample_files
jasmine --output_genotypes file_list=sample_files out_file=svision_jasmine.hg005-6-7_ONT.s10.vcf genome_file=/data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta
~/workspace/tools/SURVIVOR-master/Debug/SURVIVOR merge sample_files 1000 1 1 1 0 50 svision_survivor.hg005-6-7_ONT.s10.vcf

ls LCL5-ONT.svision.s10.vcf LCL6-ONT.svision.s10.vcf LCL7-ONT.svision.s10.vcf LCL8-ONT.svision.s10.vcf > sample_files
jasmine --output_genotypes file_list=sample_files out_file=svision_jasmine.lcl5-6-7-8_ONT.s10.vcf genome_file=/data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta
~/workspace/tools/SURVIVOR-master/Debug/SURVIVOR merge sample_files 1000 1 1 1 0 50 svision_survivor.lcl5-6-7-8_ONT.s10.vcf


========================cutesv========================

ls HG002-ONT.cutesv.s10.m50.vcf HG003-ONT.cutesv.s10.m50.vcf HG004-ONT.cutesv.s10.m50.vcf > sample_files
jasmine --output_genotypes file_list=sample_files out_file=cutesv_jasmine.hg002-3-4_ONT.s10.vcf genome_file=/data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta
~/workspace/tools/SURVIVOR-master/Debug/SURVIVOR merge sample_files 1000 1 1 1 0 50 cutesv_survivor.hg002-3-4_ONT.s10.vcf

ls HG005-ONT.cutesv.s10.m50.vcf HG006-ONT.cutesv.s10.m50.vcf HG007-ONT.cutesv.s10.m50.vcf > sample_files
jasmine --output_genotypes file_list=sample_files out_file=cutesv_jasmine.hg005-6-7_ONT.s10.vcf genome_file=/data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta
~/workspace/tools/SURVIVOR-master/Debug/SURVIVOR merge sample_files 1000 1 1 1 0 50 cutesv_survivor.hg005-6-7_ONT.s10.vcf

ls LCL5-ONT.cutesv.s10.m50.vcf LCL6-ONT.cutesv.s10.m50.vcf LCL7-ONT.cutesv.s10.m50.vcf LCL8-ONT.cutesv.s10.m50.vcf > sample_files
jasmine --output_genotypes file_list=sample_files out_file=cutesv_jasmine.lcl5-6-7-8_ONT.s10.vcf genome_file=/data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta
~/workspace/tools/SURVIVOR-master/Debug/SURVIVOR merge sample_files 1000 1 1 1 0 50 cutesv_survivor.lcl5-6-7-8_ONT.s10.vcf


========================debreak========================
ls HG002-ONT.debreak.s10.m50.vcf HG003-ONT.debreak.s10.m50.vcf HG004-ONT.debreak.s10.m50.vcf > sample_files
jasmine --output_genotypes file_list=sample_files out_file=debreak_jasmine.hg002-3-4_ONT.s10.vcf genome_file=/data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta
~/workspace/tools/SURVIVOR-master/Debug/SURVIVOR merge sample_files 1000 1 1 1 0 50 debreak_survivor.hg002-3-4_ONT.s10.vcf

ls HG005-ONT.debreak.s10.m50.vcf HG006-ONT.debreak.s10.m50.vcf HG007-ONT.debreak.s10.m50.vcf > sample_files
jasmine --output_genotypes file_list=sample_files out_file=debreak_jasmine.hg005-6-7_ONT.s10.vcf genome_file=/data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta
~/workspace/tools/SURVIVOR-master/Debug/SURVIVOR merge sample_files 1000 1 1 1 0 50 debreak_survivor.hg005-6-7_ONT.s10.vcf

ls LCL5-ONT.debreak.s10.m50.vcf LCL6-ONT.debreak.s10.m50.vcf LCL7-ONT.debreak.s10.m50.vcf LCL8-ONT.debreak.s10.m50.vcf > sample_files
jasmine --output_genotypes file_list=sample_files out_file=debreak_jasmine.lcl5-6-7-8_ONT.s10.vcf genome_file=/data/DATA/Human_genomes/reference_based_analysis/ref/GRCh38.fasta
~/workspace/tools/SURVIVOR-master/Debug/SURVIVOR merge sample_files 1000 1 1 1 0 50 debreak_survivor.lcl5-6-7-8_ONT.s10.vcf


# # # # # # denovo overlap
bedtools intersect -r -f 0.5  -wao -a collection.sniffles2.denovo.bed -b ../release_six-trios_256.svision_pro_v1.2.s10.bed > collection.sniffles2.denovo.overlap_svision_pro.bed
"""


