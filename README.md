## SVision-pro-Utils
The utilities and scripts in SVision-pro paper, including 
* Running commands 
* Comparison criteria
* Model training and interpreting 
* Figure plotting 

## License

SVision-pro and its utilities are free for non-commercial use by academic, government, and non-profit/not-for-profit institutions. A commercial version of the software is available and licensed through Xi’an Jiaotong University.
For more information, please contact with Songbo Wang (songbowang125@163.com) or Kai Ye (kaiye@xjtu.edu.cn).


## 0. Operation systems
Linux-based systems required, including but are not limited to: 
* MacOS, Big Sur
* Ubuntu (including Windows Subsystem)
* CentOS

## 1. Running commands
Dependencies:
```commandline
* SV callers: SVision-pro v1.6; SVision v1.3.9; Sniffle2 v2.0.7; cuteSV v2.0.2; pbsv v2.9.0; debrea v1.0.2; SVDSS v1.0.5; nanomonsv v0.5.0
* Merging approach: Jasmine v1.1.5; SURVIVOR v1.0.7
```

Explaination:
```commandline
 ./compare/cmds.py  Running commands of callers, mering approaches and other involved tools are listed in this file.
```

## 2. Comparison criteria
Dependencies:
```commandline
* python 3.7.9
* pysam
```

Explaination:
```commandline
./compare/run_trio_denovo.py        Extract de novo SVs from callsets generated by different approaches.
./compare/run_trio_mendelian.py     Calculate the Mendelian consistency on six familiy datasets and the twin discordancy in ChineseQuartet.
./compare/sim_csv_process.py        Calculate the CSV type/subcomponent accuracy based on the output of Turvari.
./compare/sim_denovo_generate.py    Generate the simulated trio dataset harboring de novo/inherited SSVs and CSVs. Each SV is assigend with a denovo or inherited flag and is further implanted into the two alleles of child, father and mother genome.
./compare/sim_denovo_process.py     Evaluate the  simulated trio dataset harboring de novo/inherited SSVs and CSVs. Each SV is assigend with a denovo or inherited flag and is further implanted into the two alleles of child, father and mother genome.
./compare/sim_somatic_generate_xxx.py   Assign SSVs and CSVs with low AFs and implant them into tumor genome by altering the normal reads in each SV region.
./compare/sim_somatic_process.py    Evaluate the low-AF somatic SSV/CSV detection accuracy. 
./compare/data_HG002_SVs_Tier1_v0.6.final.bed               The high confident deletions and insertion of HG002. The 9,641 calls are merged from Truvari TP-base.vcf and FN.vcf 
./compare/data_data_benchmark_origin_clone0_300.csv.bed     The simualted 3,000 CSV obtained from SVision paper
```

## 3. Model training and interpreting
Dependencies:
```commandline
numpy       1.19.5
pytorch     1.10.1
torchvision 0.11.2
sklearn     0.24.2
torchsummary
PIL
```
Explaination:
```commandline
./model/model_xxx.py    Neural network architecture of each model.
./model/op_dataset.py   Split dataset into training/validation by k-folds and load images for model.
./model/op_summary.py   Summary the parameter size of each model.
./model/op_train.py     Train and validata models. Trainning parameters can be customized.
./model/op_visulaize.py Interpret model by Layer Grad Cam and Feature Ablation. Visualize the attribution map.
```

## 4. Figure plotting
Dependencies:
```commandline
matplotlib  3.5.1
seaborn     0.11.2
pandas      1.4.4
numpy       1.23.3
scipy       1.7.1
pylab
```
Explaination:
```commandline
./plot/plot_xxx.py  Figure plotting scripts in the paper
```
## Contact
If you have any questions about this repository, please feel free to contact: songbowang125@163.com