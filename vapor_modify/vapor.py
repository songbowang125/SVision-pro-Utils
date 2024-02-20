#!/usr/bin/env python

#!python
#For debug use only
#Test deletions:
#command='Pacbio.Validator.py vcf --sv-input /mnt/EXT/Mills-scratch2/Xuefang/Pacbio.Validator/vcf_files/NA12878_Delly_DEL.vcf --output-path /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/SVelter_rec2/PacValTest/ --pacbio-input /mnt/EXT/Mills-scratch2/Xuefang/NA12878.Pacbio/alignment/sorted_final_merged.bam --reference /mnt/EXT/Mills-scratch2/reference/hg19_platinum/genome.fa'
#Test inversions:
#command='Pacbio.Validator.py vcf --sv-input /mnt/EXT/Mills-scratch2/Xuefang/Pacbio.Validator/vcf_files/NA12878_Delly_INV.vcf --output-path /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/SVelter_rec2/PacValTest/ --pacbio-input /mnt/EXT/Mills-scratch2/Xuefang/NA12878.Pacbio/alignment/sorted_final_merged.bam --reference /mnt/EXT/Mills-scratch2/reference/hg19_platinum/genome.fa'
#Test inversions:
#command='vapor.py vcf --sv-input /scratch/remills_flux/xuefzhao/SV_discovery_index/download/SVelter.version17/All.Individuals.vcf --output-path /scratch/remills_flux/xuefzhao/SV_discovery_index/download/SVelter.version17/vapor_test/  --reference /scratch/remills_flux/xuefzhao/reference/GRCh38.1KGP/GRCh38_full_analysis_set_plus_decoy_hla.fa --pacbio-input /nfs/turbo/remillsscr/scratch2_trans/datasets/1000genomes/vol1/ftp/data_collections/hgsv_sv_discovery/PacBio/alignment/HG00512.XXX.bam'
#command='vapor.py vcf --sv-input test.vcf --output-path ./test/  --reference /scratch/remills_flux/xuefzhao/reference/GRCh38.1KGP/GRCh38_full_analysis_set_plus_decoy_hla.fa --pacbio-input /nfs/turbo/remillsscr/scratch2_trans/datasets/1000genomes/vol1/ftp/data_collections/hgsv_sv_discovery/PacBio/alignment/HG00512.XXX.bam'
#Test SVelter:
#sys.argv=command.split()
from __future__ import absolute_import
from __future__ import print_function

import os
import re
import sys
import argparse
from Simple_function import *

def bed_info_readin(bed_input,out_path):
	out_path=path_modify(out_path)
	path_mkdir(out_path)
	out=[]
	#test=[]
	fin=open(bed_input)
	for line in fin:
		pin=line.strip().split()
		#test.append(pin)
		if 'DUP' in pin[4] or 'duplication' in pin[4]:	
			out.append([pin[0]]+[int(i) for i in pin[1:3]]+[pin[3]]+['a/a','a/aa'])
			#if not '_' in pin[3]:	out.append([pin[0]]+[int(i) for i in pin[1:3]]+['a/a','a/aa'])
			#else:					out.append([pin[0]]+[int(i) for i in pin[1:3]]+['a/a','/'.join([ ''.join(['a' for i in range(int(pin[3].split('_')[1].replace('<CN','').replace('>','')))]),''.join(['a' for i in range(int(pin[3].split('_')[2].replace('<CN','').replace('>','')))])])])
		elif 'DEL' in pin[4] or 'deletion' in pin[4]:	
			out.append([pin[0]]+[int(i) for i in pin[1:3]]+[pin[3]]+['a/a','/a'])
		elif 'INV' in pin[4] or 'inversion' in pin[4]:	
			out.append([pin[0]]+[int(i) for i in pin[1:3]]+[pin[3]]+['a/a','a/a^'])
		elif 'INS' in pin[4] or 'ALU' in pin[4] or 'HERVK' in pin[4] or 'LINE1' in pin[4] or 'SVA' in pin[4] or 'insertion' in pin[4]:
			if len(pin)>5:		
				out.append([pin[0],int(pin[1]),int(pin[2]),pin[3],pin[5],'INS'])
			else:				
				if not '_' in pin[4]:	
					continue
				else:					
					if pin[4].split('_')[1].isdigit():	
						out.append([pin[0],int(pin[1]),int(pin[2]),pin[3],int(pin[4].split('_')[1]),'INS'])
					else:
						out.append([pin[0],int(pin[1]),int(pin[2]),pin[3],pin[4].split('_')[1],'INS'])
	return out

def melt_info_readin(melt_input_prefix,out_path,sample_name,bam_in,ref):
	#eg of melt_input_prefix='/scratch/remills_flux/xuefzhao/VaPoR/MELT/MELT_with_seq/HG00512.CHS.melt.LINE1.pacbio.20160504.sites'
	#eg of out_path='/scratch/remills_flux/xuefzhao/VaPoR/MELT/MELT_with_seq/VaPoR_Vali_MELT'
	#eg of sample_name='HG00512'
	out_path=path_modify(out_path)
	path_mkdir(out_path)
	vcf_name=melt_input_prefix+'.vcf'
	seq_name=melt_input_prefix+'.fa'
	fin=open(vcf_name)
	out=[]
	write_output_initiate(melt_input_prefix+'.vapor.py')
	plt_li=0
	for line in fin:
		pin=line.strip().split()
		if not pin[0][0]=='#':
			key_event='_'.join(pin[:2])
			ffa=os.popen(r'''samtools faidx %s %s'''%(seq_name,key_event))
			ins_seq=''
			header=ffa.readline().strip().split()
			for pfa in ffa:
				ins_seq+=pfa.strip()
			ffa.close()
			if ins_seq=='':	ins_seq=''.join(['X' for i in range(INS_length_detect(pin))])
			if not ins_seq=='' and 'INS' in pin[3]:
				POLARITY=polarity_detect(pin)
				ins_seq=ins_seq.replace('N','X')
				plt_li+=1
				vapor_score_event=vapor_simple_ins_Vapor(num_reads_cff,plt_li,bam_in,ref,key_event,ins_seq,out_path+sample_name+'.INS.'+key_event.replace(':','__').replace(':','__')+'.png',POLARITY)
				write_output_main(melt_input_prefix+'.vapor.py',result_organize_ins([key_event,vapor_score_event]))
	fin.close()
	#command='vapor.py ins --sv-input-prefix /scratch/remills_flux/xuefzhao/VaPoR/MELT/MELT_with_seq/HG00512.CHS.melt.LINE1.pacbio.20160504.sites --output-path /scratch/remills_flux/xuefzhao/VaPoR/MELT/MELT_with_seq/VaPoR_Vali_MELT/HG00512/		--reference /scratch/remills_flux/xuefzhao/reference/GRCh38.1KGP/GRCh38_full_analysis_set_plus_decoy_hla.fa		--pacbio-input /nfs/turbo/remillsscr/scratch2_trans/datasets/1000genomes/vol1/ftp/data_collections/hgsv_sv_discovery/PacBio/alignment/HG00512.XXX.bam'

def block_reorganize(block_hash):
	#eg of block_hash={'chr12': [['chr12', 123429440, 123430134, 'del'], ['chr12', 123425669, 123429440, 'inv']]}
	if len(list(block_hash.keys()))==1:
		for k1 in list(block_hash.keys()):
			start=[i[1] for i in block_hash[k1]]
			start_new=sorted(start)
			start_rec=[start.index(i) for i in start_new]
			temp=[block_hash[k1][i] for i in start_rec]
			temp_new=[]
			for k2 in temp:
				if not k2 in temp_new:	temp_new.append(k2)
			return temp_new
	else:
		return 'error'

def del_inv_interprete(pin):
	out={}
	for x in pin[7].split(';'):
		if 'del=' in x or 'DEL=' in x:
			del_block=[x.split('=')[1].split(':')[0]]+[int(i) for i in x.split('=')[1].split(':')[1].split('-')]
			if not del_block[0] in list(out.keys()):	out[del_block[0]]=[]
			out[del_block[0]].append(del_block+['del'])
		elif 'inv=' in x or 'INV=' in x:
			inv_block=[x.split('=')[1].split(':')[0]]+[int(i) for i in x.split('=')[1].split(':')[1].split('-')]
			if not inv_block[0] in list(out.keys()):	out[inv_block[0]]=[]
			out[inv_block[0]].append(inv_block+['inv'])
	ordered_list=block_reorganize(out)
	return ordered_list

def dup_inv_interprete(pin):
	dup_seg=[pin[0],int(pin[1])]
	insert_pos=[]
	for x in pin[7].split(';'):
		if 'END=' in x:
			dup_seg.append(int(x.split('=')[1]))
		if 'insert_point' in x or 'INSERT_POINT' in x:
			insert_pos=x.split('=')[1].split(':')
	if len(insert_pos)>1:
		out=dup_seg+[insert_pos[0],int(insert_pos[1])]
	else:
		out='error'
	return out

def vcf_list_readin(file_in):
	out={}
	fin=open(file_in)
	rec=-1
	vcf_rec_hash={}
	for line in fin:
		rec+=1
		pin=line.strip().split()
		if not pin[0][0]=='#':
			#if pin[6]=='PASS':
				pin[7]=pin[7].replace('MERGE_TYPE=','SVTYPE=')
				sv_type=svtype_extract(pin)
				sv_pos=chr_start_end_extract(pin)
				#if sv_type in ['tandup','TANDUP']:
				#	genotypes=genoCN_extract(pin)
				#else:
				#	genotypes=genotype_extract(pin)
				genotypes=[0,1]
				if not sum(genotypes)==0:
					if sv_type in ['del','DEL','deletion']:
						if not 'DEL' in list(out.keys()):		out['DEL']=[]
						if not sv_pos in out['DEL']:	
							out['DEL'].append(sv_pos)
							vcf_rec_hash[rec]=':'.join([str(i) for i in sv_pos]+['DEL'])
					elif sv_type in ['inv','INV','inversion']:
						if not 'INV' in list(out.keys()):		out['INV']=[]
						if not sv_pos in out['INV']:	
							out['INV'].append(sv_pos)
							vcf_rec_hash[rec]=':'.join([str(i) for i in sv_pos]+['INV'])
					elif sv_type in ['ins','INS','insertion','LINE1','SVA','ALU','HERVK']:
						SV_len=int(sv_len_extract(pin))
						INS_Seq=sv_seq_extract(pin)
						if SV_len>0:
							if not 'INS' in list(out.keys()):		out['INS']=[]
							if not sv_pos in out['INS']:	
								out['INS'].append(sv_pos[:2]+[SV_len,INS_Seq])
								vcf_rec_hash[rec]=':'.join([str(i) for i in sv_pos[:2]+[SV_len]]+['INS'])
					elif sv_type in ['disdup','DISDUP','dis-dup']:
						insert_point=sv_insert_point_define(pin)
						if not 'DISDUP' in list(out.keys()):		out['DISDUP']=[]
						if not sv_pos in out['DISDUP']:	
							out['DISDUP'].append(sv_pos+insert_point)
							vcf_rec_hash[rec]=':'.join([str(i) for i in sv_pos+insert_point]+['DISDUP'])
					elif sv_type in ['DEL_INV','del_inv']:
						if not 'DEL_INV' in list(out.keys()):		out['DEL_INV']=[]
						del_inv_info= del_inv_interprete(pin)
						if not del_inv_info=='error':
							if not del_inv_info in out['DEL_INV']:	
								out['DEL_INV'].append(del_inv_info)
								vcf_rec_hash[rec]=':'.join(['_'.join([str(i) for i in j]) for j in del_inv_info]+['DEL_INV'])
					elif sv_type in ['DUP_INV','dup_inv']:	
						if not 'DUP_INV' in list(out.keys()):		out['DUP_INV']=[]
						dup_inv_info=dup_inv_interprete(pin)
						if not dup_inv_info=='error':
							if not dup_inv_info in out['DUP_INV']:	
								out['DUP_INV'].append(dup_inv_info)
								vcf_rec_hash[rec]=':'.join([str(i) for i in dup_inv_info+['DUP_INV']])
					elif sv_type in ['tandup','TANDUP','DUP']:
						if not 'TANDUP' in list(out.keys()):	out['TANDUP']=[]
						if not sv_pos in out['TANDUP']:	
							out['TANDUP'].append(sv_pos)
							vcf_rec_hash[rec]=':'.join([str(i) for i in sv_pos]+['TANDUP'])
					elif sv_type in ['CNV','CSV','CPX']: continue
					else:
						if 'Other=' in pin[7]:
							info=[i for i in pin[7].split(';') if i[:6]=='Other=']
						elif 'OTHER=' in pin[7]:
							info=[i for i in pin[7].split(';') if i[:6]=='OTHER=']
						else: continue
						sv_other_info=info[0].split('=')[1].split('_')
						if not 'Other' in list(out.keys()):	out['Other']=[]
						if not ['_'.join(i.split('/')) for i in sv_other_info[:2]]+sv_other_info[2].split(':') in out['Other']:	
							out['Other'].append(['_'.join(i.split('/')) for i in sv_other_info[:2]]+sv_other_info[2].split(':'))
							vcf_rec_hash[rec]=':'.join([str(i) for i in ['_'.join(i.split('/')) for i in sv_other_info[:2]]+sv_other_info[2].split(':')+['CANNOT_CLASSIFY']])
	fin.close()
	return [out,vcf_rec_hash]

def file_initiate(filename):
	fo=open(filename,'w')
	fo.close()

def Command_parametre_readin():
	global lenght_cff
	lenght_cff=100
	global dots_num_cff
	dots_num_cff=20
	global clu_dis_cff
	clu_dis_cff=5
	global point_dis_cff
	point_dis_cff=20
	global simple_dis_cff
	simple_dis_cff=100#min resolution:50bp
	global invert_base
	invert_base = { 'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C','N' : 'N','a' : 't', 't' : 'a', 'c' : 'g', 'g' : 'c','n' : 'n'}
	opts,args=getopt.getopt(sys.argv[2:],'i:',['sv-input=','reference=','pacbio-input=','output-path=','output-file='])
	dict_opts=dict(opts)
	global out_path
	out_path=path_modify(dict_opts['--output-path'])
	path_mkdir(out_path)
	global out_file_Cannot_Validate
	out_file_Cannot_Validate=out_path+'Unable.to.validate.rec'
	file_initiate(out_file_Cannot_Validate)
	global sample_name
	sample_name='.'.join(dict_opts['--sv-input'].split('/')[-1].split('.')[:-1])
	global start
	start=0
	global delta
	delta=50
	global bam_in
	bam_in=dict_opts['--pacbio-input']
	global ref
	ref=dict_opts['--reference']
	global chromos
	chromos=chromos_readin(ref)
	global flank_length
	flank_length=500
	global region_QC_Cff
	region_QC_Cff=0.4#QC score for repetitive regions; score=dots on diagnal / all dots when plotting ref vs. ref
	global min_length
	min_length=50
	global min_read_compare
	min_read_compare=20
	global case_number
	case_number=-1
	global qc_file_name
	qc_file_name=out_path+'reference_QC.rec'
	file_initiate(qc_file_name)

def svelter_readin(file_in):
	fin=open(file_in)
	header=fin.readline().strip().split()
	sample_list=header[3:]
	out_hash={}
	for line in fin:
		pin=line.strip().split()
		pin[4] = '_'.join(pin[4].split('/'))
		pin[5] = '_'.join(pin[5].split('/'))
		if not pin[4] in list(out_hash.keys()):	out_hash[pin[4]]={}
		if not pin[5] in list(out_hash[pin[4]].keys()):	out_hash[pin[4]][pin[5]]=[]
		if not pin[3].split(':') in out_hash[pin[4]][pin[5]]:	out_hash[pin[4]][pin[5]].append(pin[3].split(':'))
	fin.close()
	return out_hash

if len(sys.argv)<2: 
	import prep as prep
	prep.print_read_me()
else:
	import getopt,scipy,random,math
	from numpy import vstack,array
	from numpy.random import rand
	from scipy.cluster.vq import vq, kmeans, whiten
	from scipy import stats
	from scipy.stats import linregress
	from sklearn import cluster
	from scipy.spatial import distance
	import sklearn.datasets
	from sklearn.preprocessing import StandardScaler
	function_name=sys.argv[1]
	global dict_opts,num_reads_cff

	parser = argparse.ArgumentParser(
			description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('--sv-input', required=True, help='input file of SV calls')
	parser.add_argument('--reference', required=True, help='reference sequences')
	parser.add_argument('--pacbio-input', required=True, help='input pacbio sequences in bam format')
	parser.add_argument('--output-path', required=True, help='path of output VaPoR figures')
	parser.add_argument('--output-file', required=True, help = 'name of output file')
	parser.add_argument('--PB-supp', required=False, help='minimum number of evaluable PacBio reads')
	args = parser.parse_args(sys.argv[2:])

	if function_name=='ins':
		if len(sys.argv)==2:
			import prep as prep
			prep.readme_melt()
		else:
			def main():
				num_reads_cff=3
				if args.PB_supp:
					num_reads_cff = int(args.PB_supp)
				out_path=args.output_path
				path_mkdir(out_path)
				melt_input_prefix=args.sv_input_prefix
				sample_name=melt_input_prefix.split('/')[-1].split('.')[0]
				bam_in=args.pacbio_input
				ref=args.reference
				melt_info_readin(melt_input_prefix,out_path,sample_name,bam_in,ref)
			main()
	elif function_name=='bed':
		if len(sys.argv)==2:
			import prep as prep
			prep.readme_bed()
		else:
			#requirement for bed format: chr start end sv_type; format of sv_type: INS_280/ DEL_1_0/ DUP_<CN2>_<CN1>
			[out_name, out_path, bed_input, sample_name, bam_in, ref] = [args.output_file, args.output_path,args.sv_input,'.'.join(args.sv_input.split('/')[-1].split('.')[:-1]),args.pacbio_input,args.reference]

			num_reads_cff=3
			if args.PB_supp:
				num_reads_cff = int(args.PB_supp)

			out_path=path_modify(out_path)
			path_mkdir(out_path)
			bed_info=bed_info_readin(bed_input,out_path)
			write_output_initiate(out_name)
			plt_li=0
			for x in bed_info:
				if x[-1] in ['a/','/a','/','DEL']:
					key_event=':'.join([str(i) for i in x[:-3]]+['DEL'])
					plt_li+=1
					vapor_score_event=vapor_simple_del_Vapor(num_reads_cff,plt_li,bam_in,ref,x[:-3],out_path+sample_name+'.DEL.'+key_event.replace(':','__')+'.png')
					vapor_result=result_organize_ins([key_event,vapor_score_event])
					write_output_main(out_name,vapor_result[0].split(':')+[x[3]]+vapor_result[1:])
					print (vapor_result)
				elif x[-1] in ['a/a^','a^/a','a^/a^','INV']:	
					key_event=':'.join([str(i) for i in x[:-3]]+['INV'])
					plt_li+=1
					vapor_score_event=vapor_simple_inv_Vapor(num_reads_cff,plt_li,bam_in,ref,x[:-3],out_path+sample_name+'.INV.'+key_event.replace(':','__')+'.png')
					vapor_result=result_organize_ins([key_event,vapor_score_event])
					write_output_main(out_name,vapor_result[0].split(':')+[x[3]]+vapor_result[1:])
					print (vapor_result)
				elif x[-1] in ['INS']:
					key_event=':'.join([str(i) for i in x[:-3]+['INS']])
					plt_li+=1
					ins_pos='_'.join([str(i) for i in x[:2]])
					POLARITY='+'
					if type(x[4])==type(4):					ins_seq=''.join(['X' for i in range(x[4])])
					else:									ins_seq=x[4]

					vapor_score_event=vapor_simple_ins_Vapor(num_reads_cff,plt_li,bam_in,ref,ins_pos,ins_seq,out_path+sample_name+'.INS.'+key_event.replace(':','__')+'.png',POLARITY)
					vapor_result=result_organize_ins([key_event,vapor_score_event])
					write_output_main(out_name,vapor_result[0].split(':')+[x[3]]+vapor_result[1:])
					print (vapor_result)
				elif x[-1] in ['a/aa','aa/a','aa/aa','DUP','TANDUP']:
					key_event=':'.join([str(i) for i in x[:-3]]+['TANDUP'])
					plt_li+=1
					vapor_score_event=vapor_simple_tandup_Vapor(num_reads_cff,plt_li,bam_in,ref,x[:-3],out_path+sample_name+'.TANDUP.'+key_event.replace(':','__')+'.png')
					vapor_result=result_organize_ins([key_event,vapor_score_event])
					write_output_main(out_name,vapor_result[0].split(':')+[x[3]]+vapor_result[1:])
					print (vapor_result)
				else:	print (x)
	elif function_name=='vcf':
		if len(sys.argv)==2:
			import prep as prep
			prep.readme_vcf()
		else:
			#command='vapor.py vcf --sv-input /scratch/remills_flux/xuefzhao/SV_discovery_index/download/HG00512.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.vcf --output-path /scratch/remills_flux/xuefzhao/SV_discovery_index/download/HG00512.vcf.vapor.py/  --reference /scratch/remills_flux/xuefzhao/reference/GRCh38.1KGP/GRCh38_full_analysis_set_plus_decoy_hla.fa --pacbio-input /nfs/turbo/remillsscr/scratch2_trans/datasets/1000genomes/vol1/ftp/data_collections/hgsv_sv_discovery/PacBio/alignment/HG00512.XXX.bam'

			[ref,vcf_input,sample_name,out_path,bam_in] = [args.reference,args.sv_input,'.'.join(args.sv_input.split('/')[-1].split('.')[:-1]),path_modify(args.output_path),args.pacbio_input]

			num_reads_cff=3
			if args.PB_supp:
				num_reads_cff = int(args.PB_supp)

			path_mkdir(out_path)
			[vcf_list,vcf_rec_hash]=vcf_list_readin(vcf_input)
			vcf_rec_hash_new=vcf_rec_hash_modify(vcf_rec_hash)
			write_output_initiate(vcf_input+'.vapor.py')
			plt_li=0
			for x in list(vcf_list.keys()):
				if x=='DEL':
					for y in vcf_list[x]:
						if not 'NA' in y:
							print (y)
							#eg of y=['chr1', '91666237', '91667096']
							if y[2]-y[1]<50: 
								key_event=':'.join([str(i) for i in y]+['DEL'])
								plt_li+=1
								vapor_score_event=[]
								write_output_main(vcf_input+'.vapor.py',result_organize_ins([key_event,vapor_score_event]))
							else:
								key_event=':'.join([str(i) for i in y]+['DEL'])
								plt_li+=1
								vapor_score_event=vapor_simple_del_Vapor(num_reads_cff,plt_li,bam_in,ref,y,out_path+sample_name+'.DEL.'+key_event.replace(':','__')+'.png')
								write_output_main(vcf_input+'.vapor.py',result_organize_ins([key_event,vapor_score_event]))
				elif x=='INV':
					for y in vcf_list[x]:
						if not 'NA' in y:
							print (y)
							#eg of y=['chr1', '91666237', '91667096']
							if y[2]-y[1]<50: 
								key_event=':'.join([str(i) for i in y]+['DEL'])
								plt_li+=1
								vapor_score_event=[]
								write_output_main(vcf_input+'.vapor.py',result_organize_ins([key_event,vapor_score_event]))
							else:
								key_event=':'.join([str(i) for i in y]+['INV'])
								plt_li+=1
								vapor_score_event=vapor_simple_inv_Vapor(num_reads_cff,plt_li,bam_in,ref,y,out_path+sample_name+'.INV.'+key_event.replace(':','__')+'.png')
								write_output_main(vcf_input+'.vapor.py',result_organize_ins([key_event,vapor_score_event]))
				elif x=='INS':
					for y in vcf_list[x]:
						if not 'NA' in y:
							print (y)
							key_event=':'.join([str(i) for i in y[:3]+['INS']])
							ins_pos='_'.join([str(i) for i in y[:2]])
							POLARITY='+'
							if len(y)==4:			ins_seq=y[-1]
							else:					ins_seq=''.join(['X' for i in range(y[2])])
							plt_li+=1
							vapor_score_event=vapor_simple_ins_Vapor(num_reads_cff,plt_li,bam_in,ref,ins_pos,ins_seq,out_path+sample_name+'.INS.'+key_event.replace(':','__')+'.png',POLARITY)
							write_output_main(vcf_input+'.vapor.py',result_organize_ins([key_event,vapor_score_event]))
				elif x=='DISDUP':
					for y in vcf_list[x]:
						if not 'NA' in y:
							print (y)
							key_event=':'.join([str(i) for i in y+['DISDUP']])
							plt_li+=1
							vapor_score_event=vapor_simple_disdup_Vapor(num_reads_cff,plt_li,bam_in,ref,y,out_path+sample_name+'.DISDUP.'+key_event.replace(':','__')+'.png')
							write_output_main(vcf_input+'.vapor.py',result_organize_ins([key_event,vapor_score_event]))
				elif x=='DEL_INV':
					for y in vcf_list[x]:
						if not 'NA' in y:
							print (y)
							key_event=':'.join(['_'.join([str(i) for i in j]) for j in y]+['DEL_INV'])
							plt_li+=1
							vapor_score_event=vapor_del_inv_Vapor(num_reads_cff,plt_li,bam_in,ref,y,out_path+sample_name+'.DEL_INV.'+key_event.replace(':','__')+'.png')
							write_output_main(vcf_input+'.vapor.py',result_organize_ins([key_event,vapor_score_event]))
				elif x=='DUP_INV':
					for y in vcf_list[x]:
						if not 'NA' in y:
							print (y)
							#y=['chr7', '121522364', '121522527', 'chr7', '121522527']
							key_event=':'.join([str(i) for i in y+['DUP_INV']])
							plt_li+=1
							vapor_score_event=vapor_dup_inv_VapoR(num_reads_cff,plt_li,bam_in,ref,y,out_path+sample_name+'.DUP_INV.'+key_event.replace(':','__')+'.png')
							write_output_main(vcf_input+'.vapor.py',result_organize_ins([key_event,vapor_score_event]))
				elif x=='Other':
					for y in vcf_list[x]:
						if not 'NA' in y:
							print (y)
							#eg of y=['ab_ab', 'ab_bb^', 'chr15', '95148444', '95149304', '95149530']
							key_event=':'.join([str(i) for i in y+['CANNOT_CLASSIFY']])
							plt_li+=1
							vapor_score_event=vapor_CANNOT_CLASSIFY_VapoR(num_reads_cff,plt_li,bam_in,ref,y,out_path+sample_name+'.CANNOT_CLASSIFY.'+key_event.replace(':','__')+'.png')
							write_output_main(vcf_input+'.vapor.py',result_organize_ins([key_event,vapor_score_event]))
				else:
					print (x)
			vcf_vapor_modify(vcf_input,vcf_rec_hash_new)
	elif function_name=='svelter':
		#command='vapor.py svelter --sv-input /scratch/remills_flux/xuefzhao/SV_discovery_index/download/SVelter.version16/HG00512.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.svelter --output-path /scratch/remills_flux/xuefzhao/SV_discovery_index/download/SVelter.version16/HG00512.svelter.vapor.py/  --reference /scratch/remills_flux/xuefzhao/reference/GRCh38.1KGP/GRCh38_full_analysis_set_plus_decoy_hla.fa --pacbio-input /nfs/turbo/remillsscr/scratch2_trans/datasets/1000genomes/vol1/ftp/data_collections/hgsv_sv_discovery/PacBio/alignment/HG00512.XXX.bam'
		#command='vapor.py svelter --sv-input /scratch/remills_flux/xuefzhao/SV_discovery_index/download/HG00512.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.vcf --output-path /scratch/remills_flux/xuefzhao/SV_discovery_index/download/HG00512.vcf.vapor.py/  --reference /scratch/remills_flux/xuefzhao/reference/GRCh38.1KGP/GRCh38_full_analysis_set_plus_decoy_hla.fa --pacbio-input /nfs/turbo/remillsscr/scratch2_trans/datasets/1000genomes/vol1/ftp/data_collections/hgsv_sv_discovery/PacBio/alignment/HG00512.XXX.bam'

		[ref,svelter_input,sample_name,out_path,bam_in] = [args.reference,args.sv_input,'.'.join(args.sv_input.split('/')[-1].split('.')[:-1]),path_modify(args.output_path),args.pacbio_input]

		num_reads_cff=3
		if args.PB_supp:
			num_reads_cff = int(args.PB_supp)

		path_mkdir(out_path)
		svelter_hash=svelter_readin(svelter_input)
		#Command_parametre_readin()
		plt_li = 0
		for k1 in list(svelter_hash.keys()):
			for k2 in list(svelter_hash[k1].keys()):
				for k3 in svelter_hash[k1][k2]:
					plt_li+=1
					key_event = '.'+'_'.join(k3)
					out_figure_name = out_path+sample_name+key_event.replace(':','__')+'.png'
					#eg of sv_info=['ab_ab', 'b_b^', 'chr7', '70955990', '70961199', '70973901']
					sv_info = [k1,k2]+k3
					print(sv_info)
					vapor_score_event=vapor_CANNOT_CLASSIFY_VapoR(num_reads_cff,plt_li,bam_in,ref,sv_info,out_figure_name)
					write_output_main(args.output_file,result_organize_ins([key_event,vapor_score_event]))




