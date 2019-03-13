#!/bin/bash
# ================================================================
# Orignaal File :                 pe_map_ase_hg19_rnasys_rnamen.sh
# Parameters    : read_len                                     101
#                 samplelist                          sample_name1
#                 project_name                                test
#                 thread_use_mapsplice                          20
#                 thread_use_miso                               20
# ================================================================

#input
#read_len=150
#read_len=101
fasta_location=/picb/rnasys/share/raw_data/single_cell_rna/human/NCB2018
read_len=150
#need to change in future, 
#TODO : need disk to store user's uploading data
samplelist=(
2cell_embryo10_1
2cell_embryo1_1
2cell_embryo2_1
2cell_embryo3_1
2cell_embryo4_1
2cell_embryo4_2
2cell_embryo5_1
2cell_embryo5_2
2cell_embryo6_1
2cell_embryo6_2
2cell_embryo7_1
2cell_embryo7_2
2cell_embryo8_1
2cell_embryo9_1
4cell_embryo10_1
4cell_embryo10_2
4cell_embryo11_1
4cell_embryo1_1
4cell_embryo11_2
4cell_embryo11_3
4cell_embryo1_2
4cell_embryo2_1
4cell_embryo2_2
4cell_embryo2_3
4cell_embryo3_1
4cell_embryo3_2
4cell_embryo3_3
4cell_embryo4_1
4cell_embryo4_2
4cell_embryo4_3
4cell_embryo5_1
4cell_embryo5_2
4cell_embryo6_1
4cell_embryo7_1
4cell_embryo7_2
4cell_embryo7_3
4cell_embryo8_1
4cell_embryo8_2
4cell_embryo9_1
4cell_embryo9_2
8cell_embryo10_1
8cell_embryo10_2
8cell_embryo10_3
8cell_embryo10_4
8cell_embryo10_5
8cell_embryo10_6
8cell_embryo11_1
8cell_embryo1_1
8cell_embryo11_2
8cell_embryo11_3
8cell_embryo12_1
8cell_embryo1_2
8cell_embryo1_3
8cell_embryo2_1
8cell_embryo2_2
8cell_embryo2_3
8cell_embryo2_4
8cell_embryo2_5
8cell_embryo3_1
8cell_embryo3_2
8cell_embryo3_3
8cell_embryo3_4
8cell_embryo3_5
8cell_embryo4_1
8cell_embryo4_2
8cell_embryo4_3
8cell_embryo4_4
8cell_embryo4_5
8cell_embryo5_1
8cell_embryo5_2
8cell_embryo5_3
8cell_embryo5_4
8cell_embryo5_5
8cell_embryo6_1
8cell_embryo6_2
8cell_embryo6_3
8cell_embryo6_4
8cell_embryo6_5
8cell_embryo6_6
8cell_embryo7_1
8cell_embryo7_2
8cell_embryo8_1
8cell_embryo8_2
8cell_embryo8_3
8cell_embryo9_1
8cell_embryo9_2
8cell_embryo9_3
8cell_embryo9_4
alpha-Amanitin_6cell_embryo1_1
alpha-Amanitin_7cell_embryo1_1
alpha-Amanitin_7cell_embryo1_2
alpha-Amanitin_8cell_embryo1_3
alpha-Amanitin_8cell_embryo1_4
alpha-Amanitin_8cell_embryo1_5
alpha-Amanitin_8cell_embryo1_6
alpha-Amanitin_8cell_embryo2_1
alpha-Amanitin_8cell_embryo2_2
alpha-Amanitin_8cell_embryo2_3
alpha-Amanitin_8cell_embryo2_4
alpha-Amanitin_8cell_embryo2_5
alpha-Amanitin_8cell_embryo2_6
alpha-Amanitin_8cell_embryo2_7
alpha-Amanitin_8cell_embryo2_8
alpha-Amanitin_8cell_embryo3_1
alpha-Amanitin_8cell_embryo3_2
alpha-Amanitin_8cell_embryo3_3
alpha-Amanitin_8cell_embryo3_4
alpha-Amanitin_8cell_embryo3_5
alpha-Amanitin_8cell_embryo3_6
ICM_embryo1_1
ICM_embryo1_2
ICM_embryo1_3
ICM_embryo1_4
ICM_embryo2_1
ICM_embryo2_2
ICM_embryo2_3
ICM_embryo2_4
ICM_embryo2_5
ICM_embryo3_1
ICM_embryo3_2
ICM_embryo3_3
ICM_embryo3_4
ICM_embryo3_5
ICM_embryo3_6
ICM_embryo3_7
ICM_embryo4_10
ICM_embryo4_11
ICM_embryo4_1
ICM_embryo4_12
ICM_embryo4_13
ICM_embryo4_2
ICM_embryo4_3
ICM_embryo4_4
ICM_embryo4_5
ICM_embryo4_6
ICM_embryo4_7
ICM_embryo4_8
ICM_embryo4_9
ICM_embryo5_1
ICM_embryo5_2
ICM_embryo5_3
ICM_embryo5_4
ICM_embryo5_5
ICM_embryo6_1
ICM_embryo6_2
ICM_embryo6_3
ICM_embryo6_4
ICM_embryo6_5
ICM_embryo7_1
ICM_embryo7_2
ICM_embryo7_3
ICM_embryo7_4
ICM_embryo7_5
ICM_embryo7_6
ICM_embryo7_7
ICM_embryo7_8
ICM_embryo7_9
ICM_embryo8_1
ICM_embryo8_2
Morula_embryo1_10
Morula_embryo1_11
Morula_embryo1_1
Morula_embryo1_2
Morula_embryo1_3
Morula_embryo1_4
Morula_embryo1_5
Morula_embryo1_6
Morula_embryo1_7
Morula_embryo1_8
Morula_embryo1_9
Morula_embryo2_10
Morula_embryo2_1
Morula_embryo2_2
Morula_embryo2_3
Morula_embryo2_4
Morula_embryo2_5
Morula_embryo2_6
Morula_embryo2_7
Morula_embryo2_8
Morula_embryo2_9
Morula_embryo3_1
Morula_embryo3_2
Morula_embryo3_3
Morula_embryo3_4
Oocyte_embryo1_1
Oocyte_embryo2_1
Oocyte_embryo3_1
Oocyte_embryo4_1
Oocyte_embryo5_1
Oocyte_embryo6_1
Oocyte_embryo7_1
Oocyte_embryo8_1
Oocyte_embryo9_1
Sperm_embryo10_1
Sperm_embryo11_1
Sperm_embryo1_1
Sperm_embryo12_1
Sperm_embryo13_1
Sperm_embryo2_1
Sperm_embryo3_1
Sperm_embryo4_1
Sperm_embryo5_1
Sperm_embryo6_1
Sperm_embryo7_1
Sperm_embryo8_1
Sperm_embryo9_1
TE_embryo1_1
TE_embryo1_2
TE_embryo2_10
TE_embryo2_11
TE_embryo2_1
TE_embryo2_12
TE_embryo2_13
TE_embryo2_14
TE_embryo2_15
TE_embryo2_16
TE_embryo2_2
TE_embryo2_3
TE_embryo2_4
TE_embryo2_5
TE_embryo2_6
TE_embryo2_7
TE_embryo2_8
TE_embryo2_9
TE_embryo3_1
TE_embryo3_2
TE_embryo3_3
TE_embryo3_4
TE_embryo4_10
TE_embryo4_11
TE_embryo4_1
TE_embryo4_12
TE_embryo4_13
TE_embryo4_2
TE_embryo4_3
TE_embryo4_4
TE_embryo4_5
TE_embryo4_6
TE_embryo4_7
TE_embryo4_8
TE_embryo4_9
TE_embryo5_1
TE_embryo5_12
TE_embryo5_13
TE_embryo5_2
Zygote_embryo10_1
Zygote_embryo1_1
Zygote_embryo2_1
Zygote_embryo3_1
Zygote_embryo4_1
Zygote_embryo5_1
Zygote_embryo6_1
Zygote_embryo7_1
Zygote_embryo8_1
Zygote_embryo9_1

)
#*samplelist=(
#sample_name1
#)
project_name=human_2018
#project_name=test
#project name can be any string 
#runningsetting
thread_use_hisat2=60
thread_use_miso=60
thread_use_rsem=60



hisat=/picb/rnasys/program/src/hisat2-2.1.0/hisat2
pe_utils=/picb/rnasys/program/install/miso/lib/python2.7/site-packages/misopy/pe_utils.py
const_exons_gff=/picb/rnasys/share/database/miso_annotation/Homo_sapiens.GRCh37.65.with_chr.min_1000.const_exons.gff
#summarize=/usr/local/bin/summarize_miso
#miso=/usr/local/bin/miso

#ref_annotation_pwd
hisat2_index_mouse=/picb/rnasys/share/database/Mus_musculus/UCSC/mm10/hisat2_index/genome
hisat2_index_human=/picb/rnasys/share/database/Homo_sapiens/UCSC/hg19/hisat_index/genome

#gene_gtf=/data-2/reference/mm10/annotation/genes.gtf
#const_exons_gff=/data-2/reference/mm10/annotation/ensGene.min_500.const_exons.gff

#resm_ref=/picb/rnasys/share/database/Mus_musculus/UCSC/mm10/rsem_reference/mouse_ref_genecode
#resm_ref=/picb/rnasys/share/database/Mus_musculus/UCSC/mm10/rsem_anno/mm10
resm_ref=/picb/rnasys/share/database/Homo_sapiens/UCSC/hg19/rsem_index/hg19_genecode
#parameter

##miso_indexed
mm10_event=/picb/rnasys/share/database/miso_annotation/miso_index_mm10/indexed_
hg19_event=/picb/rnasys/share/database/miso_annotation/miso_index_hg19/indexed_
eventslist=(A3SS_events A5SS_events MXE_events RI_events SE_events)


working_pathway=`pwd`
mkdir ${working_pathway}/$project_name
for sample in ${samplelist[@]}
do
   seq1=${fasta_location}/${sample}_1.fastq
   seq2=${fasta_location}/${sample}_2.fastq
  mkdir ${sample}
  mkdir ${sample}/insert-dist
  mkdir ${sample}/out
   $hisat -x $hisat2_index_human -1 $seq1 -2 $seq2 -p $thread_use_hisat2  -S ${working_pathway}/${sample}/out.sam
	   samtools view -bS ${working_pathway}/${sample}/out.sam > ${working_pathway}/${sample}/out.bam
	 	samtools sort  -@30 ${working_pathway}/${sample}/out.bam -o ${working_pathway}/${sample}/out.sort.bam
		samtools index -@30 ${working_pathway}/${sample}/out.sort.bam 
	 rsem-calculate-expression --paired-end --ci-memory 4096 --time --estimate-rspd -p $thread_use_rsem   $seq1 $seq2 $resm_ref ${working_pathway}/${sample}/${sample}
	 cp ${working_pathway}/${sample}/${sample}.genes.results ${working_pathway}/${project_name}/ &
done




for sample in ${samplelist[@]}
do
		pe_utils --compute-insert-len ${working_pathway}/${sample}/out.sort.bam $const_exons_gff --output-dir ${working_pathway}/${sample}/insert-dist
		paired_end1=$(grep -o 'mean=\w*\.\w*' ${working_pathway}/${sample}/insert-dist/out.sort.bam.insert_len | sed 's/mean=//g')
		paired_end2=$(grep -o 'sdev=\w*\.\w*' ${working_pathway}/${sample}/insert-dist/out.sort.bam.insert_len | sed 's/sdev=//g')
	for event in ${eventslist[@]}
	do
		mkdir ${sample}/out/${event}
		echo '3 is complete'
		indexed_event=${hg19_event}${event}
		miso --run $indexed_event ${working_pathway}/${sample}/out.sort.bam --output-dir ${working_pathway}/${sample}/out/${event} - --read-len $read_len  --paired-end $paired_end1 $paired_end2 -p $thread_use_miso	
		#miso --run $indexed_event ${working_pathway}/${sample}/out.sort.bam --output-dir ${working_pathway}/${sample}/out/${event} --read-len $read_len -p $thread_use_miso
		echo '4 is complete'
	#summarize-5
		summarize_miso --summarize-samples ${working_pathway}/${sample}/out/${event} ${working_pathway}/${sample}/out 
		echo '5 is complete'
		event_name=${event}.miso_summary
		cp ${working_pathway}/${sample}/out/summary/$event_name ${working_pathway}/${project_name}/${sample}_${event_name}
	done
    cd .. 
    echo 'finished'
done


