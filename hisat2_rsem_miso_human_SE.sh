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
read_len=100
fasta_location=/picb/rnasys/share/raw_data/single_cell_rna/human/nsmb2013
#need to change in future, 
#TODO : need disk to store user's uploading data
samplelist=(
2-cell_embryo1_Cell1
2-cell_embryo1_Cell2
2-cell_embryo2_Cell1
2-cell_embryo2_Cell2
2-cell_embryo3_Cell1
2-cell_embryo3_Cell2
4-cell_embryo1_Cell1
4-cell_embryo1_Cell2
4-cell_embryo1_Cell3
4-cell_embryo1_Cell4
4-cell_embryo2_Cell1
4-cell_embryo2_Cell2
4-cell_embryo2_Cell3
4-cell_embryo2_Cell4
4-cell_embryo3_Cell1
4-cell_embryo3_Cell2
4-cell_embryo3_Cell3
4-cell_embryo3_Cell4
8-cell_embryo1_Cell1
8-cell_embryo1_Cell2
8-cell_embryo1_Cell3
8-cell_embryo1_Cell4
8-cell_embryo2_Cell1
8-cell_embryo2_Cell2
8-cell_embryo2_Cell3
8-cell_embryo2_Cell4
8-cell_embryo2_Cell5
8-cell_embryo2_Cell6
8-cell_embryo2_Cell7
8-cell_embryo2_Cell8
8-cell_embryo3_Cell1
8-cell_embryo3_Cell2
8-cell_embryo3_Cell3
8-cell_embryo3_Cell4
8-cell_embryo3_Cell5
8-cell_embryo3_Cell7
8-cell_embryo3_Cell8
ES_passage1_C1
ES_passage1_C2
hESC_passage0_Cell10
hESC_passage0_Cell4
hESC_passage0_Cell5
hESC_passage0_Cell6
hESC_passage0_Cell7
hESC_passage0_Cell8
hESC_passage10_Cell10
hESC_passage10_Cell11
hESC_passage10_Cell12
hESC_passage10_Cell13
hESC_passage10_Cell14
hESC_passage10_Cell15
hESC_passage10_Cell16
hESC_passage10_Cell17
hESC_passage10_Cell18
hESC_passage10_Cell19
hESC_passage10_Cell1
hESC_passage10_Cell21
hESC_passage10_Cell22
hESC_passage10_Cell23
hESC_passage10_Cell24
hESC_passage10_Cell25
hESC_passage10_Cell26
hESC_passage10_Cell2
hESC_passage10_Cell3
hESC_passage10_Cell4
hESC_passage10_Cell5
hESC_passage10_Cell6
hESC_passage10_Cell7
hESC_passage10_Cell8
hESC_passage10_Cell9
Late-blastocyst_embryo1_Cell10
Late-blastocyst_embryo1_Cell11
Late-blastocyst_embryo1_Cell12
Late-blastocyst_embryo1_Cell1
Late-blastocyst_embryo1_Cell2
Late-blastocyst_embryo1_Cell3
Late-blastocyst_embryo1_Cell4
Late-blastocyst_embryo1_Cell5
Late-blastocyst_embryo1_Cell6
Late-blastocyst_embryo1_Cell7
Late-blastocyst_embryo1_Cell8
Late-blastocyst_embryo1_Cell9
Late-blastocyst_embryo2_Cell10
Late-blastocyst_embryo2_Cell1
Late-blastocyst_embryo2_Cell2
Late-blastocyst_embryo2_Cell3
Late-blastocyst_embryo2_Cell4
Late-blastocyst_embryo2_Cell5
Late-blastocyst_embryo2_Cell6
Late-blastocyst_embryo2_Cell7
Late-blastocyst_embryo2_Cell8
Late-blastocyst_embryo2_Cell9
Late-blastocyst_embryo3_Cell1
Late-blastocyst_embryo3_Cell2
Late-blastocyst_embryo3_Cell3
Late-blastocyst_embryo3_Cell4
Late-blastocyst_embryo3_Cell5
Late-blastocyst_embryo3_Cell6
Late-blastocyst_embryo3_Cell7
Late-blastocyst_embryo3_Cell8
Morulae_embryo1_Cell1
Morulae_embryo1_Cell2
Morulae_embryo1_Cell3
Morulae_embryo1_Cell4
Morulae_embryo1_Cell5
Morulae_embryo1_Cell6
Morulae_embryo1_Cell7
Morulae_embryo1_Cell8
Morulae_embryo2_Cell1
Morulae_embryo2_Cell2
Morulae_embryo2_Cell3
Morulae_embryo2_Cell4
Morulae_embryo2_Cell5
Morulae_embryo2_Cell6
Morulae_embryo2_Cell7
Morulae_embryo2_Cell8
Oocyte_embryo1_Cell1
Oocyte_embryo1_Cell2
Oocyte_embryo1_Cell3
Zygote_embryo1_Cell1
Zygote_embryo1_Cell2
Zygote_embryo1_Cell3

)
#*samplelist=(
#sample_name1
#)
project_name=human_2013
#project_name=test
#project name can be any string 
#runningsetting
thread_use_hisat2=30
thread_use_miso=30
thread_use_rsem=30



hisat=/picb/rnasys/program/src/hisat2-2.1.0/hisat2
#pe_utils=/picb/rnasys/program/install/miso/lib/python2.7/site-packages/misopy/pe_utils.py
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
   seq1=${fasta_location}/${sample}.fastq
  mkdir ${sample}
  mkdir ${sample}/insert-dist
  mkdir ${sample}/out
   $hisat -x $hisat2_index_human -U $seq1 -p $thread_use_hisat2  -S ${working_pathway}/${sample}/out.sam
	   samtools view -bS ${working_pathway}/${sample}/out.sam > ${working_pathway}/${sample}/out.bam
	 	samtools sort  -@30 ${working_pathway}/${sample}/out.bam -o ${working_pathway}/${sample}/out.sort.bam
		samtools index -@30 ${working_pathway}/${sample}/out.sort.bam 
	 rsem-calculate-expression  --ci-memory 4096 --time --estimate-rspd -p $thread_use_rsem   $seq1  $resm_ref ${working_pathway}/${sample}/${sample}
	 cp ${working_pathway}/${sample}/${sample}.genes.results ${working_pathway}/${project_name}/ &
done




for sample in ${samplelist[@]}
do
	
	for event in ${eventslist[@]}
	do
		mkdir ${sample}/out/${event}
		echo '3 is complete'
		indexed_event=${hg19_event}${event}
		miso --run $indexed_event ${working_pathway}/${sample}/out.sort.bam --output-dir ${working_pathway}/${sample}/out/${event} - --read-len $read_len   -p $thread_use_miso	
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


