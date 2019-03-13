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

fasta_location=/picb/rnasys/share/raw_data/single_cell_rna/mouse/2015_genemo_bi
#need to change in future, 
#TODO : need disk to store user's uploading data
samplelist=(
2-cell_embryo1-cell1
2-cell_embryo1-cell2
2-cell_embryo2-cell1
2-cell_embryo2-cell2
2-cell_embryo3-cell1
2-cell_embryo3-cell2
4-cell_1_a
4-cell_1_b
4-cell_2_a
4-cell_2_b
8-cell_1_a
8-cell_1_b
8-cell_2_a
oocyte_embryo1_cell1
oocyte_embryo2_cell1
oocyte_embryo3_cell1
oocyte_embryo4_cell1
oocyte_embryo5_cell1
zygote_embryo1_cell1
zygote_embryo2_cell1
zygote_embryo3_cell1
zygote_embryo4_cell1
zygote_embryo5_cell1
)
#*samplelist=(
#sample_name1
#)
project_name=mouse_2015_project
#project_name=test
#project name can be any string 
#runningsetting
thread_use_hisat2=30
thread_use_miso=30
read_len=101
thread_use_rsem=30
#thread_use_mapsplice=20
#thread_use_miso=20



hisat=/picb/rnasys/program/src/hisat2-2.1.0/hisat2
#summarize=/usr/local/bin/summarize_miso
#miso=/usr/local/bin/miso

#ref_annotation_pwd
hisat2_index_mouse=/picb/rnasys/share/database/Mus_musculus/UCSC/mm10/hisat2_index/genome

#gene_gtf=/data-2/reference/mm10/annotation/genes.gtf
const_exons_gff=/picb/rnasys/share/database/miso_annotation/ensGene.min_500.const_exons.gff
#resm_ref=/picb/rnasys/share/database/Mus_musculus/UCSC/mm10/rsem_reference/mouse_ref_genecode
resm_ref=/picb/rnasys/share/database/Mus_musculus/UCSC/mm10/rsem_anno/mm10
#parameter
eventslist=(A3SS A5SS MXE RI SE)



# working_pathway=`pwd`
mkdir ${working_pathway}/$project_name
strat_time = `date`
 for sample in ${samplelist[@]}
 do
  seq1=${fasta_location}/${sample}_1.fastq
  seq2=${fasta_location}/${sample}_2.fastq
   mkdir ${sample}

  $hisat -x $hisat2_index_mouse -1 $seq1 -2 $seq2 -p $thread_use_hisat2  -S ${working_pathway}/${sample}/out.sam
	  	samtools view -bS ${working_pathway}/${sample}/out.sam > ${working_pathway}/${sample}/out.bam
	 	samtools sort  -@30 ${working_pathway}/${sample}/out.bam -o ${working_pathway}/${sample}/out.sort.bam
		samtools index -@30 ${working_pathway}/${sample}/out.sort.bam 
 	  rsem-calculate-expression --paired-end --ci-memory 4096 --time --estimate-rspd -p $thread_use_rsem  $seq1 $seq2 $resm_ref ${working_pathway}/${sample}/${sample}
	 cp ${working_pathway}/${sample}/${sample}.genes.results ${working_pathway}/${project_name}/ &
done


#
for sample in ${samplelist[@]}
do
	echo 'go!'
	pe_utils --compute-insert-len ${working_pathway}/${sample}/out.sort.bam $const_exons_gff --output-dir ${working_pathway}/${sample}
	paired_end1=$(grep -o 'mean=\w*\.\w*' ${working_pathway}/${sample}/out.sort.bam.insert_len | sed 's/mean=//g')
	paired_end2=$(grep -o 'sdev=\w*\.\w*' ${working_pathway}/${sample}/out.sort.bam.insert_len | sed 's/sdev=//g')
	echo $paired_end1
	echo $paired_end2
	echo 'miso'
	for event in ${eventslist[@]}
	do
		mkdir ${sample}/out/${event}

		echo '3 is complete'
		indexed_event=/picb/rnasys/share/database/Mus_musculus/UCSC/mm10/miso_index/indexed_${event}
		miso --run $indexed_event ${working_pathway}/${sample}/out.sort.bam --output-dir ${working_pathway}/${sample}/out/${event} --read-len $read_len --paired-end $paired_end1 $paired_end2 -p $thread_use_miso
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



