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

fasta_location=/picb/rnasys/share/raw_data/single_cell_rna/mouse/2014GB
#need to change in future, 
#TODO : need disk to store user's uploading data
samplelist=(
Four-cell_Blastomere1_Embryo1
Four-cell_Blastomere1_Embryo2
Four-cell_Blastomere1_Embryo3
Four-cell_Blastomere1_Embryo4
Four-cell_Blastomere2_Embryo1
Four-cell_Blastomere2_Embryo2
Four-cell_Blastomere2_Embryo3
Four-cell_Blastomere2_Embryo4
Four-cell_Blastomere2_Embryo5
Four-cell_Blastomere3_Embryo1
Four-cell_Blastomere3_Embryo2
Four-cell_Blastomere3_Embryo3
Four-cell_Blastomere3_Embryo4
Four-cell_Blastomere3_Embryo5
Four-cell_Blastomere4_Embryo1
Four-cell_Blastomere4_Embryo2
Four-cell_Blastomere4_Embryo3
Four-cell_Blastomere4_Embryo4
Four-cell_Blastomere4_Embryo5
Inner-cell_mass_sample1
Inner-cell_mass_sample2
Inner-cell_mass_sample3
Inner-cell_mass_sample4
Trophectoderm_mass_sample1
Trophectoderm_mass_sample2
Two-cell_Blastomere1_Embryo10
Two-cell_Blastomere1_Embryo1
Two-cell_Blastomere1_Embryo2
Two-cell_Blastomere1_Embryo3
Two-cell_Blastomere1_Embryo4
Two-cell_Blastomere1_Embryo5
Two-cell_Blastomere1_Embryo6
Two-cell_Blastomere1_Embryo7
Two-cell_Blastomere1_Embryo8
Two-cell_Blastomere1_Embryo9
Two-cell_Blastomere2_Embryo10
Two-cell_Blastomere2_Embryo1
Two-cell_Blastomere2_Embryo2
Two-cell_Blastomere2_Embryo3
Two-cell_Blastomere2_Embryo4
Two-cell_Blastomere2_Embryo5
Two-cell_Blastomere2_Embryo6
Two-cell_Blastomere2_Embryo7
Two-cell_Blastomere2_Embryo8
Two-cell_Blastomere2_Embryo9
Zygote_cell1_Embryo1
Zygote_cell1_Embryo2
Zygote_cell1_Embryo3
Zygote_cell1_Embryo4
Zygote_cell1_Embryo5
Zygote_cell1_Embryo6
Zygote_cell1_Embryo7
Zygote_cell1_Embryo8
Zygote_cell1_Embryo9

)
#*samplelist=(
#sample_name1
#)
project_name=Mouse_2014
#project_name=test
#project name can be any string 
#runningsetting
thread_use_hisat2=60
thread_use_miso=60
read_len=100
thread_use_rsem=60
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



working_pathway=`pwd`
mkdir ${working_pathway}/$project_name
strat_time = `date`
 for sample in ${samplelist[@]}
 do
  seq1=${fasta_location}/${sample}.fastq
  mkdir ${sample}

  $hisat -x $hisat2_index_mouse -U $seq1  -p $thread_use_hisat2  -S ${working_pathway}/${sample}/out.sam
	  	samtools view -bS ${working_pathway}/${sample}/out.sam > ${working_pathway}/${sample}/out.bam
	 	samtools sort  -@30 ${working_pathway}/${sample}/out.bam -o ${working_pathway}/${sample}/out.sort.bam
		samtools index -@30 ${working_pathway}/${sample}/out.sort.bam 
 	  rsem-calculate-expression  --ci-memory 4096 --time --estimate-rspd -p $thread_use_rsem  $seq1  $resm_ref ${working_pathway}/${sample}/${sample}
	  cp ${working_pathway}/${sample}/${sample}.genes.results ${working_pathway}/${project_name}/ &
done


#
for sample in ${samplelist[@]}
do

	for event in ${eventslist[@]}
	do
		mkdir ${sample}/out/${event}

		echo '3 is complete'
		indexed_event=/picb/rnasys/share/database/Mus_musculus/UCSC/mm10/miso_index/indexed_${event}
		miso --run $indexed_event ${working_pathway}/${sample}/out.sort.bam --output-dir ${working_pathway}/${sample}/out/${event} --read-len $read_len  -p $thread_use_miso
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



