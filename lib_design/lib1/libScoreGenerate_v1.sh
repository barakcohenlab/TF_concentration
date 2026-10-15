#!/bin/bash

set -e

eval $(spack load --sh  bedtools2)

function sortbed(){
        sort -k1,1 -k2,2n -o $1 $1
}

probe_length_half=100

workdir=$1
tf_name=$2
up_down=$3
inputFile="/scratch/bclab/yuning/libDesign/pfPilot/workplace/geneMatch_fimo_"$tf_name"_"$up_down".txt" ## geneMatch_fimo_foxa1_up.txt
#numberofP=$3
numberofP="all"
outputFile1=$workdir"/probeSet_"$numberofP"_"$tf_name"_"$up_down".txt"
outputFile1_noSuf=${outputFile1%.txt}
outputFile1_root=${outputFile1_noSuf##*/}

echo $outputFile1_root

perl ~/pipelines/pfTF_libDesign_score.pl $inputFile $outputFile1_noSuf
echo "score probes done"

awk 'NR>1{print $NF,$0}' $outputFile1_noSuf"_score.txt"  | sort -nr | cut -f2- -d' ' > $outputFile1_noSuf"_sorted_score.txt"
#awk -F '\t' '{if(! ($1"_"$2 in a)){a[$1"_"$2]=$3;print $1"_"$2"\t"$3}}'   $outputFile1_noSuf"_sorted_score.txt" > $outputFile1_noSuf"_loci_gene_onemap.txt"
awk -F '\t' '{if(! ($1"_"$2 in a)){a[$1"_"$2]=$NF;b[$1"_"$2]=$1"_"$2"\t"$3}else{if(a[$1"_"$2] <= $NF){b[$1"_"$2]=b[$1"_"$2]","$3;a[$1"_"$2]=$NF}}}END{for(i in b)print b[i]}'   $outputFile1_noSuf"_sorted_score.txt" > $outputFile1_noSuf"_loci_gene_onemap.txt"

awk -F '\t' -vOFS='\t' 'NR>1{print $1,$2,$2+1,$NF,$1"_"$2"\t"$3}' $outputFile1_noSuf"_score.txt" > $outputFile1_noSuf"_libGenInput.bed"
sortbed $outputFile1_noSuf"_libGenInput.bed"
bedtools merge -d 10 -c 4,5,6 -o max,distinct,distinct -i $outputFile1_noSuf"_libGenInput.bed" > $outputFile1_noSuf"_afterMerge.bed"
sort -k4nr -k1 -k2n $outputFile1_noSuf"_afterMerge.bed" > $outputFile1_noSuf"_merged_sorted.bed"

echo "start generating"

#perl ~/pipelines/generateLib.pl $outputFile1_noSuf"_merged_sorted.bed" $numberofP $outputFile1_noSuf"_loci_gene_onemap.txt" $workdir
perl ~/pipelines/generateLib_noselection.pl $outputFile1_noSuf"_merged_sorted.bed" $outputFile1_noSuf"_loci_gene_onemap.txt" $workdir $probe_length_half

echo "generating done, start getting fasta"

bedtools getfasta -fi /ref/bclab/software/yuning/refs/hg19/latest/hg19.fa -bed $workdir/"coord_"$outputFile1_root".bed" -tab | awk -F '\t' '{print $1"\t"toupper($2)}' | awk -F ':|-' '{print $1"_"$2"_"$3"\t"$4}'> $workdir/"seq_"$outputFile1_root".txt"
echo "getfasta done"

cp $workdir/"coord_"$outputFile1_root".bed" $outputFile1_noSuf"_forpileup.bed"
sortbed $outputFile1_noSuf"_forpileup.bed"
# for i in `ls /scratch/bclab/yuning/libDesign/pfPilot/all_bdgs/*bdg`; 
# do 
#  	#bash ~/pipelines/wholePeak_pileup_submit_yz.sh -f 0 -p $outputFile1_noSuf"_forpileup.bed" -b $i -o $workdir
#  	python /home/yuningz/scripts/yz_pileup_onePk_oneBdg.py -c useWholePk -f 0 -p $outputFile1_noSuf"_forpileup.bed" -b $i -o $workdir;
# done




if [[ $inputFile == *foxa1* ]]
then
	if [[ $inputFile == *_up* ]]
	then
	        echo "match foxa1 up"
	        awk -F '\t' -vOFS='\t' 'NR==FNR{a[$1]=$1"\t"$3}NR>FNR{split($NF,b,",");for(i in b){if(b[i] in a)print $1,$2,$3,a[b[i]]}}' /scratch/bclab/yuning/libDesign/pfPilot/rna_expr_data/noNA_up_id_gene_fFC_hFC_fhFC.txt $workdir/"allinfo_"$outputFile1_root".bed"  > $workdir/"fcOneMatch_"$outputFile1_root".txt"
	        awk -F '\t' -vOFS='\t' 'NR==FNR{a[$1]=$1"\t"$3}NR>FNR{split($(NF-1),b,",");for(i in b){if(b[i] in a)print $1,$2,$3,a[b[i]]}}' /scratch/bclab/yuning/libDesign/pfPilot/rna_expr_data/noNA_up_id_gene_fFC_hFC_fhFC.txt $workdir/"allinfo_"$outputFile1_root".bed"  > $workdir/"fcmultiMatch_"$outputFile1_root".txt"

	elif [[ $inputFile == *_down* ]]
	then
	        echo "match foxa1 down"
	        awk -F '\t' -vOFS='\t' 'NR==FNR{a[$1]=$1"\t"$3}NR>FNR{split($NF,b,",");for(i in b){if(b[i] in a)print $1,$2,$3,a[b[i]]}}' /scratch/bclab/yuning/libDesign/pfPilot/rna_expr_data/noNA_down_id_gene_fFC_hFC_fhFC.txt $workdir/"allinfo_"$outputFile1_root".bed"  > $workdir/"fcOneMatch_"$outputFile1_root".txt"
	        awk -F '\t' -vOFS='\t' 'NR==FNR{a[$1]=$1"\t"$3}NR>FNR{split($(NF-1),b,",");for(i in b){if(b[i] in a)print $1,$2,$3,a[b[i]]}}' /scratch/bclab/yuning/libDesign/pfPilot/rna_expr_data/noNA_down_id_gene_fFC_hFC_fhFC.txt $workdir/"allinfo_"$outputFile1_root".bed"  > $workdir/"fcmultiMatch_"$outputFile1_root".txt"	        
	else
	        echo "no up down match, no done"
	fi
elif [[ $inputFile == *hnf4a* ]]
then
	if [[ $inputFile == *_up* ]]
	then
	        echo "match hnf4a up"
	        awk -F '\t' -vOFS='\t' 'NR==FNR{a[$1]=$1"\t"$4}NR>FNR{split($NF,b,",");for(i in b){if(b[i] in a)print $1,$2,$3,a[b[i]]}}' /scratch/bclab/yuning/libDesign/pfPilot/rna_expr_data/noNA_up_id_gene_fFC_hFC_fhFC.txt $workdir/"allinfo_"$outputFile1_root".bed"  > $workdir/"fcOneMatch_"$outputFile1_root".txt"
	        awk -F '\t' -vOFS='\t' 'NR==FNR{a[$1]=$1"\t"$4}NR>FNR{split($(NF-1),b,",");for(i in b){if(b[i] in a)print $1,$2,$3,a[b[i]]}}' /scratch/bclab/yuning/libDesign/pfPilot/rna_expr_data/noNA_up_id_gene_fFC_hFC_fhFC.txt $workdir/"allinfo_"$outputFile1_root".bed"  > $workdir/"fcmultiMatch_"$outputFile1_root".txt"
	elif [[ $inputFile == *_down* ]]
	then
	        echo "match hnf4a down"
	        awk -F '\t' -vOFS='\t' 'NR==FNR{a[$1]=$1"\t"$4}NR>FNR{split($NF,b,",");for(i in b){if(b[i] in a)print $1,$2,$3,a[b[i]]}}' /scratch/bclab/yuning/libDesign/pfPilot/rna_expr_data/noNA_down_id_gene_fFC_hFC_fhFC.txt $workdir/"allinfo_"$outputFile1_root".bed"  > $workdir/"fcOneMatch_"$outputFile1_root".txt"
	        awk -F '\t' -vOFS='\t' 'NR==FNR{a[$1]=$1"\t"$4}NR>FNR{split($(NF-1),b,",");for(i in b){if(b[i] in a)print $1,$2,$3,a[b[i]]}}' /scratch/bclab/yuning/libDesign/pfPilot/rna_expr_data/noNA_down_id_gene_fFC_hFC_fhFC.txt $workdir/"allinfo_"$outputFile1_root".bed"  > $workdir/"fcmultiMatch_"$outputFile1_root".txt"
	else
	        echo "no up down match, no done"
	fi

else
        echo "no F-H match, no done"
fi


echo "whole pipeline finished"

