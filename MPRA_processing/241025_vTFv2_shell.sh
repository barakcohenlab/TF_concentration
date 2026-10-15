set -ex

dir=${1:-$PWD}

#ris ##run interactive session

#cd /scratch/bclab/jie/
#mkdir newData_8dox2
#cd /scratch/cgssb/miniseq1/Current_Runs/230912_MN00200_0714_A000H5LW7K/Data/Intensities/BaseCalls/manual/Cohen_YZ
#cp Cohen*_R* /scratch/bclab/jie/newData_8dox2/
#cd /scratch/bclab/jie/newData_8dox2
for i in `ls $dir/*R1*gz`; do gunzip -c $i > ${i%.gz}; done
for i in `ls $dir/*R1*q`; do j=${i%_mark_*};perl ~/pipelines/230803_spikein_pattern.pl $i $dir/"PM_"${j#*Cohen_};done
for i in `ls $dir/*all_ma*`; do awk -F '\t' 'NR==FNR{a[$1]=$3}NR>FNR{if($2 in a)print a[$2]"\t"$2"\t"$1}' /scratch/bclab/yuning/24_sequencing/bc_165cre_name_vTFV2.txt  <(grep hsp $i | awk -F '\t' '{print $3}' | sort | uniq -c | awk '{print $1"\t"$2}') > ${i%_all_matched.txt}"_hsp_barcode_count.txt"; done

##minP
for i in `ls $dir/*all_ma*`; do awk -F '\t' 'NR==FNR{a[$1]=$3}NR>FNR{if($2 in a)print a[$2]"\t"$2"\t"$1}' /scratch/bclab/yuning/24_sequencing/bc_165cre_name_vTFV2.txt  <(grep minP $i | awk -F '\t' '{print $3}' | sort | uniq -c | awk '{print $1"\t"$2}') > ${i%_all_matched.txt}"_minP_barcode_count.txt"; done


rm $dir/Cohen*fastq
rm $dir/*all_match*

#for i in `ls $dir/*count.txt`; do j=${i##*/}; perl ~/pipelines/231025_bc_merge.pl $i $dir/"1runC_"${j%_barcode_count.txt}; done 

