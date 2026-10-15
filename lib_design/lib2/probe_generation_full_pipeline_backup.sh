#!/bin/bash

set -e

outprefix="vtfV2"

### F

cd /scratch/bclab/yuning/libDesign/tfValidate/allinfo_files

perl ~/pipelines/tfValidate_design/ff_grammar.pl ../eachCRE_tfbs_list/f_a3_chr14_94862355_94862555.txt
perl ~/pipelines/tfValidate_design/ff_grammar.pl ../eachCRE_tfbs_list/f_a5_L1_random1.txt
perl ~/pipelines/tfValidate_design/ff_grammar.pl ../eachCRE_tfbs_list/f_a7_chr3_170720287_170720487.txt 

### H

perl ~/pipelines/tfValidate_design/hh_grammar.pl ../eachCRE_tfbs_list/h_b1_chr6_131148484_131148684.txt 
perl ~/pipelines/tfValidate_design/hh_grammar.pl ../eachCRE_tfbs_list/h_b2_chr9_104089765_104089965.txt
perl ~/pipelines/tfValidate_design/hh_grammar.pl ../eachCRE_tfbs_list/h_b3_L1_2highH.txt

### FH

perl ~/pipelines/tfValidate_design/fh_grammar.pl ../eachCRE_tfbs_list/bfh_c1_L1_highFH.txt
perl ~/pipelines/tfValidate_design/fh_grammar.pl ../eachCRE_tfbs_list/bfh_c2_chr2_231915124_231915324.txt
perl ~/pipelines/tfValidate_design/fh_grammar.pl ../eachCRE_tfbs_list/bfh_c3_chr1_57433643_57433843.txt 


### MISC
awk -F '\t' '{print $5"\tspikein_"$1}' /scratch/bclab/yuning/libDesign/pfPilot/workplace/v1_200/pft_v1/mpra_anchors_165_candidates_allInfo.txt > pftV1_mpra_anchors_165_candidates_allInfo.txt

for i in `ls ../eachCRE_tfbs_list/all_curated/cur*`; do echo $i; perl ~/pipelines/tfValidate_design/read_meme_motifs.pl $i; done

perl ~/pipelines/tfValidate_design/scramble_probe.pl

perl ~/pipelines/tfValidate_design/screen_TF_coop.pl 


cat *allInfo.txt > "tmp_cat_"$outprefix"_allInfo_tmp.txt"

perl ~/pipelines/tfValidate_design/cleanup.pl "tmp_cat_"$outprefix"_allInfo_tmp.txt"

### check if scramble is messed up

probe_num=$(cat "clean_cat_"$outprefix"_allInfo_all.txt" | wc -l)
#echo $probe_num
if [ $probe_num == 11633 ] 
then
	echo "probe number 11633 is right" 
else
	echo "probe number "$probe_num" is wrong"
	exit 1
fi


perl ~/pipelines/tfValidate_design/build_lib_v2.pl "clean_cat_"$outprefix"_allInfo_all.txt"


##### QC
outfilename=$outprefix"_agilent_submission.txt"

echo "file length:"
wc -l $outfilename ${outfilename%_agilent_submission.txt}"_agilent_alldata.txt" 

echo "barcode count:"
awk '{print substr($1,203,10)}' $outfilename | sort | uniq | wc -l

echo "number of uniq sequence:"
shortname=$outfilename
#sort $shortname | uniq | wc -l
awk '{print substr($1,1,200)}'  $shortname | sort | uniq | wc -l

echo "each seq final length:"
awk '{print length($0)}' $shortname | uniq


echo "restriction double check"
awk '{print gsub(/GCTAGC/, "")}' $shortname | uniq
awk '{print gsub(/GGTACCTGACCTC/, "")}' $shortname | uniq
awk '{print gsub(/AAGCTT/, "")}' $shortname | uniq
awk '{print gsub(/CTCGAG/, "")}' $shortname | uniq
awk '{print gsub(/GCATGC/, "")}' $shortname | uniq
awk '{print gsub(/GAGCTC/, "")}' $shortname | uniq
awk '{print gsub(/TCGCTTCGAGCA/, "")}' $shortname | uniq








echo "########## successfully reach the end ##########"
