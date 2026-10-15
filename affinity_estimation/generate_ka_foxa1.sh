
input_file=$1
sed -i 's/\s/\t/g' $input_file

prefix=${input_file%%.txt}
prefix=${prefix##foxa1_pwm_}

echo $prefix

python ~/pipelines/thermo/generate_ka.py $prefix

perl ~/pipelines/tf_newpwm_prob_perPos.pl -p $input_file -f $scr/sequencing/noRep_padded_pft_v1_name_bc.txt -o /scratch/bclab/yuning/cutTagModel/0523_redo_foxa1/trainedPWM_foxa1_perPos_score.txt -l n -r GTAAACAT

