
input_file=$1
sed -i 's/\s/\t/g' $input_file

prefix=${input_file%%.txt}
prefix=${prefix##hnf4a_pwm_}

echo $prefix

python ~/pipelines/thermo/hnf4a_generate_ka.py $prefix

perl ~/pipelines/tf_newpwm_prob_perPos.pl -p $input_file -f $scr/sequencing/noRep_padded_pft_v1_name_bc.txt -o /scratch/bclab/yuning/cutTagModel/0523_redo_hnf4a/trainedPWM_hnf4a_perPos_score.txt -l n -r CAAAGTCCA

