set -e 


input_file=$1 ### pwm file
seq_file=$2
sed -i 's/\s/\t/g' $input_file

just_file=${input_file##*/}
tf_name=${just_file%%_pwm*}
#prefix=${input_file%%.txt}
prefix=${input_file##*_pwm_}

echo $prefix $tf_name

python ~/pipelines/thermo/generate_ka_randomseq.py $seq_file  $tf_name"_kas_"$prefix $input_file

if [[ $tf_name == foxa1 ]]
then
	echo "F prob"
	perl ~/pipelines/tf_newpwm_prob_perPos.pl -p $input_file -f $seq_file -o $tf_name"_prob_"$prefix -l n -r GTAAACAT

elif [[ $tf_name == hnf4a ]]
then
	echo "hnf4a prob"
	perl ~/pipelines/tf_newpwm_prob_perPos.pl -p $input_file -f $seq_file -o $tf_name"_prob_"$prefix -l n -r CAAAGTCCA
fi
