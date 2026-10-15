#!/usr/bin/perl
#

use Getopt::Long;
use List::Util qw[max sum];
use warnings;


my $cre_motif_list = $ARGV[0];
my ($output_prefix_name) =  $cre_motif_list =~ m/f_(a[0-9]_.*).txt/;
my ($output_prefix) =  $output_prefix_name =~ m/a[0-9]_(.*)$/;


my $motif_file = "/ref/bclab/software/yuning/refs/jaspar_2024/human755_20240529015908_JASPAR2024_combined_matrices_615921_meme.txt";

open(REF1,"<$motif_file") || die "no $motif_file\n";

my %motif_hash;
my $tmp_pos_count = 0;
my $tmp_tf;

while(my $line = <REF1>){
	chomp $line;
	if($line =~ /^MOTIF\s(.*)\s/){
		$tmp_tf = $1;
		#print $tmp_tf."\n";
		$tmp_pos_count = 0;
	}
	if($line =~ /^\s/){
		my @row=split(" ",$line);
		$motif_hash{$tmp_tf}{"A"}[$tmp_pos_count] = $row[0];
		$motif_hash{$tmp_tf}{"C"}[$tmp_pos_count] = $row[1];
		$motif_hash{$tmp_tf}{"G"}[$tmp_pos_count] = $row[2];
		$motif_hash{$tmp_tf}{"T"}[$tmp_pos_count] = $row[3];
		$motif_hash{$tmp_tf}{"a"}[$tmp_pos_count] = $row[0];
		$motif_hash{$tmp_tf}{"c"}[$tmp_pos_count] = $row[1];
		$motif_hash{$tmp_tf}{"g"}[$tmp_pos_count] = $row[2];
		$motif_hash{$tmp_tf}{"t"}[$tmp_pos_count] = $row[3];

		$tmp_pos_count++;
		#print $motif_hash{$tmp_tf}[$tmp_pos_count][3]."\n";

	}
	#if($line =~ /^URL/){
		#$motif_hash{$tmp_tf} = \@tmp_matrix;
		#}
}


close(REF1);


my $sequence_file = "/scratch/bclab/yuning/sequencing/noRep_padded_pft_v1_name_bc.txt";
my %seq_hash;
open(REF3,"<$sequence_file") || die "no $sequence_file\n";
while(my $line = <REF3>){
	chomp $line;
	my @row=split("\t",$line);
	$seq_hash{$row[0]} = $row[1];
	#print $row[0]."\t".$row[1]."\n";
	#print @{$motif_hash{$row[0]}{"A"}};
}
close(REF3);


my %motif_dict;
my $ori_full_seq;
my $tfbs_number = 0;


open(REF2,"<$cre_motif_list") || die "no $cre_motif_list\n";
while(my $line = <REF2>){
	chomp $line;
	my @row=split("\t",$line);
	my $offset = $row[3]-1;
	$ori_full_seq = substr($seq_hash{$row[2]},8,165);
	#print "$ori_full_seq\n";
	if ($row[5] eq '+') {
		#print substr($seq_hash{$row[2]},$offset,$row[4]-$offset)." aaaaa\n";
	}else{
		#print rc(substr($seq_hash{$row[2]},$offset,$row[4]-$offset))." bbbbb\n";
	}
	my $tf_here = $row[0];
	my @ori_motif = split("",$row[9]);  ### motif here is the motif orientation
	my @best_motif = @ori_motif;
	my @mut_motif = @ori_motif;
	my @reallymut_motif = @ori_motif;

	$motif_dict{$tfbs_number}{"ss"} = $offset-8;
	$motif_dict{$tfbs_number}{"len"} = $row[4]-$offset;
	#$motif_dict{$tfbs_number}{"ori"} = $row[9];
	$motif_dict{$tfbs_number}{"ori"} = substr($seq_hash{$row[2]},$offset,$row[4]-$offset); ### ori here is the probe orientation
	$motif_dict{$tfbs_number}{"strand"} = $row[5];
	$motif_dict{$tfbs_number}{"tf"} = $row[1];
	$motif_dict{$tfbs_number}{"tfshort"} = substr($row[1], 0, 1);

	### Make consensus motif
	foreach my $i (0..$#ori_motif){
		my ($score, $base) = max_pos($tf_here,$i);
		if ($score > 0.5) {
			$best_motif[$i] = $base;
		}
	}
	if ($row[5] eq '+') {
		$motif_dict{$tfbs_number}{"best"} = join("",@best_motif);
	}else{
		$motif_dict{$tfbs_number}{"best"} = rc(join("",@best_motif));
	}


	### Make dead motif
	my $fold_change = 1;
	while($fold_change>0.001){
		my ($ori_s, $mut_pos) = return_maxPWM_loci($tf_here, join("",@mut_motif));
		my ($mut_s, $mut_base) = min_pos($tf_here, $mut_pos);
		$mut_motif[$mut_pos] = $mut_base;
		$fold_change = ($mut_s+1e-2)/($ori_s+1e-2) * $fold_change;
		#print "$tf_here $row[9] ".join("",@mut_motif)." $mut_pos $ori_s $mut_base $mut_s $fold_change\n";
	}
	if ($row[5] eq '+') {
		$motif_dict{$tfbs_number}{"mut"} = join("",@mut_motif);
	}else{
		$motif_dict{$tfbs_number}{"mut"} = rc(join("",@mut_motif));
	}


	### Make really dead motif (3 point mutations)
	my $mutation_count = 1;
	my $fold_change2 = 1;
	while($mutation_count<4){
		my ($ori_s, $mut_pos) = return_maxPWM_loci($tf_here, join("",@reallymut_motif));
		my ($mut_s, $mut_base) = min_pos($tf_here, $mut_pos);
		$reallymut_motif[$mut_pos] = $mut_base;
		$fold_change2 = ($mut_s+1e-2)/($ori_s+1e-2) * $fold_change2;
		$mutation_count++;
		#print "$tf_here $row[9] ".join("",@reallymut_motif)." $mut_pos $ori_s $mut_base $mut_s $mutation_count $fold_change2\n";
	}
	if ($row[5] eq '+') {
		$motif_dict{$tfbs_number}{"rmut"} = join("",@reallymut_motif);
	}else{
		$motif_dict{$tfbs_number}{"rmut"} = rc(join("",@reallymut_motif));
	}


	#print "$row[9]\t".join("",@best_motif)."\n";
	#print "$row[9]\t$motif_dict{$tfbs_number}{'ori'}\t".$motif_dict{$tfbs_number}{"mut"}."\t".$motif_dict{$tfbs_number}{"rmut"}."\n";
	#print $row[0]."\t".$row[1]."\t".consensus($row[0])."\t".$row[5]."\n";
	#print @{$motif_hash{$row[0]}{"A"}};
	$tfbs_number = $tfbs_number+1;

}

close(REF2);

my @tfbs_list = sort { $motif_dict{$a}{"ss"} <=> $motif_dict{$b}{"ss"} } keys(%motif_dict);
#print $tfbs_number."\n";
#print join(",", @tfbs_list)."\n";

### perl -le '{$a{"00"}="A";$a{"01"}="T";$a{"10"}="C";$a{"11"}="G";$n=5;for(0..(2**$n-1)){print join"::",map{$_} unpack("a1"x$n,sprintf("%0*b",$n,$_))}}' 

print "$output_prefix_name    $output_prefix\n";
open(OUT3, ">ff_$output_prefix_name"."_allInfo.txt");

### mutate/concensus/rc current BS

foreach my $n (0..(2**$tfbs_number-1)){
	my @curr_combination = unpack("a1"x$tfbs_number,sprintf("%0*b",$tfbs_number,$n));
	#print join("::",@curr_combination)."\n";
	my $edit_full_best = $ori_full_seq;
	my $edit_full_mut = $ori_full_seq;
	my $edit_full_rc = $ori_full_seq;

	my $seq_name = $output_prefix."_$motif_dict{0}{'tfshort'}_".join("",@curr_combination);

	foreach my $count (0..$#curr_combination){
		if ($curr_combination[$count] == 1) {
			#$seq_name = $seq_name.":".$motif_dict{$tfbs_list[$count]}{"tfshort"}.($count+1);

			substr($edit_full_mut, $motif_dict{$tfbs_list[$count]}{"ss"}, $motif_dict{$tfbs_list[$count]}{"len"}, lc($motif_dict{$tfbs_list[$count]}{"rmut"})); ### mut usually has two point mutations; still has some occupancy
			substr($edit_full_best, $motif_dict{$tfbs_list[$count]}{"ss"}, $motif_dict{$tfbs_list[$count]}{"len"}, lc($motif_dict{$tfbs_list[$count]}{"best"}));
			substr($edit_full_rc, $motif_dict{$tfbs_list[$count]}{"ss"}, $motif_dict{$tfbs_list[$count]}{"len"}, lc(rc($motif_dict{$tfbs_list[$count]}{"ori"})));
		}
	}
	print OUT3 "$edit_full_mut\tmut_$seq_name\n"; ### sometimes included in affinity
	print OUT3 "$edit_full_best\tcon_$seq_name\n"; ### sometimes included in affinity
	print OUT3 "$edit_full_rc\trc_$seq_name\n";
}


my $edit_mutAllBS = $ori_full_seq;

foreach my $tf (keys %motif_dict){
	substr($edit_mutAllBS, $motif_dict{$tf}{"ss"}, $motif_dict{$tf}{"len"}, ($motif_dict{$tf}{"rmut"}));
	#foreach my $j (sort keys %{$motif_dict{$tf}}){
		#print "$tf $j $motif_dict{$tf}{$j}\n";
	#}
}

#print "$ori_full_seq\n";
#print "$edit_deleteBS\n";
#print "$edit_mutAllBS\n";

print OUT3 "$ori_full_seq\toriFF_$output_prefix\n";
print OUT3 "$edit_mutAllBS\tmutAllFF_$output_prefix\n";

### move F BS around

my $tf1_bs_length = $motif_dict{0}{"len"};
my $bs_to_move = $motif_dict{0}{"best"};
my $tf_name = $motif_dict{0}{"tfshort"};

my $leftmost = 11-1;
my $rightmost = 156-1-$tf1_bs_length;

foreach my $loc ($leftmost..$rightmost){
	my $tmp_seq = $edit_mutAllBS;
	substr($tmp_seq, $loc, $tf1_bs_length, lc($bs_to_move));
	print OUT3 $tmp_seq."\t$tf_name"."move$loc"."_$output_prefix\n";
	### also move RC motif around?
	substr($tmp_seq, $loc, $tf1_bs_length, lc(rc($bs_to_move)));
	print OUT3 $tmp_seq."\t$tf_name"."moveRC$loc"."_$output_prefix\n";
}


### generate 2-way spacing

my $leftanchor = 11-1;
my $rightanchor = 135-1-8;

my $right1 = $edit_mutAllBS;
substr($right1, $rightanchor, $tf1_bs_length, lc($bs_to_move));
print OUT3 "$right1\t$tf_name"."sraNA_$output_prefix\n";

my $left1 = $edit_mutAllBS;
substr($left1, $leftanchor, $tf1_bs_length, lc($bs_to_move));
print OUT3 "$left1\t$tf_name"."slaNA_$output_prefix\n";

my $leftmost2 = $leftanchor+1; # +$tf1_bs_length
my $rightmost2 = $rightanchor-1; # -$tf1_bs_length


foreach my $loc ($leftmost2..$rightmost2){
	my $tmp_seq_right = $right1;
	my $tmp_seq_left = $left1;

	substr($tmp_seq_right, $loc, $tf1_bs_length, lc($bs_to_move));
	print OUT3 "$tmp_seq_right\t$tf_name"."sra$loc"."o1_$output_prefix\n";
	substr($tmp_seq_right, $loc, $tf1_bs_length, lc(rc($bs_to_move)));
	print OUT3 "$tmp_seq_right\t$tf_name"."sra$loc"."o2_$output_prefix\n";

	substr($tmp_seq_left, $loc, $tf1_bs_length, lc($bs_to_move));
	print OUT3 "$tmp_seq_left\t$tf_name"."sla$loc"."o1_$output_prefix\n";
	substr($tmp_seq_left, $loc, $tf1_bs_length, lc(rc($bs_to_move)));
	print OUT3 "$tmp_seq_left\t$tf_name"."sla$loc"."o2_$output_prefix\n";

}


### generate affinity combinations


### perl -le '{$a{"000"}=1;$a{"001"}=2;$a{"010"}=3;$a{"011"}=4;$a{"100"}=5;$n=9;for(0..(2**$n-1)){print join"::",map{$a{$_}} unpack("a3"x($n/3),sprintf("%0*b",$n,$_))}}'

# list_pick <- c("ATGTTTAC","ATATTTAC","ATGTGTAC","ATGTTGGC","ACGGTAAC") #acggtaac


my %affinity_hash_f ;  ### or read from file
$affinity_hash_f{"000"}= "ATGTTTAC"; # "---1----"
$affinity_hash_f{"001"}= "ATATTTAC"; # "---2----"
$affinity_hash_f{"010"}= "ATGTGTAC"; # "---3----"
$affinity_hash_f{"011"}= "ATGTTGGC"; # "---4----"
$affinity_hash_f{"100"}= "ACGGTAAC"; # "---5----"

my %affinity_hash_f_name = ("000"=>1, "001"=>2, "010"=>3, "011"=>4, "100"=>5 ) ;  ### or read from file

my $digit_c = 3;
my $total_digit = $tfbs_number*$digit_c;

foreach my $n (0..(2**($total_digit)-1)){
	my @curr_combination = unpack("a3"x$tfbs_number,sprintf("%0*b",$total_digit,$n));
	my $affinity_tmp = $edit_mutAllBS;

	my $seq_name = $output_prefix;

	my $valid_check = -1;
	foreach my $count (0..$#curr_combination){
		if (! (exists($affinity_hash_f{$curr_combination[$count]}) )) {
			$valid_check = 11111;
			last;
		}
	}

	if ($valid_check>0) {
		next;
	}
	#print join("::",@curr_combination)."\n";
	$seq_name = $seq_name."_".$tf_name;
	foreach my $count (0..$#curr_combination){
		$seq_name = $seq_name.$affinity_hash_f_name{$curr_combination[$count]};
		if ($motif_dict{$tfbs_list[$count]}{"strand"} eq "-") {
			substr($affinity_tmp, $motif_dict{$tfbs_list[$count]}{"ss"}, $motif_dict{$tfbs_list[$count]}{"len"}, lc($affinity_hash_f{$curr_combination[$count]}));
		}else{
			substr($affinity_tmp, $motif_dict{$tfbs_list[$count]}{"ss"}, $motif_dict{$tfbs_list[$count]}{"len"}, lc(rc($affinity_hash_f{$curr_combination[$count]})));
		}
	}
	print OUT3 "$affinity_tmp\taff_$seq_name\n";

}



### one BS with continuous affinity

my $center_loc = 87 - 1 - 8;

### use ATNTNTWC

my $continuous_kmer = "/scratch/bclab/yuning/libDesign/tfValidate/eachCRE_tfbs_list/ATNTNTWC.txt";
open(KMER,"<$continuous_kmer") || die "no $continuous_kmer\n";

my @ATNTNTWC;
while(my $line = <KMER>){
	chomp $line;
	my @row=split("\t",$line);
	push(@ATNTNTWC, $row[0]);
}
close(KMER);

my $tmp_con_aff_start = $ori_full_seq; 
my $tmp_con_aff_start2 = $edit_mutAllBS; 

foreach my $kmer (@ATNTNTWC){
	substr($tmp_con_aff_start, $center_loc, $tf1_bs_length, lc(($kmer)));
	print OUT3 "$tmp_con_aff_start\t$tf_name"."caffo$kmer"."_$output_prefix\n";
	substr($tmp_con_aff_start2, $center_loc, $tf1_bs_length, lc(($kmer)));
	print OUT3 "$tmp_con_aff_start2\t$tf_name"."caffm$kmer"."_$output_prefix\n";
}


### plug in more FOXA1 BS (use beta 25bp grammar)

my @plugin_sites = (62,87,112);
#print @plugin_sites;

$count_2 = scalar @plugin_sites;
#print $count_2."\n";

foreach my $n (0..(2**$count_2-1)){
	my @curr_combination = unpack("a1"x$count_2,sprintf("%0*b",$count_2,$n));
	#print join("::",@curr_combination)."\n";
	my $tmp_random_plugin = $ori_full_seq;

	my $seq_name = $output_prefix."_rpluginF";

	foreach my $count (0..$#curr_combination){
		if ($curr_combination[$count] == 1) {
			$seq_name = $seq_name.":".($count+1);

			substr($tmp_random_plugin, $plugin_sites[$count]-9, $tf1_bs_length, lc($bs_to_move)); ### mut usually has two point mutations; still has some occupancy
		}
	}
	print OUT3 "$tmp_random_plugin\t$seq_name\n";
}





############################################################################################################################

#foreach my $tf (keys %motif_hash){
	#print $tf."\t".scalar(@{$motif_hash{$tf}{"A"}})."\n";
	#print $tf."\t".consensus($tf)."\n";
	#foreach my $row (0..(scalar @{$motif_hash{$tf}{"A"}} - 1)){
		#print $row.",".$motif_hash{$tf}{"c"}[$row]."\t";
	#}
	#print "$tf\n";
#}

### print consensus and mutated TFBS;

#my ($t1, $t2) = max_pos("MA1142.2",3);
#print $t1.$t2;

sub return_maxPWM_loci{
	my ($tf, $sequence)= @_;
	my @individual_bases = split(//, $sequence);
	my $score =0;
	my $position=-1;
	foreach my $i (0..$#individual_bases){
		if ($motif_hash{$tf}{$individual_bases[$i]}[$i] > $score) {
			# body...
			$score = $motif_hash{$tf}{$individual_bases[$i]}[$i];
			$position = $i;
		}
		#print "$individual_bases[$i] $motif_hash{$tf}{$individual_bases[$i]}[$i] $score\n";
	}
	return(($score, $position));
}

sub min_pos {
	my ($tf, $pos) = @_;
#	print $pos."\t".$tf."\n";
	my $score = 111110;
	my $base = "N";
	foreach my $i ("A","C","G","T"){
		if ($motif_hash{$tf}{$i}[$pos] < $score) {
			# body...
			$base = $i;
			$score = $motif_hash{$tf}{$i}[$pos];
		}
		#print $i."\n";
	}
	return(($score, $base));
}

sub max_pos {
	my ($tf, $pos) = @_;
#	print $pos."\t".$tf."\n";
	my $score = 0;
	my $base = "N";
	foreach my $i ("A","C","G","T"){
		if ($motif_hash{$tf}{$i}[$pos] > $score) {
			# body...
			$base = $i;
			$score = $motif_hash{$tf}{$i}[$pos];
		}
		#print $i."\n";
	}
	return(($score, $base));
}


sub consensus {
	my $tf = shift;
	my $c_motif = "";
	foreach my $pos (0..(scalar @{$motif_hash{$tf}{"A"}} - 1)){
		my ($score, $base) = max_pos($tf,$pos);
		if ($score > 0.5) {
			$c_motif = $c_motif.$base;
		}else{
			$c_motif = $c_motif."-";
		}
		#print $row.",".$motif_hash{$tf}{"c"}[$row]."\t";
	}
	return($c_motif);

}


sub rc {
  my $string = shift;
  $string =~ tr/acgtACGT/tgcaTGCA/;
  $string = reverse($string);
  return($string);
}
