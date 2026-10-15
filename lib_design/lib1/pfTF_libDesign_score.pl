#!/usr/bin/perl
#
#

use List::Util qw[max];
use strict;
#use warnings;


my $inputMidCoord = $ARGV[0]; # geneMatch file, 6 cols: 0-chr, 1-mid, 2-mid+1, 3-geneID, 4-geneName, 5-strand
my $outputfile_prefix = $ARGV[1];
my $outputfile = $outputfile_prefix."_score.txt";

my ($tf_name, $up_or_down) = $inputMidCoord =~ /geneMatch_fimo_([A-Za-z0-9]*)_([A-Za-z]*).txt/;
print "$tf_name    $up_or_down\n";

open(MID,"<$inputMidCoord") || die "no such file $inputMidCoord\n";

my %coord_gene;

my %RNA_file;
$RNA_file{"up"} = "/scratch/bclab/yuning/libDesign/pfPilot/rna_expr_data/noNA_up_id_gene_fFC_hFC_fhFC.txt";
$RNA_file{"down"} = "/scratch/bclab/yuning/libDesign/pfPilot/rna_expr_data/noNA_down_id_gene_fFC_hFC_fhFC.txt";

#my @mutation_array = ("CGA", "TTA", "TGC", "GGA", "TGG", "AGA", "TGT", "TCA");

my %posneg = ("up" => 1, "down" => -1);


##### score ENCODE screen begin

open(SCREEN,"</scratch/bclab/yuning/libDesign/pfPilot/encode_screen/short_hg19_lift_simplified_removeLowDNase_hepG2.bed") || die "no such file\n";

my %screen_scores = ("CTCF-only,CTCF-bound"=>5, "dELS"=>9,"dELS,CTCF-bound"=>8,"DNase-H3K4me3"=>8,"DNase-H3K4me3,CTCF-bound"=>7,"DNase-only"=>7,"pELS"=>10,"pELS,CTCF-bound"=>9,"PLS"=>10,"PLS,CTCF-bound"=>9);
my %screen_scores_raw = ("CTCF-only,CTCF-bound"=>1, "dELS"=>2,"dELS,CTCF-bound"=>3,"DNase-H3K4me3"=>4,"DNase-H3K4me3,CTCF-bound"=>5,"DNase-only"=>6,"pELS"=>7,"pELS,CTCF-bound"=>8,"PLS"=>9,"PLS,CTCF-bound"=>10);

my %screen_coord_score_hash;
#my @tmp_midcoord_keys = keys %midcoord__scores;
while (my $curr_rna = <SCREEN>) {
        chomp $curr_rna;
        my @tmp = split /\t/, $curr_rna;
        $screen_coord_score_hash{$tmp[0]."_".$tmp[1]} = $tmp[3];
        #my $screen_anno= $tmp[3];
        #print scalar @tmp_midcoord_keys;
        #for (my $i = $tmp[1]; $i <= $tmp[2]; $i++) {
                #print "$i,";
                #my $i = $tmp[1];
                #my $tmp_pattern = $tmp[0]."_".$i."_";
                #foreach my $coord_gene (grep(/^$tmp_pattern/, @tmp_midcoord_keys)){
                #        print "\t".$coord_gene."\t".$tmp_pattern."\t\t";
                #        $midcoord__scores{$coord_gene}{"SCREEN_hepG2"} = $screen_scores{$screen_anno};

                #}
        #}
        #print "\n";
        #@tmp_midcoord_keys = grep !/^$tmp_pattern/, @tmp_midcoord_keys;
}
close(SCREEN);
##### score ENCODE screen end





##### score RNA FC begin

my $RNA_fc_penalty = -10;

my %RNA_gene_FC_F_H_FH;
open(RNA,"<$RNA_file{$up_or_down}") || die "no such file $RNA_file{$up_or_down}\n";
while (my $curr_rna = <RNA>) {
        chomp $curr_rna;
        my @tmp = split /\t/, $curr_rna;
        @{$RNA_gene_FC_F_H_FH{$tmp[0]}} = @tmp[2..4];
        #print "$tmp[0], @tmp[2..4]\n";
        #print "$RNA_gene_FC_F_H_FH{$tmp[0]}\n";
        #print "$RNA_gene_FC_F_H_FH{$tmp[0]}[0],$RNA_gene_FC_F_H_FH{$tmp[0]}[1],$RNA_gene_FC_F_H_FH{$tmp[0]}[2]\n";
}
close(RNA);

my %RNA_gene_TFscore_bothScore;
my %nameMatch = ("foxa1" => 0, "hnf4a" => 1);


if ($up_or_down eq "up") {
        my %tmp_above1;
        my %tmp_above1_both;
        #my @tmp_genelist;
        foreach my $gene (keys %RNA_gene_FC_F_H_FH){
                #print "$gene\t$RNA_gene_FC_F_H_FH{$gene}[$nameMatch{$tf_name}]\t$RNA_gene_FC_F_H_FH{$gene}[2]\n";
                if ($RNA_gene_FC_F_H_FH{$gene}[$nameMatch{$tf_name}] >= 1) {
                        $tmp_above1{$gene} = $RNA_gene_FC_F_H_FH{$gene}[$nameMatch{$tf_name}];
                }else{
                        $RNA_gene_TFscore_bothScore{$gene}{"1tf"} = $RNA_fc_penalty;
                }
                if ($RNA_gene_FC_F_H_FH{$gene}[2] >= 1) {
                        $tmp_above1_both{$gene} = $RNA_gene_FC_F_H_FH{$gene}[2];
                }else{
                        $RNA_gene_TFscore_bothScore{$gene}{"both"} = $RNA_fc_penalty;
                }
        }
        my $size_ofHash = keys %tmp_above1;
        my %rank_of;
        @rank_of{sort { $tmp_above1{$a} <=> $tmp_above1{$b} } keys %tmp_above1} = 1 .. $size_ofHash;
        foreach my $gene (keys %rank_of){
                $RNA_gene_TFscore_bothScore{$gene}{"1tf"}= int(1+ $rank_of{$gene} / ($size_ofHash+1) * 10);
                #print "1tf\t$gene\t$rank_of{$gene}\t$RNA_gene_FC_F_H_FH{$gene}[$nameMatch{$tf_name}]\t$RNA_gene_TFscore_bothScore{$gene}{'1tf'}\n";
        }
        my $size_ofHash = keys %tmp_above1_both;
        my %rank_of;
        @rank_of{sort { $tmp_above1_both{$a} <=> $tmp_above1_both{$b} } keys %tmp_above1_both} = 1 .. $size_ofHash;
        foreach my $gene (keys %rank_of){
                $RNA_gene_TFscore_bothScore{$gene}{"both"}= int(1+ $rank_of{$gene} / ($size_ofHash+1) * 10);
                #print "both\t$gene\t$rank_of{$gene}\t$RNA_gene_FC_F_H_FH{$gene}[2]\t$RNA_gene_TFscore_bothScore{$gene}{'both'}\n";
        }
}elsif($up_or_down eq "down"){
        my %tmp_above1;
        my %tmp_above1_both;
        #my @tmp_genelist;
        foreach my $gene (keys %RNA_gene_FC_F_H_FH){
                if ($RNA_gene_FC_F_H_FH{$gene}[$nameMatch{$tf_name}] <= -1) {
                        $tmp_above1{$gene} = $RNA_gene_FC_F_H_FH{$gene}[$nameMatch{$tf_name}];
                }else{
                        $RNA_gene_TFscore_bothScore{$gene}{"1tf"} = $RNA_fc_penalty;
                }
                if ($RNA_gene_FC_F_H_FH{$gene}[2] <= -1) {
                        $tmp_above1_both{$gene} = $RNA_gene_FC_F_H_FH{$gene}[2];
                }else{
                        $RNA_gene_TFscore_bothScore{$gene}{"both"} = $RNA_fc_penalty;
                }
        }
        my $size_ofHash = keys %tmp_above1;
        my %rank_of;
        @rank_of{sort { $tmp_above1{$b} <=> $tmp_above1{$a} } keys %tmp_above1} = 1 .. $size_ofHash;
        foreach my $gene (keys %rank_of){
                $RNA_gene_TFscore_bothScore{$gene}{"1tf"}= int(1+ $rank_of{$gene} / ($size_ofHash+1) * 10);
                #print "1tf\t$gene\t$rank_of{$gene}\t$RNA_gene_FC_F_H_FH{$gene}[$nameMatch{$tf_name}]\t$RNA_gene_TFscore_bothScore{$gene}{'1tf'}\n";
        }
        my $size_ofHash = keys %tmp_above1_both;
        my %rank_of;
        @rank_of{sort { $tmp_above1_both{$b} <=> $tmp_above1_both{$a} } keys %tmp_above1_both} = 1 .. $size_ofHash;
        foreach my $gene (keys %rank_of){
                $RNA_gene_TFscore_bothScore{$gene}{"both"}= int(1+ $rank_of{$gene} / ($size_ofHash+1) * 10);
                #print "both\t$gene\t$rank_of{$gene}\t$RNA_gene_FC_F_H_FH{$gene}[2]\t$RNA_gene_TFscore_bothScore{$gene}{'both'}\n";
        }
}else{
        die "error: cannot tell up or down\n";
}

# foreach my $mid (keys %midcoord__scores){
#         my $gene = $midcoord__geneID_geneName_strand{$mid}[0];
#         #print "$gene\t$RNA_gene_TFscore_bothScore{$gene}{'1tf'} \t $RNA_gene_TFscore_bothScore{$gene}{'both'}\n";
#         if ($RNA_gene_TFscore_bothScore{$gene}{"both"} > $RNA_gene_TFscore_bothScore{$gene}{"1tf"}) {
#                 $midcoord__scores{$mid}{"rnaFC"} = $RNA_gene_TFscore_bothScore{$gene}{"both"};
#         }else{
#                 $midcoord__scores{$mid}{"rnaFC"} = $RNA_gene_TFscore_bothScore{$gene}{"1tf"};
#         }
# }
##### score RNA FC end



##### score hepG2 MPRA 
open(MPRA,"</scratch/bclab/yuning/libDesign/pfPilot/hepG2_mpra_shendure/cres/short_hg19_element_quant_rep123.bed") || die "no such file\n";

my %hepG2_mpra_raw;
#my @tmp_midcoord_keys = keys %midcoord__scores;
while (my $curr_rna = <MPRA>) {
        chomp $curr_rna;
        my @tmp = split /\t/, $curr_rna;
        #my @tmp_genelist;
        if (exists($hepG2_mpra_raw{$tmp[0]."_".$tmp[1]})) {
                if ($tmp[4] > $hepG2_mpra_raw{$tmp[0]."_".$tmp[1]}) {
                        $hepG2_mpra_raw{$tmp[0]."_".$tmp[1]} = $tmp[4];
                }
        }else{
                $hepG2_mpra_raw{$tmp[0]."_".$tmp[1]} = $tmp[4];
        }
}

my %hepG2_mpra_score;
my %tmp_above0;

foreach my $gene (keys %hepG2_mpra_raw){
        if ($hepG2_mpra_raw{$gene} >= 0) {
                $tmp_above0{$gene} = $hepG2_mpra_raw{$gene};
        }else{
                $hepG2_mpra_score{$gene} = 0;
        }
}

my $size_ofHash = keys %tmp_above0;
my %rank_of;
@rank_of{sort { $tmp_above0{$a} <=> $tmp_above0{$b} } keys %tmp_above0} = 1 .. $size_ofHash;
foreach my $gene (keys %rank_of){
        $hepG2_mpra_score{$gene}= int(1+ $rank_of{$gene} / ($size_ofHash+1) * 10);
        #print "$gene\t$rank_of{$gene}\t$hepG2_mpra_raw{$gene}\t$hepG2_mpra_score{$gene}\n";
}

close(MPRA);
##### score hepG2 MPRA end


##### score K562 MPRA 
open(MPRA2,"</scratch/bclab/yuning/libDesign/pfPilot/k562_mpra_shendure//hg19/short_hg19_element_quant_rep123.bed") || die "no such file\n"; 

my %k562_mpra_raw;
#my @tmp_midcoord_keys = keys %midcoord__scores;
while (my $curr_rna = <MPRA2>) {
        chomp $curr_rna;
        my @tmp = split /\t/, $curr_rna;
        #my @tmp_genelist;
        if (exists($k562_mpra_raw{$tmp[0]."_".$tmp[1]})) {
                if ($tmp[4] > $k562_mpra_raw{$tmp[0]."_".$tmp[1]}) {
                        $k562_mpra_raw{$tmp[0]."_".$tmp[1]} = $tmp[4];
                }
        }else{
                $k562_mpra_raw{$tmp[0]."_".$tmp[1]} = $tmp[4];
        }
}

my %k562_mpra_score;
my %tmp_above0;

foreach my $gene (keys %k562_mpra_raw){
        if ($k562_mpra_raw{$gene} >= 0) {
                $tmp_above0{$gene} = $k562_mpra_raw{$gene};
        }else{
                $k562_mpra_score{$gene} = 0;
        }
}

my $size_ofHash = keys %tmp_above0;
my %rank_of;
@rank_of{sort { $tmp_above0{$a} <=> $tmp_above0{$b} } keys %tmp_above0} = 1 .. $size_ofHash;
foreach my $gene (keys %rank_of){
        $k562_mpra_score{$gene}= int(1+ $rank_of{$gene} / ($size_ofHash+1) * 10);
        #print "$gene\t$rank_of{$gene}\t$hepG2_mpra_raw{$gene}\t$hepG2_mpra_score{$gene}\n";
}

close(MPRA2);
##### score K562 MPRA end



##### score K562 ATAC before dox 
open(ATAC,"</scratch/bclab/yuning/libDesign/pfPilot/atac_Data/short_merged_all_uninduced_idr_ATACseq_005idr.bed") || die "no such file\n"; 

my %k562_atac_uninduced_score;
#my @tmp_midcoord_keys = keys %midcoord__scores;
while (my $curr_rna = <ATAC>) {
        chomp $curr_rna;
        my @tmp = split /\t/, $curr_rna;
        #my @tmp_genelist;        
        $k562_atac_uninduced_score{$tmp[0]."_".$tmp[1]} = 5;
}

close(ATAC);

##### score K562 ATAC before dox end


##### score tissue specificity begin
open(TISSUE,"</scratch/bclab/yuning/libDesign/pfPilot/proteinAtlasData/format_tissueSpecificity_proteinatlas.tsv") || die "no such file\n"; 

my %tissue_specificity_score;
#my @tmp_midcoord_keys = keys %midcoord__scores;
while (my $curr_rna = <TISSUE>) {
        chomp $curr_rna;
        my @tmp = split /\t/, $curr_rna;
        #my @tmp_genelist;
        if($tmp[2] eq "kidney" || $tmp[2] eq "pancreas") {
                if (! exists($tissue_specificity_score{$tmp[1]})) {
                        $tissue_specificity_score{$tmp[1]} = 2;
                }
        }
        if($tmp[2] eq "liver" || $tmp[2] eq "intestine"){
                $tissue_specificity_score{$tmp[1]} = 15;
        }        
        
}

close(TISSUE);
##### score tissue specificity end




##### score absolute expr after dox begin
my %RNA_tpm_file;
$RNA_tpm_file{"up"}{"foxa1"} = "/scratch/bclab/yuning/libDesign/pfPilot/rna_expr_data/foxa1_induced_TPM_mean.txt";
$RNA_tpm_file{"up"}{"hnf4a"} = "/scratch/bclab/yuning/libDesign/pfPilot/rna_expr_data/hnf4a_induced_TPM_mean.txt";

$RNA_tpm_file{"down"}{"foxa1"} = "/scratch/bclab/yuning/libDesign/pfPilot/rna_expr_data/foxa1_uninduced_TPM_mean.txt";
$RNA_tpm_file{"down"}{"hnf4a"} = "/scratch/bclab/yuning/libDesign/pfPilot/rna_expr_data/hnf4a_uninduced_TPM_mean.txt";

my $RNA_tpm_penalty = 0;

my %RNA_gene_tpm_score;
my %RNA_gene_tpm_raw;

open(RNA,"<$RNA_tpm_file{$up_or_down}{$tf_name}") || die "no such file $RNA_tpm_file{$up_or_down}{$tf_name}\n";
while (my $curr_rna = <RNA>) {
        chomp $curr_rna;
        my @tmp = split /\t/, $curr_rna;
        $RNA_gene_tpm_raw{$tmp[0]} = $tmp[1];
        if ($tmp[1] < 0.1) {
                $RNA_gene_tpm_score{$tmp[0]} = -$RNA_tpm_penalty;
        }else{
                $RNA_gene_tpm_score{$tmp[0]} = $RNA_tpm_penalty;
        }
        
        #print "$tmp[0], @tmp[2..4]\n";
        #print "$RNA_gene_FC_F_H_FH{$tmp[0]}\n";
        #print "$RNA_gene_FC_F_H_FH{$tmp[0]}[0],$RNA_gene_FC_F_H_FH{$tmp[0]}[1],$RNA_gene_FC_F_H_FH{$tmp[0]}[2]\n";
}
close(RNA);


##### score absolute expr after dox end



##### distance
open(DIST,"</scratch/bclab/yuning/libDesign/pfPilot/workplace/short_tss_gencode.bed") || die "no such file\n";
my %tss_allGenes;
while (my $curr_rna = <DIST>) {
        chomp $curr_rna;
        my @tmp = split /\t/, $curr_rna;
        $tss_allGenes{$tmp[3]}{"chr"} = $tmp[0];
        $tss_allGenes{$tmp[3]}{"tss"} = $tmp[1];
        $tss_allGenes{$tmp[3]}{"strand"} = $tmp[2];
        #print "$tmp[0], @tmp[2..4]\n";
        #print "$RNA_gene_FC_F_H_FH{$tmp[0]}\n";
        #print "$RNA_gene_FC_F_H_FH{$tmp[0]}[0],$RNA_gene_FC_F_H_FH{$tmp[0]}[1],$RNA_gene_FC_F_H_FH{$tmp[0]}[2]\n";
}
close(DIST);

##### distance end








### read all candidate file

my %midcoord__geneID_geneName_strand;
my %midcoord__scores;

while (my $curr_seq = <MID>) {
        # body...
        chomp $curr_seq;
        my @entry = split /\t/, $curr_seq;
        my $comma_gene = $entry[3];
        my ($comma_gene_beforeDot) = $entry[3] =~ /^([^\.]*)\./;
        #print "$comma_gene_beforeDot\n";
        #$comma_gene =~ tr/./_/;
        my $tmp_key_chr_ss_gene=$entry[0]."_".$entry[1]."_".$comma_gene;
        @{$midcoord__geneID_geneName_strand{$tmp_key_chr_ss_gene}}=@entry[3..5];
        if (exists($RNA_gene_tpm_score{$comma_gene})) {
                $midcoord__scores{$tmp_key_chr_ss_gene}{"rnaRPKM"}=$RNA_gene_tpm_score{$comma_gene};
                $midcoord__scores{$tmp_key_chr_ss_gene}{"rnaRPKM_raw"}=$RNA_gene_tpm_raw{$comma_gene};
        }else{
                $midcoord__scores{$tmp_key_chr_ss_gene}{"rnaRPKM"}=0;
                $midcoord__scores{$tmp_key_chr_ss_gene}{"rnaRPKM_raw"}="NaN";
        }
        $midcoord__scores{$tmp_key_chr_ss_gene}{"rnaFC"} = 0;
        #print "$tmp_key_chr_ss_gene\t$RNA_gene_TFscore_bothScore{$comma_gene}{'both'} \t$RNA_gene_TFscore_bothScore{$comma_gene}{'1tf'}\n";
        if ($RNA_gene_TFscore_bothScore{$comma_gene}{"both"} > $RNA_gene_TFscore_bothScore{$comma_gene}{"1tf"}) {
                 $midcoord__scores{$tmp_key_chr_ss_gene}{"rnaFC"} = $RNA_gene_TFscore_bothScore{$comma_gene}{"both"};
                 $midcoord__scores{$tmp_key_chr_ss_gene}{"rnaFC_raw"} = $RNA_gene_FC_F_H_FH{$comma_gene}[2];
        }else{
                 $midcoord__scores{$tmp_key_chr_ss_gene}{"rnaFC"} = $RNA_gene_TFscore_bothScore{$comma_gene}{"1tf"};
                 $midcoord__scores{$tmp_key_chr_ss_gene}{"rnaFC_raw"} = $RNA_gene_FC_F_H_FH{$comma_gene}[$nameMatch{$tf_name}];
        }
        if (exists($tissue_specificity_score{$comma_gene_beforeDot})) {
                 $midcoord__scores{$tmp_key_chr_ss_gene}{"tissueSpec"}=$tissue_specificity_score{$comma_gene_beforeDot};
        }else{
                $midcoord__scores{$tmp_key_chr_ss_gene}{"tissueSpec"}=0;
        }
        
        if (exists($screen_coord_score_hash{$entry[0]."_".$entry[1]})) {
                $midcoord__scores{$tmp_key_chr_ss_gene}{"SCREEN_hepG2"}=$screen_scores{$screen_coord_score_hash{$entry[0]."_".$entry[1]}};
                $midcoord__scores{$tmp_key_chr_ss_gene}{"SCREEN_hepG2_raw"}=$screen_scores_raw{$screen_coord_score_hash{$entry[0]."_".$entry[1]}};
        }else{
                $midcoord__scores{$tmp_key_chr_ss_gene}{"SCREEN_hepG2"}=0;
                $midcoord__scores{$tmp_key_chr_ss_gene}{"SCREEN_hepG2_raw"}="NaN";
        }
        if (exists($hepG2_mpra_score{$entry[0]."_".$entry[1]})) {
                $midcoord__scores{$tmp_key_chr_ss_gene}{"MPRA_hepG2"}=$hepG2_mpra_score{$entry[0]."_".$entry[1]} * $posneg{$up_or_down};
        }else{
                $midcoord__scores{$tmp_key_chr_ss_gene}{"MPRA_hepG2"}=0;
        }
        if (exists($hepG2_mpra_raw{$entry[0]."_".$entry[1]})) {
                $midcoord__scores{$tmp_key_chr_ss_gene}{"MPRA_hepG2_raw"}=$hepG2_mpra_raw{$entry[0]."_".$entry[1]};
        }else{
                $midcoord__scores{$tmp_key_chr_ss_gene}{"MPRA_hepG2_raw"}="NaN";
        }
        if (exists($k562_mpra_score{$entry[0]."_".$entry[1]})) {
                $midcoord__scores{$tmp_key_chr_ss_gene}{"MPRA_k562"}= -$k562_mpra_score{$entry[0]."_".$entry[1]} * $posneg{$up_or_down};
        }else{
                $midcoord__scores{$tmp_key_chr_ss_gene}{"MPRA_k562"}=0;
        }
        if (exists($k562_mpra_raw{$entry[0]."_".$entry[1]})) {
                $midcoord__scores{$tmp_key_chr_ss_gene}{"MPRA_k562_raw"}=$k562_mpra_raw{$entry[0]."_".$entry[1]};
        }else{
                $midcoord__scores{$tmp_key_chr_ss_gene}{"MPRA_k562_raw"}="NaN";
        }
        if (exists($k562_atac_uninduced_score{$entry[0]."_".$entry[1]})) {
                $midcoord__scores{$tmp_key_chr_ss_gene}{"ATAC_uninduced"}=$k562_atac_uninduced_score{$entry[0]."_".$entry[1]};
        }else{
                $midcoord__scores{$tmp_key_chr_ss_gene}{"ATAC_uninduced"}=0;
        }
        if (exists($tss_allGenes{$comma_gene})) {
                if ($entry[0] eq $tss_allGenes{$comma_gene}{"chr"}) {
                        my $tmp_dist = $tss_allGenes{$comma_gene}{"tss"} - $entry[1];
                        $midcoord__scores{$tmp_key_chr_ss_gene}{"distance"}= $tss_allGenes{$comma_gene}{"strand"}.":".$tmp_dist;
                }else{
                        $midcoord__scores{$tmp_key_chr_ss_gene}{"distance"}= "weiredGene";
                }
                
        }else{
                die "$comma_gene not in tss list\n";
        }
        
 
}
close(MID);




#################################### final output
my @title = ("rnaFC", "rnaRPKM","tissueSpec","SCREEN_hepG2","MPRA_hepG2","MPRA_k562","ATAC_uninduced");
my @title_raw = ("rnaFC_raw", "rnaRPKM_raw","SCREEN_hepG2_raw","MPRA_hepG2_raw","MPRA_k562_raw","distance");

my $outputfile_2 = $outputfile_prefix."_raw.txt";
open(OUTPUT1,">$outputfile") || die "Couldn't open 'all data' output file for $outputfile\n";
open(OUTPUT2,">$outputfile_2") || die "Couldn't open 'all data' output file for $outputfile\n";

print OUTPUT1 "chr\tmidPoint\tgeneID";
print OUTPUT2 "chr\tmidPoint\tgeneID";

foreach my $i (@title) {
        print OUTPUT1 "\t$i";
}
foreach my $i (@title_raw) {
        print OUTPUT2 "\t$i";
}
print OUTPUT1 "\ttotalScore\n";
print OUTPUT2 "\n";

foreach my $mid (keys %midcoord__scores){
        my ($chr,$ss,$gene) = split("_",$mid);
        print OUTPUT1 "$chr\t$ss\t$gene";
        print OUTPUT2 "$chr\t$ss\t$gene";
        my $totalScore=0;
        foreach my $i (@title) {
                print OUTPUT1 "\t$midcoord__scores{$mid}{$i}";
                $totalScore = $totalScore + $midcoord__scores{$mid}{$i};
        }
        foreach my $i (@title_raw) {
                print OUTPUT2 "\t$midcoord__scores{$mid}{$i}";
        }
        print OUTPUT1 "\t$totalScore\n";
        print OUTPUT2 "\n";
}



# foreach my $index (keys %index_hash){
#         if ($index_count{$index} > 1) {
#                 # body...
#                 #print "$index\t$index_count{$index}\t$index_hash{$index}\n";
#                 if(exists($mutation_count{$index_hash{$index}} )){
#                         $mutation_count{$index_hash{$index}} ++;
#                 }else{
#                         $mutation_count{$index_hash{$index}} = 1;
#                 }

#         }
# }

# foreach my $mut (@mutation_array){
#         print "$mut\t$mutation_count{$mut}\n";
# }






sub median{
        my $Aref = shift; #ref to input array;
        my $median;
        my $i; my $j; my $center;
        my @narm = grep(/[0-9]/, @{$Aref});
        my @tmpA; #will store the sorted $Aref;

        my $oddeven = (($#narm+1)%2);
        @tmpA = sort {$a <=> $b} @narm;

        if($oddeven){
                $center = $#tmpA/2;
                $median = $tmpA[$center];
        }
        else{
                $center = int($#tmpA/2);
                $median = ($tmpA[$center]+$tmpA[$center+1])/2;
        }
        return($median);
}


sub mean{
        my $target= shift @_;
        my $sum=0;
        my $count =0;
        foreach(@$target){
                $sum = $sum + $_;
                $count++;
        }
        my $aver = $sum / $count;
        return($aver);
}



