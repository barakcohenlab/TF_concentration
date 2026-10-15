#!/usr/bin/perl
#
#
my $fasta = $ARGV[0];  ### actually fastq files as input
my $outname = $ARGV[1];

open(FASTA,"<$fasta");

my @current_seq_name;


my @mutation_array = ("TGA", "CGA", "TTA", "TGC", "GGA", "TGG", "AGA", "TGT", "TCA");
my %mut_count;


foreach my $i (@mutation_array){
        $mut_count{$i}=0;
}

my $hspFile="HSP";
my $minFile="minP";
my $pilotFile="pilot";


# open($hspFile, ">$outname"."_hsp_matched.txt");
# open($minFile, ">$outname"."_minP_matched.txt");
# open($pilotFile, ">$outname"."_pilot_matched.txt");

open(OUT1, ">$outname"."_all_matched.txt");
open(OUT2, ">$outname"."_unmatched.fastq");

my $numberOfLines = 0;
my %tmp_ID_loci_array;
my $curr_4lines = "";

my $formated_4lines = "";
my $unmatchMark = 0;
my $indexBarcode ;
my $prefixPart;
while (my $curr_seq = <FASTA>) {
        # body...
        chomp $curr_seq;
        $numberOfLines++;
        #if ($curr_seq =~ /^>(.*)$/) {
        if ($numberOfLines % 4 == 1) {
                # body...
                $curr_4lines = $curr_seq;
                my @id = split(/ /, $curr_seq);
                $formated_4lines = $id[0];
                #print "$id[0]\t";
        }elsif($numberOfLines % 4 == 2){
                $curr_4lines = $curr_4lines."\n".$curr_seq;
                if($curr_seq =~ /^(.*TAGACACCGCATGC)([A-Za-z]{10})GAGCTCTCGC/i){ ### hsp
                        $prefixPart = $1;
			$indexBarcode = $2;
                        $formated_4lines = "hsp\t$formated_4lines\t$2";
                }elsif($curr_seq =~ /^(.*TAGGTTGTTGCATGC)([A-Za-z]{10})GAGCTCTCGC/i){ ### minP
                        $prefixPart = $1;
                        $indexBarcode = $2;
                        $formated_4lines = "minP\t$formated_4lines\t$2";
                }else{
                        $unmatchMark = 1;
                }

        }elsif($numberOfLines % 4 == 3){
                $curr_4lines = $curr_4lines."\n".$curr_seq;

        }elsif($numberOfLines % 4 == 0){
                $curr_4lines = $curr_4lines."\n".$curr_seq;
                if ($unmatchMark > 0) {
                        # body...
                        print OUT2 "$curr_4lines\n";
                }else{
                        #my $tmpq1 = substr $curr_seq , 0 , length($prefixPart);
                        my $tmpq2 = substr $curr_seq , length($prefixPart), length($indexBarcode);
                        $formated_4lines = "$formated_4lines\t$tmpq2";
                        print OUT1 "$formated_4lines\n";
                }


                $unmatchMark = 0;
        }
        
}


close(FASTA);





