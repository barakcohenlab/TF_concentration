#!/usr/bin/perl
#
#
my $fasta = $ARGV[0];  ### actually bc count files as input
my $outname = $ARGV[1];

open(FASTA,"<$fasta");

my %bc_count;


while (my $curr_seq = <FASTA>) {
        # body...
        chomp $curr_seq;
        my @row = split("\t",$curr_seq);
        my $tmp_name = $row[0] =~ s/_rep[0-9]*$//r;
        #print $row[0]."\t".$tmp_name."\n";
        if (exists $bc_count{$tmp_name}) {
                # body...
                $bc_count{$tmp_name} = $bc_count{$tmp_name} + $row[2];
        }else{
                $bc_count{$tmp_name} = $row[2];
        }
        
}

close(FASTA);

open(OUT1, ">$outname"."_bcrCombined.txt");

foreach my $name (sort keys %bc_count){
        print OUT1 $name."\t".$bc_count{$name}."\n";
}


close(OUT1);



