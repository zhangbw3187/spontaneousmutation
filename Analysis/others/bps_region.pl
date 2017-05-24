#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Std;
my %opt;
&getopts ('d:n:f:o:', \%opt);

my $genome=$opt{f};
my $output=$opt{o} || "result.o";
if ($opt{d}){
        my $inputdir=$opt{d};
        opendir DIR, $inputdir || die "cant open input directory $!\n";
        my @files=readdir(DIR);
        open OUTPUT, ">$output" || die "cant write into $!\n";
        foreach my $file (@files){
                next if $file =~ /^\./;
                my $input="./$inputdir/$file";
                &RUN ($input, $genome, $output);
        }
        close OUTPUT;
}
else{
exit 1 unless ($opt{n});
my $input=$opt{n};
&RUN ($input, $genome, $output);
}

sub RUN {
my ($in, $ref, $out) = @_;
open RE, "$ref" || die "can't open reference genome $!\n";
my $id;
my %seq;
while (<RE>){
        chomp;
        if (/\>/){
                s/^\>\s*//;
                s/\s+/ /g;
                $id=(split(/\s/))[0];
        }
        else{
                $seq{$id}.="$_";
        }
}
close RE;

open FI, "$in" || die "can't open file $!\n";
print OUTPUT "$in\t";
my %bps;
my %bpsb;

while(<FI>){
        chomp;
        my @item=split(/\t/);
        next if (/INDEL/);
        next if (/^#/);
        my $change=$item[1];
	my $ob=$item[3];
	my $obt=$ob;
        my $newbase=$item[4];
	my $newbaseb=$newbase;
        my $genome_id=$item[0];
	my $current_seq=$seq{$genome_id};
	my $bpsregion; ##=substr($current_seq,$change-2,2);
	my $bpsregionb;
	#strand in front of change;
	if ($ob eq "C" || $ob eq "T"){
		$bpsregion=substr($current_seq,$change-1,2);
		$bpsregion=reverse($bpsregion);
		$bpsregion=~tr/ACGT/TGCA/;
		$newbase=~tr/ACGT/TGCA/;
		$obt=~tr/ACGT/TGCA/;
	}else{
		$bpsregion=substr($current_seq,$change-2,2);
	}
	#strand back
	if ($ob eq "C" || $ob eq "T"){
                $bpsregionb=substr($current_seq,$change-2,2);
                $bpsregionb=reverse($bpsregionb);
                $bpsregionb=~tr/ACGT/TGCA/;
                $newbaseb=~tr/ACGT/TGCA/;

        }else{
                $bpsregionb=substr($current_seq,$change-1,2);
        }
	
        my $tranV=1;
        $tranV=0 if ($obt eq "A" && $newbase eq "G");
        $tranV=0 if ($obt eq "G" && $newbase eq "A");

        next if ($tranV == 0);


	my $newbps=$bpsregion;
	substr($newbps,0,1)=$newbase;
	$bps{$bpsregion}=0 unless (exists $bps{$bpsregion});
	$bps{$bpsregion}++;

	my $newbpsb=$bpsregionb;
	substr($newbpsb,0,1)=$newbaseb;
	$bpsb{$bpsregionb}=0 unless (exists $bpsb{$bpsregionb});
	$bpsb{$bpsregionb}++;

}
print OUTPUT "\nnormal_BPS\tnumber\tBPS_normal\tnumber\n";
#foreach my $key (sort keys %bps){
my $sumG=$bps{GG}+$bps{CG}+$bps{AG}+$bps{TG};
my $sumGb=$bpsb{GG}+$bpsb{GC}+$bpsb{GA}+$bpsb{GT};
print OUTPUT "GpG\t$bps{GG}\tGpG\t$bpsb{GG}\n";
print OUTPUT "CpG\t$bps{CG}\tGpC\t$bpsb{GC}\n";
print OUTPUT "ApG\t$bps{AG}\tGpA\t$bpsb{GA}\n";
print OUTPUT "TpG\t$bps{TG}\tGpT\t$bpsb{GT}\n";
print OUTPUT "sum\t$sumG\tsum\t$sumGb\n";

my $sumA=$bps{GA}+$bps{CA}+$bps{AA}+$bps{TA};
my $sumAb=$bpsb{AG}+$bpsb{AC}+$bpsb{AA}+$bpsb{AT};
print OUTPUT "GpA\t$bps{GA}\tApG\t$bpsb{AG}\n";
print OUTPUT "CpA\t$bps{CA}\tApC\t$bpsb{AC}\n";
print OUTPUT "ApA\t$bps{AA}\tApA\t$bpsb{AA}\n";
print OUTPUT "TpA\t$bps{TA}\tApT\t$bpsb{AT}\n";
print OUTPUT "sum\t$sumA\tsum\t$sumAb\n";

#}

} #end of sub

