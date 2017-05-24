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
	if ($ob eq "G" || $ob eq "T"){
		$bpsregion=substr($current_seq,$change-3,4);
		$bpsregion=reverse($bpsregion);
		$bpsregion=~tr/ACGT/TGCA/;
		$newbase=~tr/ACGT/TGCA/;
		$obt=~tr/ACGT/TGCA/;
	}else{
		$bpsregion=substr($current_seq,$change-2,4);
	}
	#strand back
	if ($ob eq "G" || $ob eq "T"){
                $bpsregionb=substr($current_seq,$change-4,5);
                $bpsregionb=reverse($bpsregionb);
                $bpsregionb=~tr/ACGT/TGCA/;
                $newbaseb=~tr/ACGT/TGCA/;

        }else{
                $bpsregionb=substr($current_seq,$change-2,5);
        }
	my $tranV=1;
	$tranV=0 if ($obt eq "A" && $newbase eq "G");
	$tranV=0 if ($obt eq "C" && $newbase eq "T");

	next if ($tranV == 1);

	my $newbps=$bpsregion;
	substr($newbps,1,1)=$newbase;
	$bps{"$bpsregion"}=0 unless (exists $bps{"$bpsregion"});
	$bps{"$bpsregion"}++;

	my $newbpsb=$bpsregionb;
	substr($newbpsb,1,1)=$newbaseb;
	$bpsb{"$bpsregionb"}=0 unless (exists $bpsb{"$bpsregionb"});
	$bpsb{"$bpsregionb"}++;

}
print OUTPUT "\n5'NpMpNpN3'\tnumber\n";
foreach my $key (sort keys %bps){
	print OUTPUT "$key\t$bps{$key}\n";
}
print OUTPUT "\n5'NpMpNpNpN3'\tnumber\n";
foreach my $key (sort keys %bpsb){
        print OUTPUT "$key\t$bpsb{$key}\n";
}


} #end of sub

