#!/usr/bin/perl
use strict;
my %opts = ('l'=> 4643117, 'c'=>4066371);
getopts('l:c:',\%opts);
die ("
Usage: 02.count_SNP.pl [options] <noted_vcf_dir> <output_dir>

This script will count the SNPs in CDS as antisence strand.

Options: -l INT length of the genome [$opts{l}]
         -c INT length of CDS in the genome [$opts{c}] 

") unless scalar(@ARGV) == 2;

#my $gff = $ARGV[0];
my $inputDir = $ARGV[0];
my $out_dir=$ARGV[1];
my $all_seq=$opts{l};
my $lengthcd=$opts{c};
my $lengthig=$all_seq-$lengthcd;
#my $lengthrp=31743;
system ("mkdir -p $out_dir");

opendir(DIR, $inputDir) || die "Can't open input directory '$inputDir'\n";
my @files = readdir(DIR);
closedir(DIR);

foreach my $vcf (@files) {
next if $vcf =~ /^\./;

my $file =$vcf;
$file=~s/\.vcf//;
$file.=".csv";

my %hash;
$hash{"A--G"}=0;
$hash{"T--C"}=0;
#
$hash{"G--A"}=0;
$hash{"C--T"}=0;
######
$hash{"A--T"}=0;
$hash{"T--A"}=0;
#
$hash{"A--C"}=0;
$hash{"T--G"}=0;
#
$hash{"G--T"}=0;
$hash{"C--A"}=0;
#
$hash{"G--C"}=0;
$hash{"C--G"}=0;
#####
my $total=0;
my $interg=0;
my $coding=0;
my $repeat=0;
my $ncRNA=0;
open(VCF,"./$inputDir/$vcf")||die "can't open $!\n";
while(<VCF>){
	chomp;
	my $currentline=$_;
	my @item=split(/\t/);
	next if ($currentline=~/^#/);
	next unless ($item[12]=~/CDS/);
	$coding++;
	if ($item[13] eq "+"){
		$item[3] =~ tr/ACGT/TGCA/;
		$item[4] =~ tr/ACGT/TGCA/;
	}

	my $change = $item[3]."--".$item[4];

	$hash{$change}++;
	$total++;
}
close(VCF);

open(OUT,">./$out_dir/$file")||die "can't create file $file\n";

my $transition=$hash{"A--G"}+$hash{"T--C"}+$hash{"G--A"}+$hash{"C--T"};
my $at2gc=$hash{"A--G"}+$hash{"T--C"};
my $gc2at=$hash{"G--A"}+$hash{"C--T"};

print OUT "$file\:\n";
print OUT "G:C->A:T trasitions\t$gc2at\n";
print OUT "C\t$hash{'C--T'}\n";
print OUT "G\t$hash{'G--A'}\n";

print OUT "A:T->G:C trasitions\t$at2gc\n";
print OUT "A\t$hash{'A--G'}\n";
print OUT "T\t$hash{'T--C'}\n";


my $at2cg=$hash{"A--C"}+$hash{"T--G"};
my $gc2ta=$hash{"G--T"}+$hash{"C--A"};
my $gc2cg=$hash{"G--C"}+$hash{"C--G"};
my $at2ta=$hash{"A--T"}+$hash{"T--A"};

print OUT "G:C->T:A trasv\t$gc2ta\n";
print OUT "C\t$hash{'C--A'}\n";
print OUT "G\t$hash{'G--T'}\n";

print OUT "A:T->C:G trasv\t$at2cg\n";
print OUT "A\t$hash{'A--C'}\n";
print OUT "T\t$hash{'T--G'}\n";

print OUT "A:T->T:A trasv\t$at2ta\n";
print OUT "A\t$hash{'A--T'}\n";
print OUT "T\t$hash{'T--A'}\n";

print OUT "G:C->C:G trasv\t$gc2cg\n";
print OUT "G\t$hash{'G--C'}\n";
print OUT "C\t$hash{'C--G'}\n";

close OUT;	
}




