#!/usr/bin/perl
use strict;
use Getopt::Std;

my %opts = ('l'=> 4643117, 'c'=>4066371);
getopts('l:c:',\%opts);
die ("
Usage: 02.count_SNP.pl [options] <noted_vcf_dir> <output_dir>

This script will count the SNP informations in the noted vcf files.

Options: -l INT length of the genome [$opts{l}]
	 -c INT length of CDS in the genome [$opts{c}] 

") unless scalar(@ARGV) == 2;
my $inputDir = $ARGV[0];
my $out_dir=$ARGV[1];
my $all_seq=$opts{l};
my $lengthcd=$opts{c};
my $lengthig=$all_seq-$lengthcd;
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
	#next if (/^#/);
	my @item=split(/\t/);
	my $change=$item[3]."--".$item[4];	
	next if ($currentline=~/^#/);
	if ($item[12]=~/\./){
		$interg++;
		}
	elsif($item[12]=~/CDS/){
		$coding++;
	}
	else{
		print STDERR "Sorry! I cannot recognize $item[12] as a feature\n";
	}
	$hash{$change}++;
	$total++;
}
close(VCF);


open(OUT,">./$out_dir/$file")||die "can't create file ./$out_dir/$file\n";

my $transition=$hash{"A--G"}+$hash{"T--C"}+$hash{"G--A"}+$hash{"C--T"};
my $at2gc=$hash{"A--G"}+$hash{"T--C"};
my $gc2at=$hash{"G--A"}+$hash{"C--T"};

my $transversion=$hash{"A--T"}+$hash{"T--A"}+$hash{"A--C"}+$hash{"T--G"}+$hash{"G--T"}+$hash{"C--A"}+$hash{"G--C"}+$hash{"C--G"};
my $at2ta=$hash{"A--T"}+$hash{"T--A"};
my $at2cg=$hash{"A--C"}+$hash{"T--G"};
my $gc2ta=$hash{"G--T"}+$hash{"C--A"};
my $gc2cg=$hash{"G--C"}+$hash{"C--G"};

print OUT "Type of substitution\,Number\,Fraction\n";
print OUT "Transitions\,$transition\,".$transition/$total."\n";
print OUT "A:T > G:C\,$at2gc\,".$at2gc/$total."\n";
print OUT "G:C > A:T\,$gc2at\,".$gc2at/$total."\n";
print OUT "Transversions\,$transversion\,".$transversion/$total."\n";
print OUT "A:T > T:A\,$at2ta\,".$at2ta/$total."\n";
print OUT "A:T > C:G\,$at2cg\,".$at2cg/$total."\n";
print OUT "G:C > T:A\,$gc2ta\,".$gc2ta/$total."\n";
print OUT "G:C > C:G\,$gc2cg\,".$gc2cg/$total."\n";
print OUT "Total\,$total\n";
print OUT "Position\,(coding:intergene\,$lengthcd\:$lengthig)\n";
print OUT "inter-gene\,$interg\,".$interg/$total."\,".$interg/$lengthig."\n";
print OUT "coding\,$coding\,".$coding/$total."\,".$coding/$lengthcd."\n";	
}




