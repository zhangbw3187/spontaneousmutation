#!/usr/bin/perl
use strict;
die ("
Usage: 02.count_strand.pl <noted_vcf_dir> <st> <ed>

This script will reverse the SNPs and count them as opposite strand of ref-genome from <st> to <ed>
SNP in other regions will be regards as normal.

") unless scalar(@ARGV) == 3;

my $vcf = $ARGV[0];
my $st=$ARGV[1];
my $ed=$ARGV[2];
my $length=$ed-$st;

my $file =$vcf;
$file=~s/\.vcf//;
$file.=".lagging.txt";

my %base;
$base{A}=0;
$base{C}=0;
$base{G}=0;
$base{T}=0;
my %hash;
$hash{"A--G"}=0;
$hash{"T--C"}=0;
#
$hash{"G--A"}=0;
$hash{"C--T"}=0;
######
$hash{"A--C"}=0;
$hash{"T--G"}=0;
#
$hash{"G--T"}=0;
$hash{"C--A"}=0;
#
open(VCF,"$vcf")||die "can't open $!\n";
while(<VCF>){
	chomp;
	next if (/^#/);
	my @item=split(/\t/);
	my $pos= $item[1];
	#my $change = $item[3]."--".$item[4];	
	next if ($currentline=~/^#/);
	if ($pos >= $st && $pos < $ed){
		$item[3]=~ tr/ACGT/TGCA/;
		$item[4]=~ tr/ACGT/TGCA/;
	}
	my $change = $item[3]."--".$item[4];
	$hash{$change}++;
	$base{$item[3]}++;
}
close(VCF);

open(OUT,">$file")||die "can't create file $file\n";

my $transition=$hash{"A--G"}+$hash{"T--C"}+$hash{"G--A"}+$hash{"C--T"};
my $at2gc=$hash{"A--G"}+$hash{"T--C"};
my $gc2at=$hash{"G--A"}+$hash{"C--T"};

print OUT "$file\:\n";
print OUT "G:C->A:T trasitions\t$gc2at\n";
print OUT "C\t$hash{C--T}\n";
print OUT "G\t$hash{G--A}\n";

print OUT "A:T->G:C trasitions\t$at2gc\n";
print OUT "A\t$hash{A--G}\n";
print OUT "T\t$hash{T--C}\n";


my $at2cg=$hash{"A--C"}+$hash{"T--G"};
my $gc2ta=$hash{"G--T"}+$hash{"C--A"};

print OUT "G:C->T:A trasv\t$gc2ta\n";
print OUT "C\t$hash{C--A}\n";
print OUT "G\t$hash{G--T}\n";

print OUT "A:T->C:G trasv\t$at2cg\n";
print OUT "A\t$hash{A--C}\n";
print OUT "T\t$hash{T--G}\n";


print OUT "Total position\tA:T:C:G\t$base{A}\t$base{T}\t$base{C}\t$base{G}\n";
#print OUT "Position\,(coding:intergene:repeats\,$lengthcd\:$lengthig\:$lengthrp)\n";




