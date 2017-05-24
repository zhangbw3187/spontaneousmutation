#!/usr/bin/perl
use strict;
die ("
Usage: 00.note_all.pl <in.gff3> <vcf_dir> <output_dir>

This script will add the CDS/non-CDS information to each variants in the vcf files in <vcf_dir>, and write the output into <output_dir>.

") unless scalar(@ARGV) == 3;
my $gff = $ARGV[0];
my $inputDir = $ARGV[1];
my $out_dir=$ARGV[2];

system ("mkdir $out_dir");

opendir(DIR, $inputDir) || die "Can't open input directory '$inputDir'\n";
my @files = readdir(DIR);
closedir(DIR);

foreach my $vcf (@files) {
next unless $vcf =~ /\.vcf/;
my $file =$vcf;
$file=~ s/vcf/note/g; 

my %hash;
open(VCF,"./$inputDir/$vcf")||die "can't open $!\n";
while(<VCF>){
	chomp;
	next if (/^#/);
	my $line=$_;
	my $pos=(split(/\t/,$line))[1];
	$hash{$pos}=$line;
	}
close(VCF);
my %coding;
open(GFF,"$gff")||die "cant open $!\n";

open(OUT,">./$out_dir/$file")||die "can't create file ./$out_dir/$file\n";

while(<GFF>){
	chomp;
	my @item=split(/\t/);
	next unless ($item[2] eq "CDS");
	while(my ($key, $value) = each(%hash))
	{
		if($key>=$item[3] && $key <=$item[4])
		{
			print OUT "$value\t$item[3]\t$item[4]\t$item[2]\t$item[6]\t$item[8]\n";
			$coding{$key}=1;
		}
	}
}
close GFF;

while(my ($key,$val) = each (%hash)){
	unless(exists $coding{$key}){
		print OUT "$val\t.\t.\t.\t.\tinter-gene\n";
		}
	}

			
}




