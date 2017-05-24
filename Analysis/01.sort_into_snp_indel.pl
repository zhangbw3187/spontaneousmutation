#!/usr/bin/perl
use strict;
die ("
Usage: 01.sort_into_snp_indel.pl <noted_vcf_dir> <snp_dir> <indel_dir>

I use this script to divide the variants in vcf files into SNP and InDel, I would analyze them separately, although this is not necessary. Good luck!

") unless scalar(@ARGV) == 3;

my $dir=$ARGV[0];
my $snp_dir=$ARGV[1];
my $indel_dir=$ARGV[2];
system ("mkdir $snp_dir") ;#|| die "cant create dir\n";
system ("mkdir $indel_dir") ;#||  die "cant create dir\n";

opendir(DIR, $dir) || die "Can't open input directory '$dir'\n";
my @files = readdir(DIR);
closedir(DIR);

foreach my $file (@files) {
next if $file =~ /^\./;

my $name=$file;
$name =~ s/vcf//g;
my $snpname=$name."snp.vcf";
my $indelname=$name."indel.vcf";

print "$snpname\t$indelname\n";

open IN,"./$dir/$file" || die;
open SNP,">./$snp_dir/$snpname"|| die;
open INDEL, ">./$indel_dir/$indelname"||die;
my %indels;
my %snps;
while (<IN>){
	chomp;
	next if (/^#/);
	my $posi=(split(/\t/))[1];
	if ($_=~/INDEL/){
		$indels{$posi}="$_";
	}
	else{
		$snps{$posi}="$_";
	}
}
foreach my $keys (sort {$a <=> $b} keys %indels){
	print INDEL "$indels{$keys}\n";
}

foreach my $keys (sort {$a <=> $b} keys %snps){
        print SNP "$snps{$keys}\n";
}

close IN;
close SNP;
close INDEL;

}

