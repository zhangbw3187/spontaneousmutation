#!/usr/bin/perl
die ("
Usage: split.pl <reference.fasta>

This script will split reference fasta sequences into two ref-files with AT or CG sites seperately.

") unless scalar(@ARGV) == 1;

my $ref=$ARGV[0];  #input reference sequence
open REF,$ref || die "cannot open ref seq $!\n";
open OUTC,">ref_CG.txt" ||die;
open OUTA,">ref_AT.txt" ||die;
my $id;
my $pos=0;
while(<REF>){
	chomp;
	if(/>(.*)/){
		$id=$1;
		$pos=0;
	}
	else{
		my $seq=$_;
		$len=length($seq);
		foreach my $i (1 .. $len){
			$pos++;
			my $base=substr($seq,$i-1,1);
			if($base =~ /[ATat]/){
				print OUTA "$id\t$pos\t$base\n";
			}elsif($base =~ /[CGcg]/){
				print OUTC "$id\t$pos\t$base\n";
			}else{
				print STDERR "error: unknown base $base\n";
			}
			
		}
	}
}
close REF;
