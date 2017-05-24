#!/usr/bin/perl
use strict;
die ("
Usage: note_list.pl <in.gff3> <input_file> <output_file>

This script will add the CDS/non-CDS information to each position in the mc ref-file.

") unless scalar(@ARGV) == 3;

my $gff = $ARGV[0];
my $inputfile = $ARGV[1];
my $out=$ARGV[2];

#my $file =$vcf;
#$file=~ s/vcf/note/g; 

my %hash;
open(IN,"./$inputfile")||die "can't open $!\n";
while(<IN>){
	chomp;
	next if (/^#/);
	my $line=$_;
	my $pos=(split(/\t/,$line))[1];
	$hash{$pos}=$line;
	}
close(IN);

my %coding;
open(GFF,"$gff")||die "cant open $!\n";

open(OUT,">$out")||die "can't create file $!\n";
my %genes;
while(<GFF>){
	chomp;
	my @item=split(/\t/);
	next unless ($item[2] eq "CDS");
	$genes{"$item[3],$item[4]"}="$item[2]\t$item[6]\t$item[8]";
	
}
close GFF;

while(my ($key, $value) = each(%hash)){
	foreach my $items (sort keys %genes){
	    my ($st,$ed)=split(/\,/,$items);
	    if($key>=$st && $key <=$ed){
		next if ($coding{$key} == 1);
		print OUT "$value\t$st\t$ed\t".$genes{"$st,$ed"}."\n";
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

			




