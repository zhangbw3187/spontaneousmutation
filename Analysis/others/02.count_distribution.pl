#!/usr/bin/perl
use strict;
die unless scalar(@ARGV) == 3;
#my $gff = $ARGV[0];
my $inputDir = $ARGV[0];
my $out_dir=$ARGV[1];
my $interval=$ARGV[2];

system ("mkdir -p $out_dir");

opendir(DIR, $inputDir) || die "Can't open input directory '$inputDir'\n";
my @files = readdir(DIR);
closedir(DIR);

foreach my $vcf (@files) {
next if $vcf =~ /^\./;
my $file =$vcf."dis"; 
my $total_pos=0;
my $max=$total_pos+$interval;
my %hash;
$hash{"$total_pos\,$max"}=0;
my %poses;
open(VCF,"./$inputDir/$vcf")||die "can't open $!\n";
while(<VCF>){
	chomp;
	next if (/^#/);
	my $line=$_;
	my $pos=(split(/\t/,$line))[1];
	my $pos1=$pos+(4643117-4600935);
	if ($pos1 > 4643117){
		$pos1 = $pos1-4643117;
	}
	
	$poses{$pos1}=0 unless (exists $poses{$pos1});
	$poses{$pos1}++;
}
close VCF;

foreach my $pos (sort {$a<=>$b} keys %poses){
	my $num=$poses{$pos};
	my $current_pos=$pos-$total_pos;
	print "$current_pos\n";
	#the problem is that they are not in order!!!
	if ($current_pos <= $interval)
		{
		$hash{"$total_pos\,$max"}=$hash{"$total_pos\,$max"}+$num;
		}
	elsif($current_pos > $interval){
		$total_pos=$total_pos+$interval;
		$max=$max+$interval;
		$hash{"$total_pos\,$max"}=$num;
		}
	else{
		print "ERROR!$current_pos\n";
	}
	}
close(VCF);

open(OUT,">./$out_dir/$file")||die "can't create file ./$out_dir/$file\n";


#while(my ($key, $value) = each(%hash))
foreach my $key ( sort { $a <=> $b } keys %hash) {
	print OUT "$key\t$hash{$key}\n";
	}
			
}




