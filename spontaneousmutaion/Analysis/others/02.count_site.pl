#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
my %opt;
&getopts ('hg:f:o:', \%opt);

unless ($opt{f} && $opt{g}) {
    &HELP_MESSAGE();
    exit 1;
}
if ($opt{h}) {
    &HELP_MESSAGE();
    exit 0;
}


my $genome=$opt{f};
my $gff=$opt{g];
my $output=$opt{o} || "result.o";
&RUN ($genome, $gff, $output);

sub RUN {
my ($ref, $gff, $out) = @_;

#########read sequence#########
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
##############read input########
#open OUTPUT, ">$out" || die "cant save to $!\n";
open FI, "$in" || die "can't open file $!\n";
print OUTPUT "$in\t";
my $totalS=0;
my $totalN=0;
my $totalT=0;
#my %bpsregions;

while(<GFF>){
	chomp;
	my @item=split(/\t/);
	next unless ($item[2] eq "CDS");
	next if (/^#/);
	my $begin=$item[3];
	my $end=$item[4];
	my $strand=$item[6];
	my $genome_id=$item[0];
	my $current_seq=$seq{$genome_id};

	#test
	my $utr5;
	#my $bpsregion;
	#test

	if ($strand eq "+"){
	  foreach my $i ($begin .. $end){

		my $diff=$i-$begin;
		$code_pos=$diff % 3;
		$utr5=substr($current_seq,$begin-1,3);#test the begin
		
		die "someting wrong with $diff % 3 = $code_pos" if ($code_pos > 3 || $code_pos < 0);
		if ($code_pos == 0){
			$codon=substr($current_seq,$i-1,3);
		}
		elsif ($code_pos == 1){
			$codon=substr($current_seq,$i-2,3);
		}
		elsif ($code_pos == 2){
			$codon=substr($current_seq,$i-3,3);
		}
		else{
			die "someting wrong with $diff % 3 = $code_pos";
		}

		my ($s,$n)= &CODE_TEST($codon,$code_pos);
		$totalS+=$s;
		$totalN+=$n;
	  }
	}
	elsif ($item[13] eq "-"){
	  foreach my $i ($begin .. $end){	
		my $diff=$end-$i;
		$code_pos=$diff % 3;
		
		#test
		$utr5=substr($current_seq,$end-3,3);
		$utr5=reverse($utr5);
		$utr5=~tr/ACGT/TGCA/;
		#test end
		die "someting wrong with $diff % 3 = $code_pos" if ($code_pos > 3 || $code_pos < 0);
		
		if ($code_pos == 0){
			$codon=substr($current_seq,$change-3,3);
			$codon= reverse ($codon);#|| print "something wrong with $codon\n";
			$codon=~tr/ACGT/TGCA/;
		}
		elsif ($code_pos == 1){
			$codon=substr($current_seq,$change-2,3);
			$codon= reverse($codon);
                        $codon=~tr/ACGT/TGCA/;
		}
		elsif ($code_pos == 2){
			$codon=substr($current_seq,$change-1,3);
			$codon= reverse($codon);
                        $codon=~tr/ACGT/TGCA/;
		}else {
			die "someting wrong with $diff % 3 = $code_pos";
			#exit 1;
		}
		#$newbase=~tr/ACGT/TGCA/;
		my ($s,$n)= &CODE_TEST($codon,$code_pos);
                $totalS+=$s;
                $totalN+=$n;	
	  }
	}#end of -
	else{
		die "not understand of the string $item[13]\n";
		#exit 1;
	}#end of else

	
}
close FI;
my $total=$totalS+$totalN;
print OUTPUT "synonymous:$totalS\tnonsynonymous:$totalN\ttotal:$total\n";
#bpsregion;
#print OUTPUT "BPS_region:\n";
#foreach my $key (sort keys %bpsregions){
#	print OUTPUT "$key\t$bpsregions{$key}\n";
#}
#close OUTPUT;
}#end of sub

sub CODE_TEST {
	my ($code,$p)= @_;
	my $c_new=$code;
	my $ob=substr($code,$p,1);
	my @acgt=("A","C","G","T");
	my $Ns=0;
	my $Ss=0;
	foreach my $newbase (each @acgt){
		next if ($newbase eq $ob);
        	substr($c_new,$p,1)=$newbase;
        	my $new_amino=&CODE_CHANGE ($c_new);
        	my $old_amino=&CODE_CHANGE ($code);
        	#my $NorS;
        	if ($new_amino eq $old_amino){
                	#$NorS="S";
                	$Ss++;
                }
        	else{
                	#$NorS="N";
                	$Ns++;
        	}
	}
	return($Ss,$Ns);

}

sub CODE_CHANGE {
	my $code=$_[0];
	my $amino;
	#my %hash;
#	$hash{""}=
my %genetic_code=(
         'TTT'=>'F',
         'TTC'=>'F',
         'TTA'=>'L',
         'TTG'=>'L',
'TCT'=>'S',
'TCC'=>'S',
'TCA'=>'S',
'TCG'=>'S',
'TAT'=>'Y',
'TAC'=>'Y',
'TAA'=>'*',
'TAG'=>'*',
'TGT'=>'C',
'TGC'=>'C',
'TGA'=>'*',
'TGG'=>'W',
         'CTT'=>'L',
         'CTC'=>'L',
         'CTA'=>'L',
         'CTG'=>'L',
'CCT'=>'P',
'CCC'=>'P',
'CCA'=>'P',
'CCG'=>'P',
'CAT'=>'H',
'CAC'=>'H',
'CAA'=>'Q',
'CAG'=>'Q',
'CGC'=>'R',
'CGA'=>'R',
'CGG'=>'R',
'CGT'=>'R',
         'ATT'=>'I',
         'ATC'=>'I',
         'ATA'=>'I',
         'ATG'=>'M',
'ACT'=>'T',
'ACC'=>'T',
'ACA'=>'T',
'ACG'=>'T',
'AAT'=>'N',
'AAC'=>'N',
'AAA'=>'K',
'AAG'=>'K',
'AGT'=>'S',
'AGC'=>'S',
'AGA'=>'R',
'AGG'=>'R',
         'GTT'=>'V',
         'GTC'=>'V',
         'GTA'=>'V',
         'GTG'=>'V',
'GCT'=>'A',
'GCC'=>'A',
'GCA'=>'A',
'GCG'=>'A',
'GAT'=>'D',
'GAC'=>'D',
'GAA'=>'E',
'GAG'=>'E',
'GGT'=>'G',
'GGC'=>'G',
'GGA'=>'G',
'GGG'=>'G',
);
$amino=$genetic_code{$code};
#print "test:$amino\n";
return ($amino);

}
#####Usage#######


sub HELP_MESSAGE {
    print <<EOD;
Usage:  
       perl prog_name -f <genome> -g <gff> [-o <result>]
             -f genome       genome sequence
             -o result       result file name default: 'result.o'
	     -g gff          query input
             -h help         show this help
EOD
}

######Usage#####







