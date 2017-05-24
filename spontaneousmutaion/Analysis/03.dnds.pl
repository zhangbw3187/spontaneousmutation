#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
my %opt;
&getopts ('hd:n:f:o:b:', \%opt);

unless ($opt{f}) {
    &HELP_MESSAGE();
    exit 1;
}
if ($opt{h}) {
    &HELP_MESSAGE();
    exit 0;
}

my $genome=$opt{f};
my $output=$opt{o} || "result.o";


my %blosum;
my @blo_head;
if($opt{b}){
my $blo=$opt{b};
open BLO,$blo || die;
while(<BLO>){
        next if (/^#/);
        chomp;
        if(/^\s/){
                @blo_head=split(/\s+/);
                next;
        }
        my @it=split(/\s+/);
        foreach my $a (1 .. $#it){
                #next if (exists $blosum{"$head[$a]-$it[0]"});
                $blosum{"$it[0]-$blo_head[$a]"}=$it[$a];
        }
}
close BLO;
}

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
open FI, "$in" || die "can't open file $!\n";
print OUTPUT "$in\t";
my $totalS=0;
my $totalN=0;
my $totalT=0;
my $con_amino=0;
my $uncon_amino=0;
#my %bpsregions;

while(<FI>){
	chomp;
	my @item=split(/\t/);
	next unless ($item[12] eq "CDS");
	next if (/INDEL/);
	next if (/^#/);
	my $begin=$item[10];
	my $end=$item[11];
	my $change=$item[1];
	my $newbase=$item[4];
	my $genome_id=$item[0];
	my $codon;
	my $codon_new;
	my $current_seq=$seq{$genome_id};
	my $code_pos;

	#test
	my $utr5;
	#my $bpsregion;
	#test

	if ($item[13] eq "+"){
		my $diff=$change-$begin;
		$code_pos=$diff % 3;
		$utr5=substr($current_seq,$begin-1,3);#test the begin
		
		die "someting wrong with $diff % 3 = $code_pos" if ($code_pos > 3 || $code_pos < 0);
		if ($code_pos == 0){
			$codon=substr($current_seq,$change-1,3);
		}
		elsif ($code_pos == 1){
			$codon=substr($current_seq,$change-2,3);
		}
		elsif ($code_pos == 2){
			$codon=substr($current_seq,$change-3,3);
		}
		else{
			die "someting wrong with $diff % 3 = $code_pos";
		}

		##BPS region###
		#$bpsregion=substr($current_seq,$change-2,4);
		##BPS region

	}
	elsif ($item[13] eq "-"){
		my $diff=$end-$change;
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
		$newbase=~tr/ACGT/TGCA/;
		#BPS region#
		#$bpsregion = substr($current_seq,$change-3,4);
		#$bpsregion = reverse($bpsregion);
		#$bpsregion =~ tr/ACGT/TGCA/;

	}#end of -
	else{
		die "not understand of the string $item[13]\n";
		#exit 1;
	}#end of else
        $codon_new=$codon;
        substr($codon_new,$code_pos,1)=$newbase;
	my $new_amino=&CODE_CHANGE ($codon_new);
	my $old_amino=&CODE_CHANGE ($codon);
	my $NorS;

	if ($new_amino eq $old_amino){
		$NorS="S";
		$totalS++;
		}
	elsif($new_amino eq '*'){
		$NorS="N";
		$totalT++;
	}else{
		$NorS="N";
		$totalN++;
	}
	
	#BPS region##
#	$bpsregions{$bpsregion}=0 unless (exists $bpsregions{$bpsregion});
#	$bpsregions{$bpsregion}++;
	#print "$utr5\t$codon\t$codon_new\t$old_amino\t$new_amino\t$NorS\n";

	#compare to BLOSUM
	if($opt{b} && $NorS eq "N"){
	my $change="$old_amino-$new_amino";
	my $score=$blosum{$change};
	if($score >=0 ){
		$con_amino++;
		}else{
		$uncon_amino++;
		}
	}
	
}
close FI;
my $total=$totalS+$totalT+$totalN;
print OUTPUT "synonymous:$totalS\tterminal:$totalT\tnonsynonymous:$totalN\ttotal:$total\n";
#bpsregion;
#print OUTPUT "BPS_region:\n";
#foreach my $key (sort keys %bpsregions){
#	print OUTPUT "$key\t$bpsregions{$key}\n";
#}
#close OUTPUT;
if($opt{b}){
	print OUTPUT "con_amino:$con_amino\tuncon_amino:$uncon_amino\n";
}

}#end of sub

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
       perl prog_name -f <genome> -n <vcf> [-o <result>]
             -f genome       genome sequence
             -o result       result file name default: 'result.o'
	     -n input        query input
             -d dir          query input is a directory
             -h help         show this help
	     -b blosum       blosum score matrix that gives the conservative information 
EOD
}

######Usage#####







