#!/usr/bin/perl
use strict;
#die unless scalar(@ARGV) == 2;
#my $gff = $ARGV[0];
my $inputDir = $ARGV[0];
#my $out_dir=$ARGV[1];
my $all_seq=4643117;
my $lengthcd=4066371;
my $lengthig=$all_seq-$lengthcd;
#my $lengthrp=31743;
#system ("mkdir -p $out_dir");

opendir(DIR, $inputDir) || die "Can't open input directory '$inputDir'\n";
my @files = readdir(DIR);
closedir(DIR);

print STDERR "file\tin run\tnot in run\t\+ 1\t\+ \>1\t\- 1\t\- \>1\t\+ 1AT\t\- 1AT\t\+ 1CG\t\- 1CG\n";

foreach my $vcf (@files) {
next unless $vcf =~ /\.vcf/;

my $file =$vcf;
$file=~s/\.vcf//;

my %hash;
#####
my $total=0;
my $in1=0;
my $del1=0;
my $dAT=0;
my $dCG=0;
my $inAT=0;
my $inCG=0;

my $ins=0;
my $dels=0;
my $run=0;
my $norun=0;

open(VCF,"./$inputDir/$vcf")||die "can't open $!\n";
while(<VCF>){
	chomp;
	my $currentline=$_;
	my @item=split(/\t/);
	$item[4] =~ s/\,\w+$//;
	my $change=length($item[4])-length($item[3]);
	next if ($change > 4 || $change < -4);
        if (length($item[3])<=2){
                $norun++;
        }else{
                $run++;
        }


	if ($change > 1 ){
		$ins++;
		print "$item[3]\t$item[4]\t$change\tin>1\n";
	}elsif($change == 1){
		$in1++;
		my $base=substr($item[4],(length($item[4])-1),1);
		if ($base =~ /[ATat]/){
			$inAT++;
		}elsif($base =~ /[CGcg]/){
			$inCG++;
		}else{
			print STDERR "Wrong: $item[3]\t$item[4]\t$base\n";
		}
		print "$item[3]\t$item[4]\t$change\tin=$base\n";
	}elsif($change == -1){
		$del1++;
		my $base=substr($item[3],(length($item[3])-1),1);
                if ($base =~ /[ATat]/){
                        $dAT++;
                }elsif($base =~ /[CGcg]/){
                        $dCG++;
                }else{
                        print STDERR "Wrong: $item[3]\t$item[4]\t$base\n";
                }
		print "$item[3]\t$item[4]\t$change\tdel=$base\n";
	}elsif($change < -1){
		$dels++;
		print "$item[3]\t$item[4]\t$change\tin<-1\n";
	}else{
		print STDERR "Wrong: $item[1]\t$item[3]\t$item[4]\t$change\n";
	}

}
print STDERR "$file\t$run\t$norun\t$in1\t$ins\t$del1\t$dels\t$inAT\t$dAT\t$inCG\t$dCG\n";

}  #end of foreach file;




