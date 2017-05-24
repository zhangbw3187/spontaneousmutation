my $fasta=$ARGV[0];
my $pat=$ARGV[1];
open FA,"$fasta" ||die;
my $pat=uc($pat);
my $len=length($pat);
my $id;
my $seq;
while(<FA>){
	chomp;
	if(/^>/){
		$id=$_;
	}
	else{
		my $line=uc($_);
		$seq.=$line;
	}
}
my $t_pat = reverse ($pat);
$t_pat =~ tr/ACTG/TGAC/;

my $count=0;
my $f_len=length($seq);
foreach my $i (0 .. ($f_len - $len)){
	my $current=substr($seq,$i,$len);
	$count++ if ($current eq $pat || $current eq $t_pat);
	}
print "num of $pat/$t_pat : $count\n";
