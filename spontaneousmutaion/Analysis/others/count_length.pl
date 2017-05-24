#!/usr/bin/perl
my $total;
#my %cds;
while(<>){
	chomp;
	my @item=split(/\t/);
	next unless ($item[2] eq "CDS");
	my $length=$item[4]-$item[3];
	$total=$total+$length;
}
print $total;

