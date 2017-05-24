#!/usr/bin/perl
my $anno=$ARGV[0];
my $vcf=$ARGV[1];

open AN,$anno || die;

my %id;
while(<AN>){
	chomp;
	next unless (/GN=/);
	my ($gid,$gnm) = $_ =~ /^(\w+)\t.*\bGN=(\w+)\s/;
	$id{$gid} = $gnm;
#	print "$gid\t$gnm\n";
}
close AN;

open VCF,$vcf || die;
while(<VCF>){
	chomp;
	my @it=split(/\t/);
	next unless ($it[12] eq "CDS");
	my ($gid) = $it[14] =~ /ID=(\w+);Name=/;
	$id{$gid}=$gid unless (exists $id{$gid});
	$it[14] = "ID=$gid;Name=$id{$gid};";
	my $line = join ("\t", @it);
	print "$line\n";
}
close VCF;


