#!/usr/bin/perl
use strict;
use Getopt::Std;

my %opts = ('T' => '2285629', 'G' => '2357488', 'n'=>'10' , 'o' => 'mc_dir', 'a' => '30', 'b' => '30', 'c' => '10','d' => '10','e' => '10', 'f' => '10');
getopts('T:G:n:o:a:b:c:d:e:f:',\%opts);
die ("
Usage: mc.pl [options] <ref_AT> <ref_CG>

This script will generate a random substitutions along genome based on a fixed mutation number and ratio given by user.

Options: -T INT total length of AT in genome [$opts{T}]
         -C INT total length of CG in genome [$opts{G}]
	 -n INT simulation times [$opts{n}]
	 -o STR output dir [$opts{o}]
	 -a INT number of A:T to G:C substitution [$opts{a}] 
	 -b INT number of G:C to A:T substitution [$opts{b}] 
	 -c INT number of A:T to T:A substitution [$opts{c}] 
	 -d INT number of A:T to C:G substitution [$opts{d}] 
	 -e INT number of G:C to T:A substitution [$opts{e}] 
	 -f INT number of G:C to C:G substitution [$opts{f}] 
") unless scalar(@ARGV) == 2;

my $out=$opts{o};

my $ref_AT=$ARGV[0];#"notes_AT.name.sg.txt";
my $ref_CG=$ARGV[1];#"notes_CG.name.sg.txt";

my $at=$opts{T};
my $cg=$opts{G};
my $mcn=$opts{n};

my $at2gc=$opts{a};
my $gc2at=$opts{b};
my $at2ta=$opts{c};
my $at2cg=$opts{d};
my $gc2ta=$opts{e};
my $gc2cg=$opts{f};

my $m_at=$at2gc+$at2ta+$at2cg;
my $m_cg=$gc2at+$gc2ta+$gc2cg;

my $r_at2gc=$at2gc/$m_at;
my $r_at2ta=$at2ta/$m_at;
my $r_at2cg=$at2cg/$m_at;

my $r_gc2at=$gc2at/$m_cg;
my $r_gc2ta=$gc2ta/$m_cg;
my $r_gc2cg=$gc2cg/$m_cg;

my $total=$m_at+$m_cg;

my $gene=0;
my $inter=0;

`mkdir -p $out`;
foreach my $i (1 .. $mcn){

	open OUT,">./$out/mc_$i.txt" || die;
	my %pos_at;
	my %pos_cg;
	my $n_at=0;
	while($n_at < $m_at){
		my $r_at=int(rand($at))+1;
		next if (exists $pos_at{$r_at});
		my $rr=rand(1);
		if ($rr<$r_at2gc){
			$pos_at{$r_at}=1;
		}elsif($rr<($r_at2gc+$r_at2ta)){
			$pos_at{$r_at}=2;
		}else{
			$pos_at{$r_at}=3;
		}
		$n_at++;
	}
	my $n_cg=0;
	while($n_cg < $m_cg){
		my $r_cg=int(rand($cg))+1;
		next if (exists $pos_cg{$r_cg});
		my $rr=rand(1);
                if ($rr<$r_gc2at){
                        $pos_cg{$r_cg}=1;
                }elsif($rr<($r_gc2at+$r_gc2ta)){
                        $pos_cg{$r_cg}=2;
                }else{
                        $pos_cg{$r_cg}=3;
                }
		$n_cg++;
	}
	
	my $p=0;
	open AT,"$ref_AT" || die "cant open $!\n";
	while(<AT>){
		chomp;
		my $line=$_;
		$p++;  #my $p=(split(/\t/))[1];
		my @item=split(/\t/);
		my $base=$item[2];
		my $newbase;
		next unless(exists $pos_at{$p});
		if($base eq "A"){
			$newbase="G" if($pos_at{$p}==1);
                        $newbase="T" if($pos_at{$p}==2);
                        $newbase="C" if($pos_at{$p}==3);
		}else{
                        $newbase="C" if($pos_at{$p}==1);
                        $newbase="A" if($pos_at{$p}==2);
                        $newbase="G" if($pos_at{$p}==3);
		}

		if($item[5] eq "gene"){
			$gene++;
			}else{
			$inter++;
		}
		print OUT "$item[0]\t$item[1]\t$base\t$newbase\t$item[3]\t$item[4]\t$item[5]\t$item[6]\t$item[7]\n";
	}
	close AT;
	
	my $p=0;
	open CG,"$ref_CG" || die "cant open $!\n";
	while(<CG>){
		chomp;
		my $line=$_;
		$p++;
		my @item=split(/\t/);
		my $base=$item[2];
		my $newbase;
		next unless(exists $pos_cg{$p});
		#my $spec=rand(1)
                if($base eq "C"){
                        $newbase="T" if($pos_at{$p}==1);
                        $newbase="A" if($pos_at{$p}==2);
                        $newbase="G" if($pos_at{$p}==3);
                }else{
                        $newbase="A" if($pos_at{$p}==1);
                        $newbase="T" if($pos_at{$p}==2);
                        $newbase="C" if($pos_at{$p}==3);
                }
                if($item[5] eq "gene"){
                        $gene++;
                        }else{
                        $inter++;
                }
                print OUT "$item[0]\t$item[1]\t$base\t$newbase\t$item[3]\t$item[4]\t$item[5]\t$item[6]\t$item[7]\n";

	}
	close CG;

	close OUT;
}

my $mean_gene=$gene/1000;
my $mean_inter=$inter/1000;
print "gene:$mean_gene\tinter:$mean_inter\n";

