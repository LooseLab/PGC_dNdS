#! /usr/bin/perl

use warnings;
use strict;
use Bio::SeqIO;
use DBI;
use Bio::AlignIO;

if (@ARGV != 3){
	die("Syntax: $0 [ortho file] [folder] [orf files,]\n");
}

my $loc = $ARGV[1];
if (! -e $loc){
	system("mkdir $loc");
}

my $seqhash;
foreach my $f (split /,/, $ARGV[2]){
	my $inseq = Bio::SeqIO->new(-file => $f, -format => 'fasta');
	while (my $seq = $inseq->next_seq){
		$seqhash->{$seq->display_id} = $seq->seq;
	}
}

open(IN, $ARGV[0]);

my $hash;

while (<IN>){
	chomp;
	my @data = split /\t/;
	my $gr = shift(@data);
	$hash->{$gr} = \@data;
}
close(IN);

foreach my $og (keys %{$hash}){
	
	my $num = $og;
	if (length($num) != 7){
		for(my $b = length($num);$b<7;$b++){
			my @line = split //, $num;
			unshift(@line,'0');
			$num = "@line";
			$num =~ s/\s//g;
		}
	}
	my $first = substr($num,0,2);
	my $second = substr($num,2,2);
	my $folder1 = $loc . "/" . $first;
	if (! -e $folder1){
		my $mk = "mkdir $folder1";
		system($mk);
	} 
	my $folder2 = $loc . "/" . $first . "/" . $second;
	if (! -e $folder2){
		my $mk = "mkdir $folder2";
		system($mk);
	}
	
	my $outfile = $folder2 . "/" . $og . ".fasta";
	
	open(OUT, ">$outfile");
	foreach my $id (@{$hash->{$og}}){
		print OUT ">" . $id . "\n" . $seqhash->{$id} . "\n";
	}
	close(OUT);
	
	my $alnfile = $folder2 . "/" . $og . "_aln";
	
	if (-e $alnfile . ".fasta"){
		system("rm " . $folder2 . "/" . $og . "_*");
	}
	
	my $dnfile = $folder2 . "/" . $og . "_dn";
	my $dsfile = $folder2 . "/" . $og . "_ds";
	
	system("megacc -a muscle_align_coding.mao -d $outfile -f Fasta -o $alnfile -n");
	
	system("megacc -a distance_estimation_pairwise_coding_jc_dn.mao -d $alnfile.fasta -o $dnfile -n");
	
	system("megacc -a distance_estimation_pairwise_coding_jc_ds.mao -d $alnfile.fasta -o $dsfile -n");
	
}
	
	