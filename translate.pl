#! /usr/bin/perl

use warnings;
use strict;
use Bio::SeqIO;

if (@ARGV == 0){
	die("Syntax: $0 [files]\n");
}

foreach my $in (@ARGV){
	my $out = $in . "_prot";
	
	my $inseq = Bio::SeqIO->new(-file => $in, -format => 'fasta');
	my $outseq = Bio::SeqIO->new(-file => ">$out", -format => 'fasta');
	
	while (my $seq = $inseq->next_seq){
		my $tmp = $seq->seq;
		$tmp =~ s/X/N/g;
		$seq->seq($tmp);
		my $trans = $seq->translate;
		$outseq->write_seq($trans);
	}
}

exit;