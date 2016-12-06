#! /usr/bin/perl

use warnings;
use strict;
use Bio::SeqIO;

if (@ARGV != 1){
	die("Syntax: $0 [fasta file]\n");
}

my $inseq = Bio::SeqIO->new(-file => $ARGV[0], -format => 'fasta');

my $out = $ARGV[0] . "_clean";
open(OUT, ">$out");

while (my $seq = $inseq->next_seq){
	if ($seq->seq =~ /N/){ ## removes sequences with ambiguous nucleotides 
		next;
	}
	my $trans = $seq->translate();
	if ($trans->seq =~ /\*[A-Z]/){ ## removes sequences with an internal STOP codon
		next;
	}
	print OUT ">" . $seq->display_id . "\n";
	print OUT $seq->seq . "\n";
}
close(OUT);

exit;