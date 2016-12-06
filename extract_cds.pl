#! /usr/bin/perl

use warnings;
use strict;
use Bio::SeqIO;

if (@ARGV != 1){
	die("Syntax: $0 [file]\n");
}

my $inseq = Bio::SeqIO->new(-file => $ARGV[0], -format => 'genbank');

my $out = $ARGV[0];
$out =~ s/\.gb$/.cds/;
open(OUT, ">$out");

while (my $seq = $inseq->next_seq){
	my @cds = grep { $_->primary_tag eq 'CDS' } $seq->get_SeqFeatures();
	foreach my $feature ( @cds ) {
		my $featureseq = $feature->spliced_seq;
		print OUT ">" . $featureseq->display_id . "\n" . $featureseq->seq . "\n";
	}
}
close(OUT);

exit;
	