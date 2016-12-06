#! /usr/bin/perl

use warnings;
use strict;

if (@ARGV != 2){
	die("Syntax: $0 [4way groups] [extra recip,]\n");
}

my $hash;

foreach my $f (split /,/, $ARGV[1]){
	print $f . "\n";
	open(IN, $f);
	while (<IN>){
		chomp;
		my @data = split /\t/;
		$hash->{$data[0]} = $data[1];
	}
	close(IN);
}

my $out = "5way.txt";
open(OUT, ">$out");

open(IN, $ARGV[0]);

while (<IN>){
	chomp;
	my @data = split /\t/;
	my $new;
	my $g = shift @data;
	foreach my $i (@data){
		if (exists $hash->{$i}){
			my $op = $hash->{$i};
			if (exists $hash->{$op} && $hash->{$op} eq $i){
				$new = $op;
				last;
			}
		}
	}
	if (!$new){
		next;
	}
	print OUT $g . "\t" . join("\t", @data) . "\t" . $new . "\n";
}

close(IN);
close(OUT);