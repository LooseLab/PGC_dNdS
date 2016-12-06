#! /usr/bin/perl

use warnings;
use strict;

if (@ARGV != 3){
	die("Syntax: $0 [file pair for 1st group,] [file pair for 2nd group,] [start seqnumber]\n");
}


my $hash;

foreach my $f (split /,/, $ARGV[0]){
	print $f . "\n";
	open(IN, $f);
	while (<IN>){
		chomp;
		my @data = split /\t/;
		$hash->{$data[0]} = $data[1];
	}
	close(IN);
}

my $hash2;

foreach my $f (split /,/, $ARGV[1]){
	print $f . "\n";
	open(IN, $f);
	while (<IN>){
		chomp;
		my @data = split /\t/;
		$hash2->{$data[0]} = $data[1];
	}
	close(IN);
}

my $donehash;
my $out = "3way.txt";
my $gr = $ARGV[2];
open(OUT, ">$out");

foreach my $id (keys %{$hash}){
	my $op = $hash->{$id};
	if (exists $donehash->{$id} || exists $donehash->{$op}){
		next;
	}
	if (exists $hash->{$op} && $hash->{$op} eq $id){
		if (exists $hash2->{$id}){
			#print "First ID is in other pair!\n";
			my $op2 = $hash2->{$id};
			if (exists $hash2->{$op2} && $hash2->{$op2} eq $id){
				$gr++;
				print OUT $gr . "\t" . $id . "\t" . $op . "\t" . $op2 . "\n";
				$donehash->{$id}++;
				$donehash->{$op}++;
			}
		}elsif (exists $hash2->{$op}){
			#print "Second ID is in other pair!\n";
			my $op2 = $hash2->{$op};
			if (exists $hash2->{$op2} && $hash2->{$op2} eq $op){
				$gr++;
				print OUT $gr . "\t" . $op . "\t" . $id . "\t" . $op2 . "\n";
				$donehash->{$id}++;
				$donehash->{$op}++;
			}
		}
	}
}
close(OUT);

exit;