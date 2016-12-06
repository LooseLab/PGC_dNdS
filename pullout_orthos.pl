#! /usr/bin/perl

use warnings;
use strict;

if (@ARGV != 3){
	die("Syntax: $0 [file pair for 1st group,] [file pair for 2nd group,] [file pair for between groups,]\n");
}

my $hash1;
my @gr1;
my $gr1hash;
my $fcount = 0;

foreach my $f (split /,/, $ARGV[2]){
	$fcount++;
	print $f . "\n";
	open(IN, $f);
	while (<IN>){
		chomp;
		my @data = split /\t/;
		if ($fcount == 1){
			$hash1->{$data[1]} = $data[0];
		}elsif ($fcount == 2){
			if (exists $hash1->{$data[0]} && $hash1->{$data[0]} eq $data[1]){
				push(@gr1, $data[0] . ":" . $data[1]);
				$gr1hash->{$data[0] . ":" . $data[1]} = 1;
			}
		}
	}
	close(IN);
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

my $donehash;
my $out = "Pairs_1.txt";
my $gr = 0;
open(OUT, ">$out");

foreach my $id (keys %{$hash}){
	my $op = $hash->{$id};
	if (exists $donehash->{$id} || exists $donehash->{$op}){
		next;
	}
	if (exists $hash->{$op} && $hash->{$op} eq $id){
		$gr++;
		print OUT $gr . "\t" . $id . "\t" . $op . "\n";
		$donehash->{$id}++;
		$donehash->{$op}++;
	}
}
close(OUT);

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

$out = "Pairs_2.txt";
open(OUT, ">$out");

foreach my $id (keys %{$hash}){
	my $op = $hash->{$id};
	if (exists $donehash->{$id} || exists $donehash->{$op}){
		next;
	}
	if (exists $hash->{$op} && $hash->{$op} eq $id){
		$gr++;
		print OUT $gr . "\t" . $id . "\t" . $op . "\n";
		$donehash->{$id}++;
		$donehash->{$op}++;
	}
}
close(OUT);

$out = "4way.txt";
open(OUT, ">$out");

foreach my $p (keys %{$gr1hash}){
	my @pairs = split /:/, $p;
	my $n1;
	if (exists $hash->{$pairs[0]}){
		$n1 = $hash->{$pairs[0]};
		if (! exists $hash->{$n1} || $hash->{$n1} ne $pairs[0]){
			next;
		}
	}
	my $n2;
	if (exists $hash->{$pairs[1]}){
		$n2 = $hash->{$pairs[1]};
		if (! exists $hash->{$n2} || $hash->{$n2} ne $pairs[1]){
			next;
		}
	}
	if (!$n1 || !$n2){
		next;
	}
	$gr++;
	print OUT $gr . "\t" . $pairs[0] . "\t" . $n1 . "\t" . $pairs[1] . "\t" . $n2 . "\n";
}
close(OUT);

exit;
			
