#! /usr/bin/perl

use warnings;
use strict;

if (@ARGV != 4){
	die("Syntax: $0 [groupname] [start seqnum] [end seqnum] [folder]\n");
}

my $loc = $ARGV[3];
my $spe = $ARGV[0];
$spe =~ s/_/ /g;

my $out = "summary.txt";
open(SUM, ">$out");
print SUM "Division\tSeqNumber\tID1\tID2\tNoSites\tdn\tds\tdnds\n";

for (my $a = $ARGV[1]; $a <= $ARGV[2]; $a++){
	
	my $num = $a;
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
	my $folder2 = $loc . "/" . $first . "/" . $second;
	
	my $dnfile = $folder2 . "/" . $a . "_dn.meg";
	my $dsfile = $folder2 . "/" . $a . "_ds.meg";
	
	if (! -e $dnfile){
		next;
	}
	
	my $dndsfile = $folder2 . "/" . $a . "_dnds.txt";
	
	my @dn;
	my @ds;
	my @ids;
	my $nsites;
	
	my $fail = 0;
	
	open(IN, $dnfile);
	while (<IN>){
		chomp;
		if (/^\[[1-4]/){
			if (! /\#/){
				my @data = split /\s+/;
				foreach my $d (@data){
					if ($d =~ /\[/){
						next;
					}
#					if ($d !~ /^\d/){
#						$fail++;
#					}
					push(@dn, $d);
				}
			}else{
				my @data = split /\s+/;
				foreach my $d (@data){
					if ($d =~ /\[/){
						next;
					}
					$d =~ s/^#//;
					push (@ids, $d);
				}
			}
		}
		if (/^\s+No. of Sites/){
			my @data = split /:/;
			$nsites = $data[1];
			$nsites =~ s/^\s+//;
		}
	}
	close(IN);
	
	open(IN, $dsfile);
	while (<IN>){
		chomp;
		if (/^\[[2-4]/ && ! /\#/){
			my @data = split /\s+/;
			foreach my $d (@data){
				if ($d =~ /\[/){
					next;
				}
#				if ($d !~ /^\d/ || $d == 0){
#					$fail++;
#				}
				push(@ds, $d);
			}
		}
	}
	close(IN);
	
	foreach my $d ($ds[0], $dn[0]){
		if ($d !~ /^\d/){
			$fail++;
		}
	}
	if ($fail > 0){
		print "$a contains a non-number\n";
		next;
	}
	if ($ds[0] == 0){
		print "$a contains dS == 0\n";
		next;
	}

	
	open(OUT, ">$dndsfile");
	print OUT "\t" . join("\t", @ids) . "\n";
	print OUT $ids[0] . "\n";
	print OUT $ids[1] . "\t" . ($dn[0]/$ds[0]) . "\n";
	close(OUT);
	
	my $ingroup1 = $dn[0]/$ds[0];
	
	print SUM $ARGV[0] . "\t" . $a . "\t" . $ids[0] . "\t" . $ids[1] . "\t" . $nsites . "\t" . $dn[0] . "\t" . $ds[0] . "\t" . $ingroup1 . "\n";

}

close(SUM);
exit;