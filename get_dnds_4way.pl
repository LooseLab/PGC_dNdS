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
print SUM "Division\tSeqNumber\tI.ID1\tI.ID2\tP.ID1\tP.ID2\tNoSites\tI.dn\tI.ds\tI.dnds\tP.dn\tP.ds\tP.dnds\n";

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
	
	foreach my $d ($ds[0], $ds[5], $dn[0], $dn[5]){
		if ($d !~ /^\d/){
			$fail++;
		}
	}
	if ($fail > 0){
		print "$a contains a non-number\n";
		next;
	}
	foreach my $d ($ds[0], $ds[5]){
		if ($d == 0){
			$fail++;
		}
	}	
	
	if ($fail > 0){
		print "$a contains dS == 0\n";
		next;
	}
	
	my $dnds1 = "NA";
	if ($dn[1] =~ /^\d/ && $ds[1] =~ /^\d/){
		$dnds1 = ($dn[1]/$ds[1]);
	}
	my $dnds2 = "NA";
	if ($dn[2] =~ /^\d/ && $ds[2] =~ /^\d/){
		$dnds2 = ($dn[2]/$ds[2]);
	}
	my $dnds3 = "NA";
	if ($dn[3] =~ /^\d/ && $ds[3] =~ /^\d/){
		$dnds3 = ($dn[3]/$ds[3]);
	}
	my $dnds4 = "NA";
	if ($dn[4] =~ /^\d/ && $ds[4] =~ /^\d/){
		$dnds4 = ($dn[4]/$ds[4]);
	}
	
	open(OUT, ">$dndsfile");
	print OUT "\t" . join("\t", @ids) . "\n";
	print OUT $ids[0] . "\n";
	print OUT $ids[1] . "\t" . ($dn[0]/$ds[0]) . "\n";
	if (@dn > 1){
		print OUT $ids[2] . "\t" . $dnds1 . "\t" . $dnds2 . "\n";
		print OUT $ids[3] . "\t" . $dnds3 . "\t" . $dnds4 . "\t" . ($dn[5]/$ds[5]) . "\n";
	}
	close(OUT);
	
	my $ingroup1 = $dn[0]/$ds[0];
	my $ingroup2;
	if (@dn > 1){
		$ingroup2 = $dn[5]/$ds[5];
	}
	
	my $mode1 = "I";
	if ($ids[0] =~ /Xenopus/ || $ids[0] =~ /^\d*$/ || $ids[0] =~ /Xentr4/ || $ids[0] =~ /_FC$/ || $ids[0] =~ /_FP$/){
		$mode1 = "P";
	}
	
	if (@dn > 1){
		if ($mode1 eq "I"){
			print SUM $ARGV[0] . "\t" . $a . "\t" . $ids[0] . "\t" . $ids[1] . "\t" . $ids[2] . "\t" . $ids[3] . "\t" . $nsites . "\t" . $dn[0] . "\t" . $ds[0] . "\t" . $ingroup1 . "\t" . $dn[5] . "\t" . $ds[5] . "\t" . $ingroup2 . "\n";
		}else{
			print SUM $ARGV[0] . "\t" . $a . "\t" . $ids[2] . "\t" . $ids[3] . "\t" . $ids[0] . "\t" . $ids[1] . "\t" . $nsites . "\t" . $dn[5] . "\t" . $ds[5] . "\t" . $ingroup2 . "\t" . $dn[0] . "\t" . $ds[0] . "\t" . $ingroup1 . "\n";
		}
	}else{
		print SUM $ARGV[0] . "\t" . $a . "\t" . $ids[0] . "\t" . $ids[1] . "\t" . $nsites . "\t" . $dn[0] . "\t" . $ds[0] . "\t" . $ingroup1 . "\n";
	}

}

close(SUM);
exit;