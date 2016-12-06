#! /usr/bin/perl

## Load all modules, hashes and databases ##

use warnings;
use strict;
use DBI;
use Bio::SeqIO;
use Bio::Index::Genbank;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::SearchIO;

## Now have major program ##

if (@ARGV != 2){
	die("This program will run the relative rate test\nSyntax: $0 [reciprocal file] [ORF files,]\n");
}

my $seqhash;

foreach my $f (split /,/, $ARGV[1]){
	my $inseq = Bio::SeqIO->new(-file => $f, -format => 'fasta');
	while (my $seq = $inseq->next_seq){
		$seqhash->{$seq->display_id} = $seq->translate()->seq;
	}
}

my $batch = "paupbatch.txt";
my $error = "error.txt";
open(ERR, ">$error");
	
my $loc = "Alignment_files_prot/";
if (! -e $loc){
	system("mkdir $loc");
}
my $loc2 = "Hyphy_prot/";

my $commandfile = "hyphy.txt";

my $outfile = "rrtest_prot.txt";
open(FINALOUT, ">$outfile");

print FINALOUT "seq_number\tprob\tquery_distance\ttest_distance\tmam_distance\tln_diff\tdistinct_sites\ttotal_sites\tave_percent_identity\toverall_percent_identity\tmin_undist\tmax_undist\tave_undist\n";

open(RECIP, $ARGV[0]);

while (<RECIP>){
	chomp;
	
	my @data = split /\t/;
	my $a = $data[0];
	
	my $qid = $data[1];
	my $qseq = $seqhash->{$data[1]};
	my $mamid = $data[2];
	my $mamseq = $seqhash->{$data[2]};
	my $testid = $data[3];
	my $testseq = $seqhash->{$data[3]};
		
	
	## Replace stops with X
	
	$qseq =~ s/\*/X/g;
	$mamseq =~ s/\*/X/g;
	$testseq =~ s/\*/X/g;
	

	# Need to alter the id's to remove any '.' #
	
	$qid =~ s/\.//g;
	$mamid =~ s/\.//g;
	$testid =~ s/\.//g;
	
	## We now want to create the multiple alignment ##

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
	if (! -e $folder1){
		my $mk = "mkdir $folder1";
		system($mk);
	} 
	my $folder2 = $loc . "/" . $first . "/" . $second;
	if (! -e $folder2){
		my $mk = "mkdir $folder2";
		system($mk);
	}
	
	my $outfile = $folder2 . "/" . $a . ".fasta";
	
	open (OUT, ">$outfile");
	print OUT ">" . $qid . "\n";
	print OUT $qseq . "\n\n";
	print OUT ">" . $mamid . "\n";
	print OUT $mamseq . "\n\n";
	print OUT ">" . $testid . "\n";
	print OUT $testseq . "\n\n";
	close(OUT);
	
	my $alnfile = $folder2 . "/" . $a . "_aln.fasta";
	
	my $muscle = "muscle -in $outfile -out $alnfile";
	system($muscle);
	
	if (! -e $alnfile){
		print ERR $a . " failed to align\n";
		next;
	}
	
	## Want to go through alignment and remove any positions with Xs
	
	my %lochash;
	my $lochash;
	my $alnfile2 = $folder2 . "/" . $a . "_aln2.fas";
	open (OUT, ">$alnfile2");
	
	my $exit = 0;
	
	my $inialnlen;
	my $alnin1 = Bio::AlignIO->new(-file => $alnfile, -format => 'fasta');
	while (my $aln = $alnin1->next_aln){
		
		$inialnlen = $aln->length;
		foreach my $seq ($aln->each_seq){
			my $str = $seq->seq;
			my @chars = split //, $str;
			my $pos = 0;
			foreach my $c (@chars){
				$pos++;
				if ($c eq 'X'){
					$lochash->{$pos}++;
				}
			}
		}
		foreach my $seq ($aln->each_seq){
			my $str = $seq->seq;
			my @chars = split //, $str;
			my $pos = 0;
			my @new;
			foreach my $c (@chars){
				$pos++;
				if (!exists ($lochash->{$pos})){
					push(@new, $c);
				}
			}
			my $n = "@new";
			$n =~ s/\s//g;
			print OUT ">" . $seq->display_id . "\n";
			print OUT "$n\n";
		}
	}
	
	close(OUT);
	
	
	## Run Gblocks with new parameters
	
	my $gblocks = "Gblocks $alnfile2 -t=p -b4=20 -b2=3";
	system($gblocks);
	
	my $gb = $alnfile2 . "-gb";
	if (!-e $gb || -z $gb){
		print ERR "$a: Gblocks Alignment failed\n";
		next;
	}
	
	## We now need to find which sequences are where in the alignment ##
	
	my $nex = $folder2 . "/" . $a . ".nex";
	
	my $str = Bio::AlignIO->new(-file => $gb, -format => 'fasta');
	my $alnout = Bio::AlignIO->new(-file => ">$nex", -format => 'nexus'); 
	my $aln = $str->next_aln;
	my $len = $aln->length;
	if ($len == 0){
		print ERR "$a: Gblocks gave aln length of 0\n";
		next;
	}
	$alnout->write_aln($aln);
	my $apid = $aln->percentage_identity;
	my $opid = $aln->overall_percentage_identity;
	my $seqcount = 0;
	my %hash;
	my $hash;
	my $failcount = 0;
	foreach my $seq ($aln->each_seq){
		$seqcount++;
		if ($seq->display_id eq $qid){
			$hash->{$qid} = $seqcount;
			$hash->{$seqcount}->{'id'} = $qid;
		}elsif($seq->display_id eq $testid){
			$hash->{$testid} = $seqcount;
			$hash->{$seqcount}->{'id'} = $testid;
		}elsif($seq->display_id eq $mamid){
			$hash->{$mamid} = $seqcount;
			$hash->{$seqcount}->{'id'} = $mamid;
		}else{
			print "This sequences doesn't match!\n";
		}
		#$seqcount++;
	}
	
	if ($failcount > 0){
		next;
	}
	
	
	## We can now build the batch file and run the relative rate test ##
	
	
	
	my $folder3 = $loc2 . "/" . $first;
	if (! -e $folder3){
		my $mk = "mkdir $folder3";
		system($mk);
	} 
	my $folder4 = $loc2 . "/" . $first . "/" . $second;
	if (! -e $folder4){
		my $mk = "mkdir $folder4";
		system($mk);
	}
	
	my $hyphyfile = $folder4 . "/" . $num . ".txt";
	
	if (-e $hyphyfile){
		my $rm = "rm $hyphyfile";
		system($rm);
	}
	
	my $cmd = "echo \'13\n1\n1\n/DataStore/LabData/Teri/Take2/RefuteExtavour/$gb\n" . $hash->{$mamid} . "\nJONES\n3\n1\n4\n\' | HYPHYMP | less > $hyphyfile";
	system($cmd);

	## We will now process the hyphy results into a table ##

	
	open (IN, $hyphyfile);
	
	my $tot = 0;
	my $dis = 0;
	my $kh = 0;
	my $k1 = 0;
	my $k2 = 0;
	my $huloc = $hash->{$mamid} - 1;
	my $spe1loc = $hash->{$qid} - 1;
	my $spe2loc = $hash->{$testid} - 1;
	my $lndif = 0;
	my $p = 1;
	my $huid = 0;
	my $spe1id = 0;
	my $spe2id = 0;
	
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^Total/){
			#print $line . "\n";
			my @stuff = split /:/, $line;
			$tot = $stuff[1];
			$tot =~ s/\;//;
		}
		if ($line =~ /^Distinct/){
			my @stuff = split /:/, $line;
			$dis = $stuff[1];
			$dis =~ s/\;//;
		}
		if ($line =~ /^Tree three/){
			#print $line . "\n";
			$line =~ s/Tree threeTaxaTree=\(//;
			$line =~ s/\)\;//;
			#print $line . "\n";
			my @stuff = split /,/, $line;
			$kh = $stuff[$huloc];
			$kh =~ s/^\w+://;
			$k1 = $stuff[$spe1loc];
			$k1 =~ s/^\w+://;
			$k2 = $stuff[$spe2loc];
			$k2 =~ s/^\w+://;
		}
		if ($line =~ /^-2\(Ln Like/){
			my @stuff = split /=/, $line;
			$lndif = $stuff[1];
		}
		if ($line =~ /^P-Val/){
			my @stuff = split /:/, $line;
			$p = $stuff[1];
		}
	}
	
	close(IN);
	
	my $undistfile = $folder2 . "/" . $a . "_undist.txt";
	
	open (BATCH, ">$batch");
		
	select BATCH;
	
	print "#NEXUS\n";
	print "execute $nex;\n";
	print "savedist file=$undistfile format=onecolumn;\n";
	print "quit;\n";
	print "end;\n";
	
	select STDOUT;
	
	my $paup = "paup -n $batch";
	system($paup);
	
	
	my $unmindist = 99999999999;
	my $unmaxdist = 0;
	my $untotdist = 0;
	
	open(IN, $undistfile);
	
	while (my $line = <IN>){
		chomp $line;
		my @data = split /\t/, $line;
		my $id1 = $data[0];
		my $id2 = $data[1];
		my $dist = $data[2];
		if ($dist > $unmaxdist){
			$unmaxdist = $dist;
		}
		if ($dist < $unmindist){
			$unmindist = $dist;
		}
		$untotdist = $untotdist + $dist;
	}
	
	close(IN);
	
	my $unavedist = $untotdist/3;
	
	#print "Min: $unmindist\tMax: $unmaxdist\tAve: $unavedist\n";
	
	#print "$a,$p,$k1,$k2,$kh,$lndif,$apid,$dis,$tot\n";
	
	print FINALOUT "$a\t$p\t$k1\t$k2\t$kh\t$lndif\t$dis\t$tot\t$apid\t$opid\t$unmindist\t$unmaxdist\t$unavedist\n";

}

close(FINALOUT);
close(RECIP);

exit;				