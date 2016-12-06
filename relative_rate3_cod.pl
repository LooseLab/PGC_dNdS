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
use Cwd;


## Now have major program ##

if (@ARGV != 3){
	die("This program will run the relative rate test\nSyntax: $0 [reciprocal file] [type(orf/first/third)] [fasta files]\n");
}

my $seqhash;

foreach my $f (split /,/, $ARGV[2]){
	my $inseq = Bio::SeqIO->new(-file => $f, -format => 'fasta');
	while (my $seq = $inseq->next_seq){
		$seqhash->{$seq->display_id} = $seq->seq;
	}
}

my $t = $ARGV[1];

my $batch = "paupbatch.txt";
my $hyphyfile = "hyphytmp";
my $error = "error.txt";
open(ERR, ">$error");
	
my $loc = "Alignment_files_$t/";
if (! -e $loc){
	system("mkdir $loc");
}
my $loc2 = "Alignment_files_prot/";

my $pwd = cwd();
print $pwd . "\n";

my $table;
if ($t eq 'orf'){
	$table = "rrtest_orf.txt";
}elsif ($t eq 'first'){
	$table = "rrtest_first.txt";
}elsif ($t eq 'third'){	
	$table = "rrtest_third.txt";
}else{
	die("Please provide either 'orf', 'first' or 'third'\n");
}

my $commandfile = "hyphy.txt";

open(FINALOUT, ">$table");

print FINALOUT "seq_number\tprob\tquery_distance\ttest_distance\tmam_distance\tln_diff\tdistinct_sites\ttotal_sites\tave_percent_identity\toverall_percent_identity\tquery_gc\ttest_gc\tmam_gc\tave_gc\tmin_dist\tmax_dist\tave_dist\tmin_undist\tmax_undist\tave_undist\n";

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
		
	#print $testseq . "\n";
	
	# Remove '.' from ID and make sequences uppercase
		
	$qid =~ s/\.//g;
	$qseq =~ s/-//g;
	$qseq = uc($qseq);
	my $oriqid = $qid;
	$oriqid =~ s/_O$/_P/;
	
	$mamid =~ s/\.//g;
	$mamseq =~ s/-//g;
	$mamseq = uc($mamseq);
	#print $mamseq . "\n";
	
	$testid =~ s/\.//g;
	$testseq =~ s/-//g;
	$testseq = uc($testseq);
	my $oritid = $testid;
	$oritid =~ s/_O$/_P/;

	
	# We will now create the file in the correct subdirectory
	
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
	
	# Now need to create codon alignment - or find file if already exists 
	
	my $gb = "Alignment_files_orf/" . $first . "/" . $second . "/" . $a . "_alncod.fas";
	
	if (!-e $gb || -z $gb){
		
		my $progb = $loc2 . "/" . $first . "/" . $second . "/" . $a . "_aln2.fas-gb.htm";
		my $proaln = $loc2 . "/" . $first . "/" . $second . "/" . $a . "_aln.fasta";
	
		my %gbhash;
		my $gbhash;
		
		open (IN, $progb);
		while (my $line = <IN>){
			chomp $line;
			if ($line =~ /^Flanks: /){
				#print $line . "\n";
				$line =~ s/^Flanks: //;
				my @data = split /\]\s+\[/, $line;
				foreach my $pair (@data){
					$pair =~ s/\]//g;
					$pair =~ s/\[//g;
					my @sp = split /\s+/, $pair;
					my $st = $sp[0];
					my $en = $sp[1];
					for (my $x = $st;$x<=$en;$x++){
						$gbhash->{$x} = 1;
						#print "$x ";
					}
					#print "\n";
				}
			}
		}
		close(IN);	
		
		
		my %lochash;
		my $lochash;
		
		my @qaln;
		my @mamaln;
		my @testaln;
		
		my $alnin1 = Bio::AlignIO->new(-file => $proaln, -format => 'fasta');
		while (my $aln = $alnin1->next_aln){
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
				if ($seq->display_id eq $qid){
					#print "Have query seq\n";
					#print $seq->seq  . "\n";
					my $string = $seq->seq;
					my @data = split //, $string;
					my $count = 0;
					my $pos = 0;
					my $gb = 0;
					foreach my $aa (@data){
						$pos++;
						$gb++;
						#print $count . ": $aa\n";
						if (exists($lochash->{$pos})){
							$gb--;
							if ($aa ne '-'){
								$count += 3;
							}
							#print "Deleted this place!\n";
							next;
						}
						
						my $cod;
						if ($aa eq '-'){
							$cod = '---';
						}else{
							$cod = substr($qseq,$count,3);
							$count += 3;
						}
						if (exists($gbhash->{$gb})){
							push(@qaln, $cod);
						}
					}
				}
				if ($seq->display_id eq $mamid){
					#print "Have query seq\n";
					#print $seq->seq  . "\n";
					my $string = $seq->seq;
					my @data = split //, $string;
					my $count = 0;
					my $pos = 0;
					my $gb = 0;
					foreach my $aa (@data){
						$pos++;
						$gb++;
						#print $count . ": $aa\n";
						if (exists($lochash->{$pos})){
							$gb--;
							if ($aa ne '-'){
								$count += 3;
							}
							#print "Deleted this place!\n";
							next;
						}
						
						my $cod;
						if ($aa eq '-'){
							$cod = '---';
						}else{
							$cod = substr($mamseq,$count,3);
							$count += 3;
						}
						if (exists($gbhash->{$gb})){
							push(@mamaln, $cod);
						}
					}
				}
				if ($seq->display_id eq $testid){
					#print "Have query seq\n";
					#print $seq->seq  . "\n";
					my $string = $seq->seq;
					my @data = split //, $string;
					my $count = 0;
					my $pos = 0;
					my $gb = 0;
					foreach my $aa (@data){
						$pos++;
						$gb++;
						#print $count . ": $aa\n";
						if (exists($lochash->{$pos})){
							$gb--;
							if ($aa ne '-'){
								$count += 3;
							}
							#print "Deleted this place!\n";
							next;
						}
						
						my $cod;
						if ($aa eq '-'){
							$cod = '---';
						}else{
							$cod = substr($testseq,$count,3);
							$count += 3;
						}
						if (exists($gbhash->{$gb})){
							push(@testaln, $cod);
						}
					}
				}
			}
		}
		
		my $newqseq = "@qaln";
		$newqseq =~ s/\s//g;
		#print $newqseq . "\n";
		my $newmamseq = "@mamaln";
		$newmamseq =~ s/\s//g;
		#print $newmamseq . "\n";
		my $newtestseq = "@testaln";
		$newtestseq =~ s/\s//g;
		#print $newtestseq . "\n";
		
		if (!@mamaln || !$newqseq || length($newqseq) <= 1){
			print "No suitable ORF alignment\n";
			next;
		}
		
		open(OUT, ">$gb");
		print OUT ">$qid\n$newqseq\n";
		print OUT ">$testid\n$newtestseq\n";
		print OUT ">$mamid\n$newmamseq\n";
		close(OUT);
	}
	

	## Can now alter Gblocks alignment based on type of entry
	
	my $finalalnfile = $folder2 . "/" . $a . "_" . $t . ".fas";
	
	my $alnlen;
	
	if ($t ne 'orf'){
		open (OUT, ">$finalalnfile");
		my $seqin = Bio::SeqIO->new(-file => $gb, -format => 'fasta');
		while (my $seq = $seqin->next_seq){
			#print $seq->display_id . "\t" . $seq->seq . "\n";
			my $sequence = $seq->seq;
			my $total = length($sequence)/3;
			my $count = 0;
			my @data;
			for (my $c = 0;$c<length($sequence);$c+=3){
				if ($t eq 'first'){
					my $str = substr($sequence,$c, 2);
					push(@data, $str);
				}else{
					my $str = substr($sequence,$c+2, 1);
					push(@data, $str);
				}
			}
			print OUT ">" . $seq->display_id . "\n";
			my $tot = "@data";
			$tot =~ s/\s+//g;
			$alnlen = length($tot);
			print OUT "$tot\n";
		}
		close(OUT);
	}else{
		my $cp = "cp $gb $finalalnfile";
		system($cp);
	}
	
		
	## We now need to find which sequences are where in the alignment ##
	
	my $nex = $folder2 . "/" . $a . "_" . $t . ".nex";
	
	my $str = Bio::AlignIO->new(-file => $finalalnfile, -format => 'fasta');
	my $alnout = Bio::AlignIO->new(-file => ">$nex", -format => 'nexus'); 
	my $aln = $str->next_aln;
	$alnout->write_aln($aln);
	my $len = $aln->length;
	if ($len == 0){
		print ERR "$a: Gblocks gave aln length of 0\n";
		next;
	}
	my $apid = $aln->percentage_identity;
	my $opid = $aln->overall_percentage_identity;
	my $seqcount = 0;
	my %hash;
	my $hash;
	my $qgc;
	my $tgc;
	my $mgc;
	my $failcount = 0;
	foreach my $seq ($aln->each_seq){
		$seqcount++;
		my $string = $seq->seq;
		my @array = split //, $string;
		my $gcount = grep /G/, @array;
		my $ccount = grep /C/, @array;
		#print "G: $gcount\tC: $ccount\t";
		my $gc = (($gcount+$ccount)/length($string))*100;
		#print "%GC: $gc\n";
		if ($seq->display_id eq $qid){
			$qgc = $gc;
			$hash->{$qid} = $seqcount;
			$hash->{$seqcount}->{'id'} = $qid;
		}elsif($seq->display_id eq $testid){
			$tgc = $gc;
			$hash->{$testid} = $seqcount;
			$hash->{$seqcount}->{'id'} = $testid;
		}elsif($seq->display_id eq $mamid){
			$mgc = $gc;
			$hash->{$mamid} = $seqcount;
			$hash->{$seqcount}->{'id'} = $mamid;
		}else{
			$failcount++;
			print ERR "$a: Found unknown sequence in alignment\n";
		}
		#$seqcount++;
	}
	
	if ($failcount > 0){
		next;
	}
	my $avegc = ($qgc+$tgc+$mgc)/3;
	
	## We can now build the batch file and run the relative rate test ##
	
	
	my $cmd = "echo \'13\n1\n1\n$pwd/$finalalnfile\n" . $hash->{$mamid} . "\nCUSTOM\n012345\n\n3\n1\n4\n1\n1\n\' | HYPHYMP | less > $hyphyfile";
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
	
	my $undistfile = $folder2 . "/" . $a . "_" . $t . "_undist.txt";
	my $distfile = $folder2 . "/" . $a . "_" . $t . "_dist.txt";
	
	open (BATCH, ">$batch");
		
	select BATCH;
	
	print "#NEXUS\n";
	print "execute $nex;\n";
	print "savedist file=$undistfile format=onecolumn;\n";
	print "nj;\n";
	print "lscores 1/nst=6 rmatrix=estimate rates=gamma shape=estimate ncat=8;\n";
	print "lset rmatrix=previous shape=previous;\n";
	print "dset distance=ml;\n";
	print "savedist file=$distfile format=onecolumn;\n";
	print "quit;\n";
	print "end;\n";
	
	select STDOUT;
	
	my $paup = "paup -n $batch";
	system($paup);
	
	my $mindist = 99999999999;
	my $maxdist = 0;
	my $totdist = 0;
	
	open(IN, $distfile);
	
	while (my $line = <IN>){
		chomp $line;
		my @data = split /\t/, $line;
		my $id1 = $data[0];
		my $id2 = $data[1];
		my $dist = $data[2];
		if ($dist > $maxdist){
			$maxdist = $dist;
		}
		if ($dist < $mindist){
			$mindist = $dist;
		}
		$totdist = $totdist + $dist;
	}
	
	close(IN);
	
	my $avedist = $totdist/3;
	
	
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
	
	print FINALOUT "$a\t$p\t$k1\t$k2\t$kh\t$lndif\t$dis\t$tot\t$apid\t$opid\t$qgc\t$tgc\t$mgc\t$avegc\t$mindist\t$maxdist\t$avedist\t$unmindist\t$unmaxdist\t$unavedist\n";

}

close(FINALOUT);
close(RECIP);

exit;				