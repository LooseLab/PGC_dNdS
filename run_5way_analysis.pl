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

if (@ARGV != 2){
	die("This program will run the relative rate test\nSyntax: $0 [reciprocal file] [ORF files,]\n");
}

my $seqhash;

foreach my $f (split /,/, $ARGV[1]){
	my $inseq = Bio::SeqIO->new(-file => $f, -format => 'fasta');
	while (my $seq = $inseq->next_seq){
		my $id = $seq->display_id;
		$id =~ s/\./_/g;
		$seqhash->{$id} = $seq;
	}
}
	
my $loc = "5way_orf/";
if (! -e $loc){
	system("mkdir $loc");
}

my $batch = "paupbatch.txt";
my $commandfile = "hyphy.txt";
my $ctl = "code.ctl";

my $pwd = cwd();
print $pwd . "\n";

my $allout = "5way_results.txt";
open(ALLOUT, ">$allout");

print ALLOUT "SeqNumber\tid1\tid2\tid3\tid4\toutgroup\tAlnLen\tAve.pid\tOve.pid\tAve.gc\trrt_1v2.p\trrt_1v2.d1\trrt_1v2.d2\trrt_1v3.p\trrt_1v3.d1\trrt_1v3.d2\trrt_1v4.p\trrt_1v4.d1\trrt_1v4.d2\trrt_2v3.p\trrt_2v3.d1\trrt_2v3.d2\trrt_2v4.p\trrt_2v4.d1\trrt_2v4.d2\trrt_3v4.p\trrt_3v4.d1\trrt_3v4.d2\tb1.dn\tb1.ds\tb1.dnds\tb2.dn\tb2.ds\tb2.dnds\tb3.dn\tb3.ds\tb3.dnds\tb4.dn\tb4.ds\tb4.dnds\tbI.dn\tbI.ds\tbI.dnds\tbP.dn\tbP.ds\tbP.dnds\tbO.dn\tbO.ds\tbO.dnds\tp1v2.dn\tp1v2.ds\tp1v2.dnds\tp1v3.dn\tp1v3.ds\tp1v3.dnds\tp1v4.dn\tp1v4.ds\tp1v4.dnds\tp2v3.dn\tp2v3.ds\tp2v3.dnds\tp2v4.dn\tp2v4.ds\tp2v4.dnds\tp3v4.dn\tp3v4.ds\tp3v4.dnds\n";

open(RECIP, $ARGV[0]);

while (<RECIP>){
	chomp;
	
	$_ =~ s/\./_/g;
	
	my @data = split /\t/;
	my $a = $data[0];
	
	my $i1id = $data[1];
	my $i1orf = $seqhash->{$data[1]}->seq;
	my $i1prot = $seqhash->{$data[1]}->translate()->seq;
	my $i2id = $data[2];
	my $i2orf = $seqhash->{$data[2]}->seq;
	my $i2prot = $seqhash->{$data[2]}->translate()->seq;
	
	my $p1id = $data[3];
	my $p1orf = $seqhash->{$data[3]}->seq;
	my $p1prot = $seqhash->{$data[3]}->translate()->seq;
	my $p2id = $data[4];
	my $p2orf = $seqhash->{$data[4]}->seq;
	my $p2prot = $seqhash->{$data[4]}->translate()->seq;
	
	my $mamid = $data[5];
	my $mamorf = $seqhash->{$data[5]}->seq;
	my $mamprot = $seqhash->{$data[5]}->translate()->seq;
	
	my $idorder;
	
	$idorder->{$i1id} = "id1";
	$idorder->{$i2id} = "id2";	
	$idorder->{$p1id} = "id3";
	$idorder->{$p2id} = "id4";
	$idorder->{$mamid} = "outid";
	
	## Replace stops with X
	
	$i1prot =~ s/\*/X/g;
	$i2prot =~ s/\*/X/g;
	$p1prot =~ s/\*/X/g;
	$p2prot =~ s/\*/X/g;
	$mamprot =~ s/\*/X/g;
	

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
	
	### First want to make protein alignment
	
	my $protfile = $folder2 . "/" . $a . ".faa";
	
	open (OUT, ">$protfile");
	print OUT ">" . $i1id . "\n" . $i1prot . "\n";
	print OUT ">" . $i2id . "\n" . $i2prot . "\n";
	print OUT ">" . $p1id . "\n" . $p1prot . "\n";
	print OUT ">" . $p2id . "\n" . $p2prot . "\n";
	print OUT ">" . $mamid . "\n" . $mamprot . "\n";
	close(OUT);
	
	my $protalnfile = $folder2 . "/" . $a . "_aln.faa";
	
	my $muscle = "muscle -in $protfile -out $protalnfile";
	system($muscle);
	
	if (! -e $protalnfile){
		print $a . " failed to align\n";
		next;
	}
	
	## Remove any positions with X
	
	my $lochash;
	my $protalnfile2 = $folder2 . "/" . $a . "_aln2.faa";
	open (OUT, ">$protalnfile2");
	
	my $exit = 0;
	
	my $inialnlen;
	my $alnin1 = Bio::AlignIO->new(-file => $protalnfile, -format => 'fasta');
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
	
	## Run Gblocks
	
	my $gblocks = "Gblocks $protalnfile2 -t=p -b4=20 -b2=5";
	system($gblocks);
	
	my $gb = $protalnfile2 . "-gb";
	if (!-e $gb || -z $gb){
		print "$a: Gblocks Alignment failed\n";
		next;
	}
	
	### Now need to create codon alignment
	
	my $codalnfile = $folder2 . "/" . $a . "_alncod.fas";
	open(OUT, ">$codalnfile");
	
	my $gbhash;
	
	open (IN, $gb . ".htm");
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
	
	
	$lochash = ();
	my $fail = 0;
	
	
	my $alnin2 = Bio::AlignIO->new(-file => $protalnfile, -format => 'fasta');
	while (my $aln = $alnin2->next_aln){
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
			my $id = $seq->display_id;
			#print "Have query seq\n";
			#print $seq->seq  . "\n";
			my $string = $seq->seq;
			my @data = split //, $string;
			my $count = 0;
			my $pos = 0;
			my $gb = 0;
			my @qaln;
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
					$cod = substr($seqhash->{$id}->seq,$count,3);
					$count += 3;
				}
				if (exists($gbhash->{$gb})){
					push(@qaln, $cod);
				}
			}
			if (!@qaln || @qaln == 0){
				$fail++;
				last;
			}
				
			print OUT ">$id\n" . join("", @qaln) . "\n";
		}
	}
	close(OUT);
	
	if ($fail > 0){
		next;
	}
	
	## We now need to find which sequences are where in the alignment ##
	
	my $nex = $folder2 . "/" . $a . ".nex";
	
	my $str = Bio::AlignIO->new(-file => $codalnfile, -format => 'fasta');
	my $alnout = Bio::AlignIO->new(-file => ">$nex", -format => 'nexus'); 
	my $aln = $str->next_aln;
	$alnout->write_aln($aln);
	my $len = $aln->length;
	if ($len == 0){
		print "$a: Gblocks gave aln length of 0\n";
		next;
	}
	my $apid = $aln->percentage_identity;
	my $opid = $aln->overall_percentage_identity;
	my $seqcount = 0;
	my $hash;
	my $gchash;
	my $totgc = 0;
	my $failcount = 0;
	foreach my $seq ($aln->each_seq){
		$seqcount++;
		my $string = $seq->seq;
		my @array = split //, $string;
		my $gcount = grep /G/, @array;
		my $ccount = grep /C/, @array;
		#print "G: $gcount\tC: $ccount\t";
		my $gc = (($gcount+$ccount)/length($string))*100;
		$gchash->{$seq->display_id} = $gc;
		$totgc += $gc;
		#print "%GC: $gc\n";
		$hash->{$seq->display_id} = $seqcount;
		$hash->{$seqcount}->{'id'} = $seq->display_id;
		#$seqcount++;
	}
	
	if ($failcount > 0){
		next;
	}
	my $avegc = $totgc/$seqcount;
	
	print $avegc . "\n";
	
	my $treefile = $folder2 . "/" . $a . "_tree.txt";
	open(OUT, ">$treefile");
	print OUT "(($i1id,$i2id),($p1id,$p2id),$mamid);\n";
	close(OUT);
	
	my $rrtout = $folder2 . "/" . $a . "_rrt.txt";
	
	my $cmd = "echo \'13\n2\n1\n$pwd/$codalnfile\n" . $hash->{$mamid} . "\nCUSTOM\n012345\n3\n1\n4\n1\n1\n$pwd/$rrtout\n\' | HYPHYMP ";
	system($cmd);
	
	## Now run codeml to find branch dN/dS
	
	my $codout1 = $folder2 . "/" . $a . "_branchmodel.txt";
	
	open(OUT, ">$ctl");
	
	print OUT "seqfile = $codalnfile * sequence data filename\n";
	print OUT "treefile = $treefile      * tree structure file name\n";
	print OUT "outfile = $codout1           * main result file name\n";
	print OUT "noisy = 9  * 0,1,2,3,9: how much rubbish on the screen\n";
    print OUT "verbose = 1  * 0: concise; 1: detailed, 2: too much\n";
    print OUT "runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic\n";
    print OUT "seqtype = 1 * 1:codons; 2:AAs; 3:codons-->AAs\n";
    print OUT "CodonFreq = 1  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table\n";
    print OUT "clock = 0  * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis\n";
    print OUT "aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a\n";
    print OUT "aaRatefile = dat/jones.dat  * only used for aa seqs with model=empirical(_F)\n";
    print OUT "model = 1\n";
    print OUT "NSsites = 0  \n";
    print OUT "icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below\n";
    print OUT "Mgene = 0\n";
    print OUT "fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated\n";
    print OUT "kappa = 2  * initial or fixed kappa\n";
    print OUT "fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate \n";
    print OUT "omega = .4 * initial or fixed omega, for codons or codon-based AAs\n";
    print OUT "fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha\n";
    print OUT "alpha = 0. * initial or fixed alpha, 0:infinity (constant rate)\n";
    print OUT "Malpha = 0\n";
    print OUT "ncatG = 8  * # of categories in dG of NSsites models\n";
    print OUT "getSE = 0  * 0: don't want them, 1: want S.E.s of estimates\n";
    print OUT "RateAncestor = 1  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)\n";
    print OUT "Small_Diff = .5e-6\n";
    print OUT "cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?\n";
    print OUT "*  fix_blength = -1  * 0: ignore, -1: random, 1: initial, 2: fixed\n";
    print OUT "method = 0\n";
    
    close(OUT);
    
    system("codeml $ctl");
    
    ## Now run codeml to find pairwise dN/dS
	
	my $codout2 = $folder2 . "/" . $a . "_pairwise.txt";
	
	open(OUT, ">$ctl");
	
	print OUT "seqfile = $codalnfile * sequence data filename\n";
	print OUT "treefile = $treefile      * tree structure file name\n";
	print OUT "outfile = $codout2           * main result file name\n";
	print OUT "noisy = 9  * 0,1,2,3,9: how much rubbish on the screen\n";
    print OUT "verbose = 1  * 0: concise; 1: detailed, 2: too much\n";
    print OUT "runmode = -2  * 0: user tree;  1: semi-automatic;  2: automatic\n";
    print OUT "seqtype = 1 * 1:codons; 2:AAs; 3:codons-->AAs\n";
    print OUT "CodonFreq = 1  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table\n";
    print OUT "clock = 0  * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis\n";
    print OUT "aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a\n";
    print OUT "aaRatefile = dat/jones.dat  * only used for aa seqs with model=empirical(_F)\n";
    print OUT "model = 1\n";
    print OUT "NSsites = 0  \n";
    print OUT "icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below\n";
    print OUT "Mgene = 0\n";
    print OUT "fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated\n";
    print OUT "kappa = 2  * initial or fixed kappa\n";
    print OUT "fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate \n";
    print OUT "omega = .4 * initial or fixed omega, for codons or codon-based AAs\n";
    print OUT "fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha\n";
    print OUT "alpha = 0. * initial or fixed alpha, 0:infinity (constant rate)\n";
    print OUT "Malpha = 0\n";
    print OUT "ncatG = 8  * # of categories in dG of NSsites models\n";
    print OUT "getSE = 0  * 0: don't want them, 1: want S.E.s of estimates\n";
    print OUT "RateAncestor = 1  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)\n";
    print OUT "Small_Diff = .5e-6\n";
    print OUT "cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?\n";
    print OUT "*  fix_blength = -1  * 0: ignore, -1: random, 1: initial, 2: fixed\n";
    print OUT "method = 0\n";
    
    close(OUT);
    
    system("codeml $ctl");
	
	my $mlds = $folder2 . "/" . $a . "_ML.dS";
	my $mldn = $folder2 . "/" . $a . "_ML.dN";
	system("mv 2ML.dS $mlds; mv 2ML.dN $mldn");
	
	my $rrtres;
	my $pos1;
	my $pos2;
	
	open(IN, $rrtout);
	
	while (<IN>){
		chomp;
		if (/^,/){
			next;
		}
		if (/^Outgroup/){
			my @data = split /,/;
			for (my $a = 0; $a < @data; $a++){
				if ($data[$a] eq 'threeTaxaTree.FirstSpecies.t'){
					$pos1 = $a;
				}elsif ($data[$a] eq 'threeTaxaTree.SecondSpecies.t'){
					$pos2 = $a;
				}
			}
			next;
		}
		my @data = split /,/;
		my $comp;
		my $d1;
		my $d2;
		my $i1 = $idorder->{$data[1]};
		my $i2 = $idorder->{$data[2]};
		if ($i1 le $i2){
			$comp = $i1 . "@" . $i2;
			$d1 = $data[$pos1];
			$d2 = $data[$pos2];
		}else{
			$comp = $i2 . "@" . $i1;
			$d2 = $data[$pos1];
			$d1 = $data[$pos2];
		}
		
		$rrtres->{$comp}->{'prob'} = $data[24];
		$rrtres->{$comp}->{'d1'} = $d1;
		$rrtres->{$comp}->{'d2'} = $d2;
		
		print $comp . "\t" . $data[24] . "\t" . $data[$pos1] . "\t" . $data[$pos2] . "\n";
		
	}
	close(IN);
	
	my $branchdnds;
	
	open(IN, $codout1);
	my $gods = 0;
	my $godn = 0;
	
	while (<IN>){
		chomp;
		if (/^dS tree/){
			$gods++;
			next;
		}
		if ($gods > 0){
			print $_ . "\n";
			$_ =~ s/[(),:;]//g;
			my @data = split /\s+/;
			print "@data\n";
			$gods = 0;
			$branchdnds->{$data[0]}->{'ds'} = $data[1];
			$branchdnds->{$data[2]}->{'ds'} = $data[3];
			$branchdnds->{'Ibranch'}->{'ds'} = $data[4];
			$branchdnds->{$data[5]}->{'ds'} = $data[6];
			$branchdnds->{$data[7]}->{'ds'} = $data[8];
			$branchdnds->{'Pbranch'}->{'ds'} = $data[9];
			$branchdnds->{'Outgroup'}->{'ds'} = $data[11];
		}
		if (/^dN tree/){
			$godn++;
			next;
		}
		if ($godn > 0){
			print $_ . "\n";
			$_ =~ s/[(),:;]//g;
			my @data = split /\s+/;
			print "@data\n";
			$godn = 0;
			$branchdnds->{$data[0]}->{'dn'} = $data[1];
			$branchdnds->{$data[2]}->{'dn'} = $data[3];
			$branchdnds->{'Ibranch'}->{'dn'} = $data[4];
			$branchdnds->{$data[5]}->{'dn'} = $data[6];
			$branchdnds->{$data[7]}->{'dn'} = $data[8];
			$branchdnds->{'Pbranch'}->{'dn'} = $data[9];
			$branchdnds->{'Outgroup'}->{'dn'} = $data[11];
		}
	}
	close(IN);
	
	foreach my $branch (keys %{$branchdnds}){
		if ($branchdnds->{$branch}->{'ds'} > 0){
			$branchdnds->{$branch}->{'dnds'} = $branchdnds->{$branch}->{'dn'}/$branchdnds->{$branch}->{'ds'};
		}else{
			$branchdnds->{$branch}->{'dnds'} = -1;
		}
		print $branch . "\t" . $branchdnds->{$branch}->{'dn'} . "\t" . $branchdnds->{$branch}->{'ds'} . "\t" . $branchdnds->{$branch}->{'dnds'} . "\n";
	}
	
	my $pairwise;
	
	open(IN, $mlds);
	
	my @ids;
	
	while (<IN>){
		chomp;
		if (/^\s/){
			next;
		}
		my @data = split /\s+/;
		for (my $a = 0; $a < @data; $a++){
			if ($a == 0){
				push(@ids, $data[$a]);
			}else{
				my $comp;
				if ($data[0] eq $mamid || $ids[$a-1] eq $mamid){
					next;
				}
				if ($data[0] =~ /_A/ && $ids[$a-1] =~ /_A/){
					$comp = "betweenI";
				}elsif ($data[0] =~ /_F/ && $ids[$a-1] =~ /_F/){
					$comp = "betweenP";
				}else{
					my $i1 = $idorder->{$data[0]};
					my $i2 = $idorder->{$ids[$a-1]};
					if ($i1 le $i2){
						$comp = $i1 . "@" . $i2;
					}else{
						$comp = $i2 . "@" . $i1;
					}
				}
				$pairwise->{$comp}->{'ds'} = $data[$a];
			}
		}
	}
	close(IN);
	
	open(IN, $mldn);
	
	@ids = ();
	
	while (<IN>){
		chomp;
		if (/^\s/){
			next;
		}
		my @data = split /\s+/;
		for (my $a = 0; $a < @data; $a++){
			if ($a == 0){
				push(@ids, $data[$a]);
			}else{
				my $comp;
				if ($data[0] eq $mamid || $ids[$a-1] eq $mamid){
					next;
				}
				if ($data[0] =~ /_A/ && $ids[$a-1] =~ /_A/){
					$comp = "betweenI";
				}elsif ($data[0] =~ /_F/ && $ids[$a-1] =~ /_F/){
					$comp = "betweenP";
				}else{
					my $i1 = $idorder->{$data[0]};
					my $i2 = $idorder->{$ids[$a-1]};
					if ($i1 le $i2){
						$comp = $i1 . "@" . $i2;
					}else{
						$comp = $i2 . "@" . $i1;
					}
				}
				$pairwise->{$comp}->{'dn'} = $data[$a];
			}
		}
	}
	close(IN);
	
	foreach my $comp (keys %{$pairwise}){
		if ($pairwise->{$comp}->{'ds'} > 0){
			$pairwise->{$comp}->{'dnds'} = $pairwise->{$comp}->{'dn'}/$pairwise->{$comp}->{'ds'};
		}else{
			$pairwise->{$comp}->{'dnds'} = -1;
		}
		print $comp . "\t" . $pairwise->{$comp}->{'dn'} . "\t" . $pairwise->{$comp}->{'ds'} . "\t" . $pairwise->{$comp}->{'dnds'} . "\n";
	} 
	
	#"SeqNumber\tid1\tid2\tid3\tid4\toutgroup\tAlnLen\tAve.pid\tOve.pid\tAve.gc\trrt_1v2.p\trrt_1v2.d1\trrt_1v2.d2\trrt_1v3.p\trrt_1v3.d1\trrt_1v3.d2\trrt_1v4.p\trrt_1v4.d1\trrt_1v4.d2\trrt_2v3.p\trrt_2v3.d1\trrt_2v3.d2\trrt_2v4.p\trrt_2v4.d1\trrt_2v4.d2\trrt_3v4.p\trrt_3v4.d1\trrt_3v4.d2\tb1.dn\tb1.ds\tb1.dnds\tb2.dn\tb2.ds\tb2.dnds\tb3.dn\tb3.ds\tb3.dnds\tb4.dn\tb4.ds\tb4.dnds\tbI.dn\tbI.ds\tbI.dnds\tbP.dn\tbP.ds\tbP.dnds\tbO.dn\tbO.ds\tbO.dnds\tp1v2.dn\tp1v2.ds\tp1v2.dnds\tp1v3.dn\tp1v3.ds\tp1v3.dnds\tp1v4.dn\tp1v4.ds\tp1v4.dnds\tp2v3.dn\tp2v3.ds\tp2v3.dnds\tp2v4.dn\tp2v4.ds\tp2v4.dnds\tp3v4.dn\tp3v4.ds\tp3v4.dnds\n";
	
	print ALLOUT "$a\t$i1id\t$i2id\t$p1id\t$p2id\t$mamid\t$len\t$apid\t$opid\t$avegc\t" . $rrtres->{'id1@id2'}->{'prob'} . "\t" . $rrtres->{'id1@id2'}->{'d1'} . "\t" . $rrtres->{'id1@id2'}->{'d2'} . "\t" . $rrtres->{'id1@id3'}->{'prob'} . "\t" . $rrtres->{'id1@id3'}->{'d1'} . "\t" . $rrtres->{'id1@id3'}->{'d2'} . "\t" . $rrtres->{'id1@id4'}->{'prob'} . "\t" . $rrtres->{'id1@id4'}->{'d1'} . "\t" . $rrtres->{'id1@id4'}->{'d2'} . "\t" . $rrtres->{'id2@id3'}->{'prob'} . "\t" . $rrtres->{'id2@id3'}->{'d1'} . "\t" . $rrtres->{'id2@id3'}->{'d2'} . "\t" . $rrtres->{'id2@id4'}->{'prob'} . "\t" . $rrtres->{'id2@id4'}->{'d1'} . "\t" . $rrtres->{'id2@id4'}->{'d2'} . "\t" . $rrtres->{'id3@id4'}->{'prob'} . "\t" . $rrtres->{'id3@id4'}->{'d1'} . "\t" . $rrtres->{'id3@id4'}->{'d2'} . "\t" . $branchdnds->{$i1id}->{'dn'} . "\t" . $branchdnds->{$i1id}->{'ds'} . "\t" . $branchdnds->{$i1id}->{'dnds'} . "\t" . $branchdnds->{$i2id}->{'dn'} . "\t" . $branchdnds->{$i2id}->{'ds'} . "\t" . $branchdnds->{$i2id}->{'dnds'} . "\t" . $branchdnds->{$p1id}->{'dn'} . "\t" . $branchdnds->{$p1id}->{'ds'} . "\t" . $branchdnds->{$p1id}->{'dnds'} . "\t" . $branchdnds->{$p2id}->{'dn'} . "\t" . $branchdnds->{$p2id}->{'ds'} . "\t" . $branchdnds->{$p2id}->{'dnds'} . "\t" . $branchdnds->{'Ibranch'}->{'dn'} . "\t" . $branchdnds->{'Ibranch'}->{'ds'} . "\t" . $branchdnds->{'Ibranch'}->{'dnds'} . "\t" . $branchdnds->{'Pbranch'}->{'dn'} . "\t" . $branchdnds->{'Pbranch'}->{'ds'} . "\t" . $branchdnds->{'Pbranch'}->{'dnds'} . "\t" . $branchdnds->{'Outgroup'}->{'dn'} . "\t" . $branchdnds->{'Outgroup'}->{'ds'} . "\t" . $branchdnds->{'Outgroup'}->{'dnds'} . "\t" . $pairwise->{'betweenI'}->{'dn'} . "\t" . $pairwise->{'betweenI'}->{'ds'} . "\t" . $pairwise->{'betweenI'}->{'dnds'} . "\t" . $pairwise->{'id1@id3'}->{'dn'} . "\t" . $pairwise->{'id1@id3'}->{'ds'} . "\t" . $pairwise->{'id1@id3'}->{'dnds'} . "\t" . $pairwise->{'id1@id4'}->{'dn'} . "\t" . $pairwise->{'id1@id4'}->{'ds'} . "\t" . $pairwise->{'id1@id4'}->{'dnds'} . "\t" . $pairwise->{'id2@id3'}->{'dn'} . "\t" . $pairwise->{'id2@id3'}->{'ds'} . "\t" . $pairwise->{'id2@id3'}->{'dnds'} . "\t" . $pairwise->{'id2@id4'}->{'dn'} . "\t" . $pairwise->{'id2@id4'}->{'ds'} . "\t" . $pairwise->{'id2@id4'}->{'dnds'} . "\t" . $pairwise->{'betweenP'}->{'dn'} . "\t" . $pairwise->{'betweenP'}->{'ds'} . "\t" . $pairwise->{'betweenP'}->{'dnds'} . "\n";
	
	

}

close(ALLOUT);
close(RECIP);
exit;
