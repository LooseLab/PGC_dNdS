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

my $commandfile = "hyphy.txt";

my $pwd = cwd();
print $pwd . "\n";

my $allout = "5way_results_RRT.txt";
open(ALLOUT, ">$allout");

print ALLOUT "SeqNumber\trrtf_1v2.p\trrtf_1v2.d1\trrtf_1v2.d2\trrtf_1v3.p\trrtf_1v3.d1\trrtf_1v3.d2\trrtf_1v4.p\trrtf_1v4.d1\trrtf_1v4.d2\trrtf_2v3.p\trrtf_2v3.d1\trrtf_2v3.d2\trrtf_2v4.p\trrtf_2v4.d1\trrtf_2v4.d2\trrtf_3v4.p\trrtf_3v4.d1\trrtf_3v4.d2\n";

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
	
	### First want to find orf alignment
	
	my $codalnfile = $folder2 . "/" . $a . "_alncod.fas";
	
	if (!-e $codalnfile){
		next;
	}
	open(IN, "$codalnfile");
	
	my $alnfile = $folder2 . "/" . $a . "_first.fas";
	
	open(OUT, ">$alnfile");
	
	
	my $failcount = 0;
	
	my $alnin2 = Bio::AlignIO->new(-file => $codalnfile, -format => 'fasta');
	while (my $aln = $alnin2->next_aln){
		if (length($aln) == 0){
			$failcount++;
			last;
		}
		foreach my $seq ($aln->each_seq){
			my $str = $seq->seq;
			my @new;
			for (my $c = 0;$c<length($str);$c+=3){
				my $sub = substr($str,$c, 2);
				push(@new, $sub);
			}
			print OUT ">" . $seq->display_id . "\n";
			print OUT join("", @new) . "\n";
		}
	}
	close(IN);
	close(OUT);
	
	if ($failcount > 0){
		next;
	}
	
	## We now need to find which sequences are where in the alignment ##
		
	my $str = Bio::AlignIO->new(-file => $alnfile, -format => 'fasta');
	my $aln = $str->next_aln;
	if (!$aln){
		next;
	}
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
	
	my $rrtout = $folder2 . "/" . $a . "_first_rrt.txt";
	
	my $cmd = "echo \'13\n2\n1\n$pwd/$alnfile\n" . $hash->{$mamid} . "\nCUSTOM\n012345\n3\n1\n4\n1\n1\n$pwd/$rrtout\n\' | HYPHYMP ";
	system($cmd);
	
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
	
	#"SeqNumber\tid1\tid2\tid3\tid4\toutgroup\tAlnLen\tAve.pid\tOve.pid\tAve.gc\trrt_1v2.p\trrt_1v2.d1\trrt_1v2.d2\trrt_1v3.p\trrt_1v3.d1\trrt_1v3.d2\trrt_1v4.p\trrt_1v4.d1\trrt_1v4.d2\trrt_2v3.p\trrt_2v3.d1\trrt_2v3.d2\trrt_2v4.p\trrt_2v4.d1\trrt_2v4.d2\trrt_3v4.p\trrt_3v4.d1\trrt_3v4.d2\tb1.dn\tb1.ds\tb1.dnds\tb2.dn\tb2.ds\tb2.dnds\tb3.dn\tb3.ds\tb3.dnds\tb4.dn\tb4.ds\tb4.dnds\tbI.dn\tbI.ds\tbI.dnds\tbP.dn\tbP.ds\tbP.dnds\tbO.dn\tbO.ds\tbO.dnds\tp1v2.dn\tp1v2.ds\tp1v2.dnds\tp1v3.dn\tp1v3.ds\tp1v3.dnds\tp1v4.dn\tp1v4.ds\tp1v4.dnds\tp2v3.dn\tp2v3.ds\tp2v3.dnds\tp2v4.dn\tp2v4.ds\tp2v4.dnds\tp3v4.dn\tp3v4.ds\tp3v4.dnds\n";
	
	print ALLOUT "$a\t" . $rrtres->{'id1@id2'}->{'prob'} . "\t" . $rrtres->{'id1@id2'}->{'d1'} . "\t" . $rrtres->{'id1@id2'}->{'d2'} . "\t" . $rrtres->{'id1@id3'}->{'prob'} . "\t" . $rrtres->{'id1@id3'}->{'d1'} . "\t" . $rrtres->{'id1@id3'}->{'d2'} . "\t" . $rrtres->{'id1@id4'}->{'prob'} . "\t" . $rrtres->{'id1@id4'}->{'d1'} . "\t" . $rrtres->{'id1@id4'}->{'d2'} . "\t" . $rrtres->{'id2@id3'}->{'prob'} . "\t" . $rrtres->{'id2@id3'}->{'d1'} . "\t" . $rrtres->{'id2@id3'}->{'d2'} . "\t" . $rrtres->{'id2@id4'}->{'prob'} . "\t" . $rrtres->{'id2@id4'}->{'d1'} . "\t" . $rrtres->{'id2@id4'}->{'d2'} . "\t" . $rrtres->{'id3@id4'}->{'prob'} . "\t" . $rrtres->{'id3@id4'}->{'d1'} . "\t" . $rrtres->{'id3@id4'}->{'d2'} . "\n";
	
	

}

close(ALLOUT);
close(RECIP);
exit;
