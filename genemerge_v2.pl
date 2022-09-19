#!/usr/bin/env perl

#[ HISTORY
#
#	v2.1 (4-Mar-2011)
#		Correct minor errors and output format.
#]

use strict;
use warnings;
use POSIX;
use Getopt::Std;

# -a association file
# -d description file
# -p population list
# -s study list
# -o output base name

our $VERSION = '2.1';

our($opt_a,$opt_d,$opt_p,$opt_s,$opt_o);
getopts('a:d:p:s:o:');

my $assoc_file = $opt_a; #80100011\sGO:0009987;GO:0048468;GO:0008219;GO:0012501\n
my $desc_file  = $opt_d; #GO:0000012\ssingle strand break repair\n
my $pop_file   = $opt_p; #80100006\n
my $study_file = $opt_s; #80100006\n
my $out_file   = $opt_o;

my %HoAssoc               = ();
my %HoPopAssocCount       = ();
#my %HoPopGenes           = ();
my %HoDesc                = ();
my %HoPopAssocFreq        = ();
#my %HoStudyGenes         = ();
my %HoStudyGeneAssocCount = ();
my %HoAssocStudyGene      = ();
my %HoStudyGeneAssocPVal  = ();
my $PopGeneNo             = '0';
my $StudyGeneNo           = '0';
my $StudyGeneNoAssoc      = '0';
my $StudyGeneUniqAssoc    = '0';
my $BonferroniCorr        = '0';

#gene->GO_term associations
open ASSOC, $assoc_file or die "$0 : can't open file $assoc_file : $!\n";
while (<ASSOC>) {
	chomp;
	my @assoc_line = split;
	my $assoc_gene = $assoc_line[0];
	my @assoc_go   = ();
	if ($assoc_line[1]) {
		@assoc_go = split /;/, $assoc_line[1];
		@{$HoAssoc{$assoc_gene}} = @assoc_go;
	}
}
close ASSOC;

#population, count GO_term, population number
open POP, $pop_file or die "$0 : can't open file $pop_file : $!\n";
while (<POP>) {
	chomp;
	my $PopGene = $_;
	$PopGeneNo++;
	#$HoPopGenes{$PopGene}++;
	if (exists $HoAssoc{$PopGene}) {
		foreach my $AssocGO (@{$HoAssoc{$PopGene}})  {
			$HoPopAssocCount{$AssocGO}++;
		}
	}
}
close POP;

#frequency of GO_term in population
foreach my $PopAssocCountKey (keys %HoPopAssocCount) {
	my $freq = $HoPopAssocCount{$PopAssocCountKey} / $PopGeneNo;
	$HoPopAssocFreq{$PopAssocCountKey} = $freq;
}

#description
open DESC, $desc_file or die "$0 : can't open file $desc_file : $!\n";
while (<DESC>) {
	chomp;
	my @desc_line = split /\s/, $_, 2;
	#my @desc_line = split /\t/, $_, 2;
	$HoDesc{$desc_line[0]} = $desc_line[1];
}
close DESC;

#study, count study_genes, count study_gene->GO_term associations, count study_genes with no GO_term
open STUDY, $study_file or die "$0 : can't open file $study_file  : $!\n";
while (<STUDY>) {
	chomp;
	my @study_line = split;
	my @StudyGenes = @study_line;
	foreach my $StudyGene (@StudyGenes) {
		$StudyGeneNo++;
		#$HoStudyGenes{$StudyGene}++;
		if (exists $HoAssoc{$StudyGene}) {
			foreach my $StudyGeneGO (@{$HoAssoc{$StudyGene}}) {
				$HoStudyGeneAssocCount{$StudyGeneGO}++;
				push @{$HoAssocStudyGene{$StudyGeneGO}}, $StudyGene;
			}
		} else {
			$StudyGeneNoAssoc++;
		}
	}
}
close STUDY;

#bonferroni correction
foreach my $StudyGeneAssocCountKey (keys %HoStudyGeneAssocCount) {
	$StudyGeneUniqAssoc++;
	if ($HoPopAssocFreq{$StudyGeneAssocCountKey} > (1 / $PopGeneNo)) {
		$BonferroniCorr++;
	}
}

#study_gene P-values based on hypergeometric distribution
foreach my $StudyGeneAssocCountKey2 (keys %HoStudyGeneAssocCount) {
	my $PVal  = '0';
	my $PValC = '0';
	my $N = $PopGeneNo;
	my $P = $HoPopAssocFreq{$StudyGeneAssocCountKey2};
	my $K = $StudyGeneNo;
	my $R = $HoStudyGeneAssocCount{$StudyGeneAssocCountKey2};
	if ($R != '1') {
		$PVal  = &hypergeometric($N,$P,$K,$R);
		$PValC = ($PVal * $BonferroniCorr >= 1) ? 1 : $PVal * $BonferroniCorr;
	} else {
		$PVal  = 'NA';
		$PValC = 'NA';
	}
	${$HoStudyGeneAssocPVal{$StudyGeneAssocCountKey2}}[0] = $PVal;
	${$HoStudyGeneAssocPVal{$StudyGeneAssocCountKey2}}[1] = $PValC;
}

#print all the shit
open OUT, ">$out_file.txt" or die "$0 : can't open file $out_file.txt : $!\n";
foreach my $result (sort keys %HoStudyGeneAssocCount) {
	my $GOterm       = $result;
	my $PopFreq      = $HoPopAssocFreq{$result};
	my $PopFrac      = $HoPopAssocCount{$result};
	my $PopFracAll   = $PopGeneNo;
	my $StudyFrac    = $HoStudyGeneAssocCount{$result};
	my $StudyFracAll = $StudyGeneNo;
	my $RawEs        = ${$HoStudyGeneAssocPVal{$result}}[0];
	my $EScore       = ${$HoStudyGeneAssocPVal{$result}}[1];
	my $Desc         = $HoDesc{$result};
	my @Contrib      = @{$HoAssocStudyGene{$result}};
	unless ($Desc) { $Desc = 'no description' }


	#printf OUT "%10s %-25s %5d %5d %5d %5d %-25s %-25s %-50s\n", $GOterm, $PopFreq, $PopFrac, $PopFracAll, $StudyFrac, $StudyFracAll, $RawEs, $EScore, $Desc;
	print OUT "$GOterm\t$PopFreq\t$PopFrac\t$PopFracAll\t$StudyFrac\t$StudyFracAll\t$RawEs\t$EScore\t$Desc\n";
	
	#pop_gene no
	#study_gene no
	#study_gene GO_term no (singletons)
	#study_genes with GO_term
	#study_genes with no GO_term
}
close OUT;
open OUT2, ">$out_file-contrib.txt" or die "$0 : can't open file $out_file-contrib.txt : $!\n";
foreach my $result2 (sort keys %HoStudyGeneAssocCount) {
	my $GOterm2     = $result2;
	my @Contrib     = @{$HoAssocStudyGene{$result2}};
	my $ContribNo   = scalar @Contrib;
	my $ContribJoin = join(';',@Contrib);
	my $Desc2       = $HoDesc{$GOterm2};
	print OUT2 "$GOterm2\t$Desc2\t$ContribNo\t$ContribJoin\n";
}
close OUT2;

sub hypergeometric {
	my $n = shift;
	my $p = shift;
	my $k = shift;
	my $r = shift;

	my $i    = '0';
	my $q    = '0';
	my $np   = '0';
	my $nq   = '0';
	my $top  = '0';
	my $sum  = '0';
	my $lfoo = '0';

	my $logNchooseK = '0';

	$q = 1 - $p;

	$np = floor( $n * $p + 0.5 );
	$nq = floor( $n * $q + 0.5 );

	$logNchooseK = &logNchooseK( $n, $k );

	$top = ($np < $k) ? $np : $k;

	$lfoo = &logNchooseK($np, $top) + &logNchooseK($n * (1 - $p), $k - $top);

	for ($i = $top ; $i >= $r ; $i--) {
		$sum += exp($lfoo - $logNchooseK);
		if ($i > $r) { $lfoo = $lfoo + log($i / ($np - $i + 1)) + log(($nq - $k + $i) / ($k - $i + 1)) }
	}
	return $sum;
}

sub logNchooseK {
	my $n = shift;
	my $k = shift;

	my $i = '0';
	my $result = '0';

	$k = ($k > ($n - $k)) ? $n - $k : $k;

	for ($i = $n ; $i > ($n - $k) ; $i--) { $result += log($i) }

	$result -= &lFactorial($k);

	return $result;
}

sub lFactorial {
	my $number = shift;
	my $result = '0';
	my $i      = '0';

	for ($i = 2 ; $i <= $number ; $i++) { $result += log($i) }

	return $result;
}
#v0.1 bioscripts