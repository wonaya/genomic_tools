#/usr/bin/perl

#by Ethan Ford (with code lifed from the O'Reilly Perl Cookbook)
#Usage: sh randomLines.pl <fileName> <NumberOfLines>

use strict;
use warnings;

my $bedfilename = $ARGV[0];
my $numberoflines = $ARGV[1];

my $outfile = $ARGV[0];
$outfile =~ s/\.bed$/\.$numberoflines\.bed/;

open(TEXTFILE, $ARGV[0]);
open(TEXTOUT, ">$outfile");

my @textfile = <TEXTFILE>;

sub fisher_yates_shuffle {
	my $textfile = shift;
	my $i;
	for ($i = @$textfile; --$i; ) {
		my $j = int rand ($i+1);
		next if $i == $j;
		@$textfile[$i,$j] = @$textfile[$j,$i];
	}
}
fisher_yates_shuffle( \@textfile );


my @subsample = @textfile[0..$numberoflines-1];


foreach my $line (@subsample) {
  chomp($line);
  print TEXTOUT ("$line\n");
}


close(TEXTFILE); close(TEXTOUT);








