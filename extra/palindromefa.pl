#! /usr/bin/perl -w
use strict;
die "Usage: $0 fafile\n" unless $ARGV[0];
my ($infile, $ct, $head, $seq, @lines) = ($ARGV[0], 0);
open IN, $infile;
unlink "$infile.pal";
open OUT, "> $infile.pal";
while (<IN>) {
 if (/^>(\S+)/) {
  if ($seq) {push @lines, ">$head\n$seq\n"; $ct ++; $seq = ''}
  $head = $1;
  next unless $ct == 1000;
  Blast();
  @lines = ();
  $ct = 0;
 } else {chomp; $seq .= $_}
}
close IN;
push @lines, ">$head\n$seq\n" if $seq;
exit unless $lines[0];
Blast();
close OUT;

sub Blast {
 open TEST, ">$infile.paltest"; for (@lines) {print TEST $_} close TEST;
 system "makeblastdb -in $infile.paltest -out $infile.paltest -dbtype nucl > /dev/null";
 my @out = `blastn -perc_identity 95 -query $infile.paltest -db $infile.paltest -outfmt 6 | awk '(\$1==\$2)&&(\$9>\$10){print \$1}'`;
 system "rm $infile.paltest*";
 my %uniq; 
 for (@out) {$uniq{$_} ++}
 for (@out) {next if $uniq{$_} > 1; print OUT $_}
}
