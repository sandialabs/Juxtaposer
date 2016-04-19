#! /usr/bin/perl -w
use strict;

die "Usage $0 infile permutelength outfile\n" unless @ARGV == 3;
my ($infile, $permlen, $outfile) = @ARGV;
my ($head, $seq);

open IN, $infile or die "No infile $infile\n";
open OUT, ">$outfile" or die "Can't write outfile $outfile\n";
while (<IN>) {
 chomp;
 if (/^>/) {
  print OUT Perm() if $seq;
  $head = "$_ perm$permlen";
  $seq = '';
  next;
 } 
 $seq .= $_;
}
close IN;
print OUT Perm() if $seq;
close OUT;

sub Perm {
 my $append = '';
 while (length $append < $permlen) {$append .= $seq}
 $append = substr $append, 0, $permlen;
 return "$head\n$seq$append\n";
}
