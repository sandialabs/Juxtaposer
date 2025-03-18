#! /usr/bin/perl -w
use strict;

die "Usage $0 prefix readlength\n" unless @ARGV == 2;
my ($prefix, $readlen) = @ARGV;
my ($head, $seq);

open IN, "$prefix.fa" or die "No infile $prefix.fa\n";
open OUT, ">${prefix}_clos_$readlen.fa" or die "Can't write outfile ${prefix}_clos_$readlen.fa\n";
open GNM, ">$prefix.gnm" or die "Can't write genomefile $prefix.gnm\n";
while (<IN>) {
 chomp;
 if (/^>/) {
  print OUT Jxn() if $seq;
  print GNM "$head\t", length($seq), "\n" if $seq;
  $head = $_;
  $seq = '';
  next;
 } 
 $seq .= $_;
}
close IN;
print OUT Jxn() if $seq;
print GNM "$head\t", length($seq), "\n" if $seq;
close OUT;
close GNM;

sub Jxn {
 my $headclos = $head; $headclos =~ s/^(>\S+)/${1}_CLOSURE_$readlen/;
 return "$head\n$seq\n$headclos\n" . substr($seq, -1 * $readlen) . substr($seq, 0, $readlen) . "\n";
}
