#! /usr/bin/perl -w
use strict;

die "Usage: $0 .covFile .faFile outFile\n" unless @ARGV == 3;
my (%kmers, $dna, %seqs, %refmers, %out);
my $circ = 1;
my $k = 21;

open IN, $ARGV[0] or die "No .covFile $ARGV[0]\n";
while (<IN>) {
 next unless /(\S+)\t(\S+)\t(\d+)/;
 $kmers{$1} = $3;
 $kmers{$2} = $3;
}
close IN;

open IN, $ARGV[1] or die "No .faFile $ARGV[1]\n";
while (<IN>) {
 if (/^>(\S+)/) {$dna = $1; next}
 chomp;
 $seqs{$dna} .= $_;
}
close IN;

for my $dna (sort {length $seqs{$b} <=> length $seqs{$a}} keys %seqs) {
 my $seq = uc $seqs{$dna};
 my $front = substr $seq, 0, ($k-1)/2;
 my $back = substr $seq, (1-$k)/2;
 $seq = $back . $seq . $front; # permute circle
 my $ct = 0;
 for my $mer ($seq =~ /(?=(.{$k}))/g) {
  $refmers{$mer} ++;
  my $cov = 0;
  $cov = $kmers{$mer} if $kmers{$mer};
  $ct ++;
  %{$out{$dna}{$ct}} = (cov => $cov, mer => $mer);
 }
}
open OUT, ">$ARGV[2]" or die "No outfile $ARGV[2]\n";
for my $dna (sort {length $seqs{$b} <=> length $seqs{$a}} keys %seqs) {
 for my $ct (sort {$a <=> $b} keys %{$out{$dna}}) {
  print OUT join("\t", $dna, $ct, $out{$dna}{$ct}{cov}, $refmers{$out{$dna}{$ct}{mer}}), "\n";
 }
}
close OUT;
