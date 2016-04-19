#! /usr/bin/perl -w
use strict;

unless (@ARGV == 1) {die "Usage: $0 ref\n"}
my $k = 21;
my (%dnas, $dna, %mers, %merlocs, %revs, %cts);
my ($files) = (@ARGV);

for (split ',', $files) {
 if (/\.gz$/) {open IN, "zcat $_ |"}
 else {open IN, $_ or die "No file $_\n"}
 while (<IN>) {
  if (/^>(\S+)/) {$dna = $1; next}
  chomp;
  $dnas{$dna} .= uc $_;
 }
 close IN;
}

for my $dna (keys %dnas) {
 my $i = int($k/2);
 my $seq = substr($dnas{$dna}, -1*$i) . $dnas{$dna} . substr($dnas{$dna}, 0, $i);
 for ($seq =~ /(?=(.{$k}))/g) {my $lo = (sort $_, Revcomp($_))[0]; push @{$merlocs{$dna}}, $lo; $mers{$lo} ++}
}
#print scalar (keys %mers), " diff mers\n";

for (sort keys %merlocs) {
 print "#$_ ", scalar(@{$merlocs{$_}}), "\n";
 for (@{$merlocs{$_}}) { 
  if ($mers{$_} == 1) {print "$_\n"} else {print "\n"}
 }
}

sub Revcomp {my $ret = reverse $_[0]; $ret =~ tr/ACGT/TGCA/; return $ret}
