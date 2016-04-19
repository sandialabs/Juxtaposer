#! /usr/bin/perl -w
use strict;

unless (@ARGV == 3) {die "Usage: $0 comma-sep'd-readfiles ref.uniq21mers k\n"}
my (%mers, %mers2, %revs, %refmers, %pos, %cts, $dna, $i);
my ($files, $reffile, $k) = (@ARGV);

open IN, $reffile or die "No ref uniq21mers file $reffile\n";
while (<IN>) {
 if (/^#(\S+)/) {$dna = $1; $i = 0; next}
 chomp;
 my $rev = Revcomp($_);
 push @{$pos{$dna}}, [$_, $rev];
 $i ++;
 next unless $_;
 @{$refmers{$_}} = ($dna, $i);
 @{$refmers{$rev}} = ($dna, $i);
}
close IN;

for (split ',', $files) {
 if (/\.gz$/) {open IN, "zcat $_ |"}
 else {open IN, $_ or die "No file $_\n"}
 while (<IN>) {
  next if /^>/;
  for (/(?=(.{$k}))/g) {$mers{$_} ++ if $refmers{$_}}
 }
 close IN;
}

for my $dna (sort {scalar @{$pos{$b}} <=> scalar @{$pos{$a}}} keys %pos) {
 print "#$dna\n";
 for my $i (0 .. $#{$pos{$dna}}) {
  if ($pos{$dna}[$i][0] eq '') {print "\n"; next}
  for (@{$mers{$pos{$dna}[$i]}}) {$_ = 0 unless $_}
  print $mers{$pos{$dna}[$i][0]} + $mers{$pos{$dna}[$i][1]}, "\n";
 }
}

sub Revcomp {my $ret = reverse $_[0]; $ret =~ tr/ACGT/TGCA/; return $ret}
