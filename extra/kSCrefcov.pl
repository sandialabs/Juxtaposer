#! /usr/bin/perl -w
use strict;

#die "Usage: $0 k repliconFasta\n" unless scalar @ARGV >= 2;
die "Usage: $0 k repliconFasta readCovFile\n" unless scalar @ARGV >= 3;
my ($k, $replfile, $covfile) = @ARGV;
die "Use odd k instead of $k\n" unless $k % 2;

my (%mers, %dnas, %repl, %copy, $dna, %revcomps);
sub Revcomp {my $ret = reverse $_[0]; $ret =~ tr/ACGT/TGCA/; return $ret}

open IN, $replfile;
while (<IN>) {
 if (/^>(\S+)/) {$dna = $1; next}
 chomp;
 $repl{$dna} .= uc $_;
}
close IN;

for my $dna (sort keys %repl) {
 $repl{$dna} .= substr($repl{$dna}, 0, $k-1); # Permute circles to get circle-closing kmers
 for ($repl{$dna} =~ /(?=(.{$k}))/g) {
  my $lower = (sort $_, Revcomp($_))[0];
  $copy{$lower} ++;
  $mers{$lower} = $dna;
 }
}

for my $mer (sort keys %mers) {
 next unless $copy{$mer} == 1;
 $dnas{$mers{$mer}} ++;
}
print "dnas"; for (sort keys %dnas) {print "\t$_"} print "\n";
print "uniq$k.mers"; for (sort keys %dnas) {print "\t", $dnas{$_}/2} print "\n";

#for my $file (glob "exp[0-6]/*/qf.fq.fa.cov") {
for my $file ($covfile) {
 my %cts;
 my @max = (0);
 my %hits;
 open IN, $file;
 while (<IN>) {
  chomp;
  next unless /^(\S+)\t\S+\t(\d+)/;
  next unless $copy{$1} and $copy{$1} == 1;
  $hits{$1} ++;
  push @{$cts{$mers{$1}}}, $2;
  @max = ($2, $1) if $2 > $max[0];
 }
 close IN;
 for my $mer (sort keys %mers) {
  next unless $copy{$mer} == 1;
  next if $hits{$mer};
  push @{$cts{$mers{$mer}}}, 0;
 }
 print $file;
 for my $dna (sort keys %dnas) {
  @{$cts{$dna}} = sort {$a <=> $b} @{$cts{$dna}};
  my $n = scalar @{$cts{$dna}};
  my ($sum, $med, $min, $max, $sd) = (0, $cts{$dna}[int($n/2)], $cts{$dna}[0], $cts{$dna}[-1], 0);
  for (@{$cts{$dna}}) {$sum += $_}
  my $avg = $sum/$n;
  for (@{$cts{$dna}}) {$sd += ($_ - $avg)**2};
  $sd = ($sd/$n)**0.5;
  print "\t", join("\t", $min, $med, $max, $avg, $sd);
 }
 print "\t$max[0]\t$max[1]\n";
}

