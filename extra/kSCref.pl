#! /usr/bin/perl -w
use strict;

die "Usage: $0 k repliconFasta\n" unless scalar @ARGV >= 2;
my ($k, $replfile) = @ARGV;
die "Use odd k instead of $k\n" unless $k % 2;

my (%mers, %dnas, %repl, %copy, $dna, %revcomps);

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
  for ($_, Revcomp($_)) {
   $copy{$_} ++;
   $mers{$_} = $dna;
  }
 }
}

for my $mer (sort keys %mers) {
 next unless $copy{$mer} == 1;
 $dnas{$mers{$mer}} ++;
 #print "$mer\t$revcomps{$mer}\t$mers{$mer}\n";
}
print "dnas"; for (sort keys %dnas) {print "\t$_"} print "\n";
print "uniq$k.mers"; for (sort keys %dnas) {print "\t", $dnas{$_}/2} print "\n";
exit;

for my $file (glob "exp[0-6]/*/qf.fq.fa") {
 my %cts;
 open IN, $file;
 while (<IN>) {
  next if /^>/;
  chomp;
  for (/(?=(.{$k}))/g) {
   next unless $copy{$_} and $copy{$_} == 1;
   $cts{$mers{$_}} ++;
  }
 }
 close IN;
 print $file; for (sort keys %dnas) {if ($cts{$_}) {print "\t$cts{$_}"} else {print "\t0"}} print "\n";
}

sub Revcomp {my $ret = reverse $_[0]; $ret =~ tr/ACGT/TGCA/; return $ret}
