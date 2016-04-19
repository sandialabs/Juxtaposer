#! /usr/bin/perl -w
use strict;

my (%mobs, %pos, $dna, @samples);

my %sizes = (qw/chromosome 5435746 pNDM-US 140825 pCuAs 117755 pHg 85163 pMYS 2014/);
open IN, "mobpref.bed"; 
while (<IN>) {next if s/^#//; chomp; my @f = split "\t"; push @{$mobs{$f[0]}}, {name => $f[3], L => $f[1]+1, R => $f[2]}}
close IN; # Mobility genes

my %out;
my @files = glob "exp*/*/qf.fq.fa.uniq21mers";
for my $file (@files) {
 %pos = ();
 push @samples, $file; $samples[-1] =~ s/\/qf.fq.fa.uniq21mers//;
 open IN, $file or die "No uniq21mers data file $file\n"; while (<IN>) {if (/^#(\S+)/) {$dna = $1} else {chomp; push @{$pos{$dna}}, $_}} close IN;
 for my $dna (sort {$sizes{$b} <=> $sizes{$a}} keys %sizes) {
  for my $mob (@{$mobs{$dna}}) {
   my $len = $$mob{R}-$$mob{L};
   my ($L, $R) = ($$mob{L}-$len+10-1, $$mob{R}+$len-10-1);
   ($L, $R) = (1+10-1, $len-10-1) if $R > $sizes{$dna};
   for ($L .. $R) {push @{$out{$$mob{name}}{$file}}, $pos{$dna}[$_]}
  }
 }
}

for my $mob (keys %out) {
 open OUT, ">mobs/$mob.txt";
 print OUT join("\t", @samples), "\n";
 for my $i (0..$#{$out{$mob}{$files[0]}}) {
  for my $file (@files) {
   print OUT "$out{$mob}{$file}[$i]\t";
  }
  print OUT "\n";
 }
 close OUT;
}
