#! /usr/bin/perl -w
use strict;

unless (@ARGV == 3) {die "Usage: $0 comma-sep'd-readfiles k outprefix\n"}
my (%mers, %mers2, %revs, %cts);
my ($files, $k, $outp) = (@ARGV);
for (split ',', $files) {
 if (/\.gz$/) {open IN, "zcat $_ |"}
 else {open IN, $_ or die "No file $_\n"}
 while (<IN>) {
  my $seq = <IN>;
  chomp $seq;
  for ($seq =~ /(?=(.{$k}))/g) {$mers{$_} ++}
  <IN>; <IN>;
 }
 close IN;
}
print scalar (keys %mers), " ungrouped mers\n";

for my $mer (keys %mers) {
 my $rev = Revcomp($mer);
 if (($rev cmp $mer) == -1) {$mers2{$rev} += $mers{$mer}; $revs{$rev}=$mer}
 else {$mers2{$mer} += $mers{$mer}; $revs{$mer} = $rev}
}
print scalar (keys %mers2), " diff mers\n";

for (keys %mers2) {$cts{$mers2{$_}} ++}

open COV, ">$outp.cov";
for (sort {$mers2{$a} <=> $mers2{$b}} keys %mers2) {print COV "$_\t$revs{$_}\t$mers2{$_}\n"}
close COV;

open PRF, ">$outp.profile";
for (sort {$a <=> $b} keys %cts) {print PRF "$_\t$cts{$_}\n"}
close PRF;

sub Revcomp {my $ret = reverse $_[0]; $ret =~ tr/ACGT/TGCA/; return $ret}
