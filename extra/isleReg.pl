#! /usr/bin/perl -w
use strict;

my %centers = qw/TandL 599436 49R 1321377 38RybB 1804319 29S 1980411 23SapB 2298106 42GloC 2366912 40GuaA 3989971 37X 4147737 16Fis 4927117/;
my (%covs, %isles);
for (keys %centers) {
 @{$isles{$_}} = ($centers{$_}-40000, $centers{$_}+39999);
}

open IN, $ARGV[0] or die "No infile $ARGV[0]\n";
while (<IN>) {
 #pCuAs	1	115	1
 chomp;
 my @f = split "\t";
 next unless $f[3] == 1;
 next unless $f[0] eq 'chromosome';
 my $flag;
 for (keys %isles) {
  next unless $f[1] >= $isles{$_}[0];
  $flag = $_ if $f[1] <= $isles{$_}[1];
 }
 next unless $flag;
 $covs{$flag}{$f[1] - $isles{$flag}[0]} = $f[2];
}
close IN;

print "coord"; for my $isle (sort {$centers{$a} <=> $centers{$b}} keys %centers) {print "\t$isle"} print "\n";
for my $i (0..79999) {
 print $i;
 for my $isle (sort {$centers{$a} <=> $centers{$b}} keys %centers) {
  if (defined $covs{$isle}{$i}) {print "\t$covs{$isle}{$i}"} else {print "\t"}
 }
 print "\n";
}
