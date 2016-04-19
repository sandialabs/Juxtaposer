#! /usr/bin/perl -w
use strict;

my @tiles;
while (<>) {
 chomp;
 my @f = split "\t";
 my ($L, $R) = ($f[8], $f[9]);
 if ($L > $R) {($L, $R) = ($R, $L)}
 for my $tile (@tiles) {
  next if $f[1] ne $$tile{dna} or $L > $$tile{R} or $R < $$tile{L};
  print "$_ overlaps $$tile{line}\n";
 }
 push @tiles, {dna => $f[1], L => $L, R => $R, line => $_};
 #die "$tiles[0]{dna} $tiles[0]{L} $tiles[0]{R} $tiles[0]{line}\n";
}
