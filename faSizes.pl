#! /usr/bin/perl -w
use strict;

my ($entry, $dna);
while (<>) {
 chomp;
 if ( /^>(.+)/ ) {
  if ($entry) { Load() }
  $entry = $1;
  $dna = '';
 } else {$dna .= $_}
}
if ($entry) { Load() }

sub Load {
 $dna =~ s/[^a-zA-Z]//g;
 print "$entry\t", length($dna), "\n";
}



