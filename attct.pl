#! /usr/bin/perl
use strict;
use warnings;
# Not paying attention to partner reads

die "Usage: $0 atts-fastA-file read-fastQ-file(s)\n" unless @ARGV > 1;
my (%queries, %reads, %regexs, $island, $att, %hits);
open IN, $ARGV[0] or die "No atts.fa file\n";
while (<IN>) {
 chomp;
 if (/^>/) {
  die "Att.fa headers require an internal period, eg >EcoX.L\n" unless /^>(\S+)\.(\S+)$/;
  ($island, $att) = ($1, $2);
  next
 }
 s/-//g;
 $queries{$island}{$att}{seq} .= uc $_;
}
close IN;

for my $island (sort keys %queries) {
 for my $att (sort keys %{$queries{$island}}) {
  $regexs{"$island.$att.fwd"} = qr/$queries{$island}{$att}{seq}/;
  my $revcomp = Revcomp($queries{$island}{$att}{seq});
  $regexs{"$island.$att.rev"} = qr/$revcomp/;
 }
}

my $seq;
for (1..$#ARGV) {
 open IN, "$ARGV[$_]" or die "No fastq file $ARGV[$_]\n";
 while (<IN>) {
  my $seq = <IN>;
  <IN>; <IN>;
  chomp $seq;
  for my $att (keys %regexs) {
   next unless $seq =~ $regexs{$att};
   $att =~ /^(.*)\.(rev|fwd)$/;
   if ($2 eq 'fwd') {push @{$hits{$1}}, $seq} else {push @{$hits{$1}}, Revcomp($seq)}
  }
 }
 close IN;
}

for my $island (sort keys %queries) {
 my (%ct, @atts);
 for my $att (sort keys %{$queries{$island}}) {
  push @atts, $att;
  if ($hits{"$island.$att"}) {$ct{$att} = scalar(@{$hits{"$island.$att"}})}
  else {$ct{$att} = 0}
 }
 print "## $island ", join(',', @atts); for (@atts) {print "\t$ct{$_}"} print "\n";
 for my	$att (@atts) {if ($hits{"$island.$att"}) {print join("\n", "# $island $att", @{$hits{"$island.$att"}}), "\n\n"}}
}

sub Revcomp {my $ret = reverse $_[0]; $ret =~ tr/ACGT/TGCA/; return $ret}
