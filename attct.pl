#! /usr/bin/perl -w
use strict;
# Not paying attention to partner reads

my (%queries, %reads, %regexs, @islands, $island, $att, %hits);
open IN, "/data1/users/kpwilli/mitomycin/atts.fa" or die "No atts.fa file\n";
while (<IN>) {
 chomp;
 if (/^>(\S+)\.([LRBP])$/) {$island = $1; $att = $2; next}
 s/-//g;
 #print "'$_'\n";
 $queries{$island}{$att}{seq} .= uc $_;
 #print "$island $att '", uc($_), "'\n";
}
close IN;
#die scalar(keys %queries), " att sets\n";

for my $island (sort keys %queries) {
 unless (keys %{$queries{$island}} == 4) {warn "Not all 4 atts given for island $island; skip\n"; next}
 push @islands, $island;
 for my $att (qw/L R B P/) {
  $regexs{"$island.$att.fwd"} = qr/$queries{$island}{$att}{seq}/;
  my $revcomp = Revcomp($queries{$island}{$att}{seq});
  $regexs{"$island.$att.rev"} = qr/$revcomp/;
 }
}

my $seq;
open IN, "$ARGV[0]" or die "No fastq file $ARGV[0]\n";
while (<IN>) {
 my $seq = <IN>;
 <IN>; <IN>;
 chomp $seq;
 #print "$seq\n";
 for my $att (keys %regexs) {
  #print "$att $regexs{$att}\n";
  next unless $seq =~ $regexs{$att};
  #print "hi\n";
  $att =~ /^(.*)\.(rev|fwd)$/;
  #print "$1 $2\n";
  $seq = Revcomp($seq) if $2 eq 'rev';
  push @{$hits{$1}}, $seq;
  #print "$hits{$1}[0]\n";
 }
 #print "yes\n" if $seq =~ /CAAGTAATCTTCGGCATAA/;
 #exit;
}
close IN;

open OUT, ">$ARGV[0].atts";
for my $island (@islands) {
 my %ct;
 for my $att (qw/P B L R/) {if ($hits{"$island.$att"}) {$ct{$att} = scalar(@{$hits{"$island.$att"}})} else {$ct{$att} = 0}}
 print OUT "# $island PBLR\t$ct{P}\t$ct{B}\t$ct{L}\t$ct{R}\n";
 for my	$att (qw/P B L R/) {if ($hits{"$island.$att"}) {print OUT join("\n", "# $island $att", @{$hits{"$island.$att"}}), "\n\n"}}
}
close OUT;

sub Revcomp {my $ret = reverse $_[0]; $ret =~ tr/ACGT/TGCA/; return $ret}
