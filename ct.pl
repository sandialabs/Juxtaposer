#! /usr/bin/perl -w
use strict;

my (%raw, %regexs, $entry);
my @rawcats = qw/raw meanLen/;
for (`grep exp sumq`) {next unless /^(\S+)\t(\d+)\t\d+\t([0-9\.]+)/; %{$raw{$1}} = (raw => $2, meanLen => $3)}
my @ctscats = qw/reads nonstandard genomic multi double shift readspal nonstandardpal genomicpal multipal doublepal shiftpal/;

my @files = glob "exp[0-6]/*/qf.fq.fa";
print scalar(@files), " files\n";

my %uniqs;
open IN, "/data1/users/kpwilli/mitomycin/regexall.txt" or die "No regexall.txt file\n";
while (<IN>) {
 chomp;
 if (/^>(\S+)/) {$entry = $1; print $entry unless $entry =~ /\.[BPLR]$/; next}
 $_ = uc $_;
 push @{$regexs{$entry}{fwd}}, qr/($_)/;
 my $revcomp = Revcomp($_);
 $revcomp =~ s/\*\./.*/g;
 print "$entry same as $uniqs{$_}\n" if $uniqs{$_};
 print "$entry same as $uniqs{$revcomp}\n" if $uniqs{$revcomp};
 $uniqs{$_} = $entry; $uniqs{$revcomp} = $entry;
 push @{$regexs{$entry}{rev}}, qr/($revcomp)/;
}
close IN;
exit;

for my $readfile (@files) {
#my $readfile = 'nonstandard.fa'; $readfile = $ARGV[0] if $ARGV[0];
my (%cts, %hits, %lens);
open IN, $readfile or die "No fasta file $readfile\n";
while (<IN>) {
 next if /^>/;
 $cts{reads}++;
 chomp $_;
 my %top;
 for my $re (keys %regexs) {
  for my $i (0..$#{$regexs{$re}{fwd}}) {
   $top{$1}{$re}          ++ if $_ =~ $regexs{$re}{fwd}[$i];
   $top{Revcomp($1)}{$re} ++ if $_ =~ $regexs{$re}{rev}[$i];
  }
 }
 next unless %top;
 #warn scalar(keys %top), " regions hit in read $_\n" unless scalar(keys %top) == 1;
 for my $seq (keys %top) {
  warn scalar(keys %{$top{$seq}}), " regexs hit $seq\n" unless scalar(keys %{$top{$seq}}) == 1;
  for my $re (keys %{$top{$seq}}) {$hits{$re}{$seq} ++;}
 }
}
close IN;
print "$cts{reads} reads\n";

open OUT, ">$readfile.hits";
for my $re (sort keys %regexs) {
 my %sizes;
 if ($hits{$re}) {for (keys %{$hits{$re}}) {$sizes{length($_)} += $hits{$re}{$_}; $sizes{00} += $hits{$re}{$_}}}
 print "# $re"; for (sort {$a <=> $b} keys %sizes) {print "\t$_=$sizes{$_}"} print "\n";
 my $ct = 0;
 if ($hits{$re}) {for (sort {$hits{$re}{$b} <=> $hits{$re}{$a} || $a cmp $b} keys %{$hits{$re}}) {$ct ++; print OUT ">$ct\_$hits{$re}{$_}\_", length($_), "_$re\n$_\n"}}
}
close OUT;
}

sub Revcomp {my $ret = reverse $_[0]; $ret =~ tr/ACGT/TGCA/; return $ret}
