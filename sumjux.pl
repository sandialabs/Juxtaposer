#! /usr/bin/perl -w
use strict;

my (%jux, @samples);
my @folders = ('.');
#die "Usage: perl $0 <folders with juxtas.txt files>\n" unless @ARGV;
if (@ARGV) {for (@ARGV) {push @folders, $_}}

for my $folder (@folders) {
my @files = (glob("$folder/*juxtas.txt"), glob("$folder/*/*juxtas.txt"));
for my $sample (@files) {
 unless (open IN, $sample) {warn "Can't read $sample\n"; next}
 $sample =~ s/\/*juxtas.txt//;
 $sample =~ s/.*\///;
 push @samples, $sample;
 while (<IN>) {
  # scar    chromosome      18776   4099628 1       0       seq
  chomp;
  my @f = split "\t";
  my $name = join "\t", @f[0..3];
  #die "$folder $sample $name\n";
  %{$jux{$name}{$sample}} = (tot => $f[4], twohit => $f[5]);
  $jux{$name}{tot}{tot} += $f[4];
  $jux{$name}{tot}{twohit} += $f[5];
  #push @{$jux{$name}{eg}}, $f[6];
  $f[6] =~ s/,.*//;
  $jux{$name}{eg} = $f[6];
 }
 close IN;
}
}

#for (@samples) {print "$_\n"} exit;
die "No data, so no sumj.txt file written\n" unless %jux;
print "juxtaClass\tdna\tL\tR\tlength\tall tot\tall 2hit";
for (@samples) {print "\t$_ tot\t$_ 2hit"}
print "\teg\n";
for my $name (sort {$jux{$b}{tot}{tot} <=> $jux{$a}{tot}{tot} || $jux{$b}{tot}{twohit} <=> $jux{$a}{tot}{twohit}}keys %jux) {
 $name =~ /\t(\d+)\t(\d+)$/;
 my $len = $2 - $1 + 1;
 print "$name\t$len\t$jux{$name}{tot}{tot}\t$jux{$name}{tot}{twohit}";
 for (@samples) {print "\t", $jux{$name}{$_}{tot} // 0, "\t", $jux{$name}{$_}{twohit} // 0}
 #print "\t", join(',', @{$jux{$name}{eg}}), "\n";
 print "\t$jux{$name}{eg}";
 print "\n";
}
