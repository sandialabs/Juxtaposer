#! /usr/bin/perl -w
use strict;

my (%mobs, %pos, $dna);

my %sizes = (qw/chromosome 5435746 pNDM-US 140825 pCuAs 117755 pHg 85163 pMYS 2014/);
open IN, "mobpref.bed"; 
while (<IN>) {next if s/^#//; chomp; my @f = split "\t"; push @{$mobs{$f[0]}}, {name => $f[3], L => $f[1]+1, R => $f[2]}}
close IN; # Mobility genes

print "#sample";
for my $dna (sort {$sizes{$b} <=> $sizes{$a}} keys %sizes) {
 for my $mob (@{$mobs{$dna}}) {
  my $mobname = "$$mob{name}.$dna.$$mob{L}.$$mob{R}";
  for (qw/island Lflank Rflank/) {
   print "\tn_$_\_$mobname\tmean_$_\_$mobname\tmed_$_\_$mobname\tsd_$_\_$mobname";
  }
 }
}
print "\n";

my @files = glob "exp*/*/qf.fq.fa.uniq21mers";
#die scalar(@files), " files\n";
for my $file (@files) {
 print $file;
 %pos = ();
 open IN, $file or die "No uniq21mers data file $file\n"; while (<IN>) {if (/^#(\S+)/) {$dna = $1} else {chomp; push @{$pos{$dna}}, $_}} close IN;
 for my $dna (sort {$sizes{$b} <=> $sizes{$a}} keys %sizes) {
  #next unless $dna eq 'pHg';
  for my $mob (@{$mobs{$dna}}) {
   #next unless $$mob{name} eq 'IS26'; #print "$dna $$mob{L} $$mob{R}\n"; system "collectSeq.pl -i FINAL9.fa -e pHg -L $$mob{L} -R $$mob{R}";
   my ($L, $R, $len, @vals) = ($$mob{L}+10-1, $$mob{R}-10-1, $$mob{R}-$$mob{L}+1, ());
   for ($L .. $R) {push @vals, $pos{$dna}[$_] unless $pos{$dna}[$_] eq ''}
   my ($n, $mean, $med, $stdev) = MeanMedSd(@vals);
   print "\t$n\t$mean\t$med\t$stdev";
   #print "$dna\t$$mob{L}\t$$mob{R}\t$$mob{name}\t$len\t", scalar(@vals);
   ($L, $R, @vals) = ($$mob{L}-$len+10-1, $$mob{R}-$len-10-1, ());
   for ($L .. $R) {push @vals, $pos{$dna}[$_] unless $pos{$dna}[$_] eq ''}
   ($n, $mean, $med, $stdev) = MeanMedSd(@vals);
   print "\t$n\t$mean\t$med\t$stdev";
   #print "\t", scalar(@vals);
   ($L, $R, @vals) = ($$mob{L}+$len+10-1, $$mob{R}+$len-10-1, ());
   ($L, $R) = (1+10-1, $len-10-1) if $R > $sizes{$dna};
   for ($L .. $R) {push @vals, $pos{$dna}[$_] unless $pos{$dna}[$_] eq ''}
   ($n, $mean, $med, $stdev) = MeanMedSd(@vals);
   print "\t$n\t$mean\t$med\t$stdev";
   #print "\t", scalar(@vals), "\n";
  }
 }
 print "\n";
}

sub MeanMedSd {
 my ($mean, $stdev, $med, $n, @vals) = (0, 0, 0, scalar(@_), @_);
 if ($n == 0) {return (0, 0, 0, 0)}
 #for my $i (0..$#vals) {die "$i of $n\n" unless defined $vals[$i]} die "OK $n\n";
 for (@vals) {$mean += $_}               $mean /= $n;
 for (@vals) {$stdev += ($mean - $_)**2} $stdev = ($stdev/$n)**0.5;
 @vals = sort {$a <=> $b} @vals; 
 $med = @vals[int($n/2)];
 #return (sprintf("%.0f", $mean), sprintf("%.0f", $stdev));
 return ($n, $mean, $med, $stdev);
}

