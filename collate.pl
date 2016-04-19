#! /usr/bin/perl -w
use strict;

my (%excludes, %cts);
open IN, "/data1/users/kpwilli/mitomycin/exp5/excludes";
while (<IN>) {
 chomp;
 next unless /^\S+ (\S+) (\d+)-(\d+)/;
 $excludes{$1}{$2} = $3;
}
close IN;
#die scalar keys(%excludes), "\n";

open IN, $ARGV[0];
while (<IN>) {
 #pCuAs	1	115	1
 chomp;
 my @f = split "\t";
 next unless $f[3] == 1;
 my $flag;
 if ($excludes{$f[0]}) {
  for (keys %{$excludes{$f[0]}}) {
   next unless $f[1] >= $_;
   $flag ++ if $f[1] <= $excludes{$f[0]}{$_}
  }
 }
 next if $flag;
 push @{$cts{$f[0]}}, $f[2];
}
close IN;

for (sort keys %cts) {
 my ($n, $mean, $sd) = NMeanSD(@{$cts{$_}});
 print "$_\t$n\t$mean\t$sd\n";
}


sub NMeanSD { # returns mean and SD for list of values
 my $n = @_;
 if ($n < 1) {return 0, 0, 0}
 my $mean;
 for (@_) {$mean += $_}
 $mean /= $n;
 if ($n <= 1) {return $mean , 0}
 my $sumsquares = 0;
 for (@_) {$sumsquares += ($_ - $mean) ** 2}
 return $n, $mean , sqrt($sumsquares/($n - 1)) ;
}

