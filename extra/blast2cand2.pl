#! /usr/bin/perl -w
use strict;

my ($min, $max) = (1000, 200000);
my (%lens, %reads, %cts, %seqs, $head);
open IN, "nonstandard.fa";
while (<IN>) {
 chomp;
 if (s/^>//) {$head = $_; next}
 $seqs{$head} .= $_;
}
close IN;

open IN, "nonstandard.blast";
while (<IN>) { # Group Blast hits from each read, by DNA & orientation
 chomp;
 my @f = split "\t";
 my ($read, $dna, $id, $qL, $qR, $sL, $sR, $sign) = (@f[0,1,2,6,7,8,9], '+');
 if ($id < 94) {next} 
 $read =~ /_(\d+)$/; 
 my $len = $1;
 $lens{$read} = $len;
 my $seq = $seqs{$read};
 $seq = lc(substr($seq, 0, $qL-1)) . uc(substr($seq, $qL-1, $qR-$qL+1)) . lc(substr($seq, $qR));
 if ($sL > $sR) {
  for ($qL, $qR) {$_ = $len -$_ +1} # Reorient query to sense of DNA hit
  ($qL, $qR, $sL, $sR, $sign) = ($qR, $qL, $sR, $sL, '-');
  $seq = Revcomp($seq);
 }
 my $readset = "$read\t$dna\t$sign";
 push @{$reads{$readset}}, {dna => $dna, qL =>$qL, qR =>$qR, sL =>$sL, sR =>$sR, sign => $sign, id => $id, len => $len, seq => $seq};
 #die "$_\n$seq\n$seqs{$read}\n", length $seqs{$read};
}
close IN;

for my $readset (keys %reads) {
 @{$reads{$readset}} = sort {$$a{sL} <=> $$b{sL}} @{$reads{$readset}};
 my $hits = $reads{$readset};
 my $ct = scalar @{$reads{$readset}};
 for my $i (0 .. $ct - 2) {
  for my $j ($i + 1 .. $ct - 1){
   my ($class, $L, $R);
   if ($$hits[$i]{qL} < $$hits[$j]{qL} and $$hits[$i]{qR} < $$hits[$j]{qR}) {($class, $L, $R) = ('scar', $$hits[$i]{sR}, $$hits[$j]{sL})}
   if ($$hits[$i]{qL} > $$hits[$j]{qL} and $$hits[$i]{qR} > $$hits[$j]{qR}) {($class, $L, $R) = ('circle', $$hits[$i]{sL}, $$hits[$j]{sR})}
   next unless $class;
   my $seq = $$hits[$i]{seq};
   $seq = substr($seq, 0, $$hits[$j]{qL}-1) . 'x' . substr($seq, $$hits[$j]{qL}-1, $$hits[$j]{qR}-$$hits[$j]{qL}+1) . 'x' . substr($seq, $$hits[$j]{qR});
   $cts{$class}{$$hits[$i]{dna}}{$L}{$R}{$ct} ++;
   $cts{$class}{$$hits[$i]{dna}}{$L}{$R}{tot} ++;
   $readset =~ /^(\S+)\t.*(.)$/;
   push @{$cts{$class}{$$hits[$i]{dna}}{$L}{$R}{seqs}}, "$1:$seq";
  }
 }
}

for my $class (qw/scar circle/) {
 for my $dna (sort {$a cmp $b} keys %{$cts{$class}}) {
  for my $L (sort {$a <=> $b} keys %{$cts{$class}{$dna}}) {
   for my $R (sort {$a <=> $b} keys %{$cts{$class}{$dna}{$L}}) {
    unless ($cts{$class}{$dna}{$L}{$R}{2}) {$cts{$class}{$dna}{$L}{$R}{2} = 0}
    my $size = $R - $L +1; #print "$size\n";
    next if $size < $min or $size > $max;
    print join("\t", $class, $size, $dna, $L, $R, $cts{$class}{$dna}{$L}{$R}{tot}, $cts{$class}{$dna}{$L}{$R}{2},
     (join(',', @{$cts{$class}{$dna}{$L}{$R}{seqs}}))), "\n";
    #die;
   }
  }
 }
}

sub Revcomp {my $ret = reverse $_[0]; $ret =~ tr/ACGTactg/TGCAtgca/; return $ret}
