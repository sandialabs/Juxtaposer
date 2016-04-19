#! /usr/bin/perl -w
use strict;

# Output note: Blast hit 1 in upper case, hit 2 flanked by x's

die "Usage: $0 reflenfile\n" unless $ARGV[0];
my $lenfile = $ARGV[0]; 
die "No reflenfile $lenfile\n" unless -f $lenfile;
my ($min, $max) = (0, 20000000000);
#my ($min, $max) = (1000, 200000);
my (%lens, %reads, %cts, %seqs, $head, %pals, %reflens);
open IN, "nonstandard.fa";
while (<IN>) {
 chomp;
 if (s/^>//) {$head = $_; next}
 $seqs{$head} .= $_;
}
close IN;

open IN, "nonstandard.fa.pal"; while (<IN>) {chomp; $pals{$_} ++} close IN; #die scalar(keys %pals), " pals\n";
open IN, $lenfile; while (<IN>) {$reflens{$1} = $2 if /^(\S+)\s+(\S+)/} close IN; 

open IN, "nonstandard.blast";
while (<IN>) { # Group Blast hits from each read, by DNA & orientation
 chomp;
 my @f = split "\t";
 my ($read, $dna, $id, $qL, $qR, $sL, $sR, $sco, $sign) = (@f[0,1,2,6,7,8,9,11], '+');
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
 push @{$reads{$read}{$dna}{$sign}}, {qL => $qL, qR => $qR, sL => $sL, sR => $sR, id => $id, len => $len, seq => $seq, sco => $sco};
}
close IN;

for my $read (keys %reads) {
 my (@out, $circjunxn);
 for my $dna (keys %{$reads{$read}}) {
  for my $sign (keys %{$reads{$read}{$dna}}) {
   @{$reads{$read}{$dna}{$sign}} = sort {$$a{sL} <=> $$b{sL}} @{$reads{$read}{$dna}{$sign}};
   my $hits = $reads{$read}{$dna}{$sign};
   my $ct = scalar @{$reads{$read}{$dna}{$sign}};
   for my $i (0 .. $ct - 2) {
    for my $j ($i + 1 .. $ct - 1){
     my ($class, $L, $R);
     if ($$hits[$i]{qL} < $$hits[$j]{qL} and $$hits[$i]{qR} < $$hits[$j]{qR}) {($class, $L, $R) = ('scar', $$hits[$i]{sR}, $$hits[$j]{sL})}
     if ($$hits[$i]{qL} > $$hits[$j]{qL} and $$hits[$i]{qR} > $$hits[$j]{qR}) {($class, $L, $R) = ('circle', $$hits[$i]{sL}, $$hits[$j]{sR})}
     next unless $class;
     my $sco = $$hits[$i]{sco} + $$hits[$j]{sco};
     my $seq = $$hits[$i]{seq};
     $seq = substr($seq, 0, $$hits[$j]{qL}-1) . '-' . substr($seq, $$hits[$j]{qL}-1, $$hits[$j]{qR}-$$hits[$j]{qL}+1) . '-' . substr($seq, $$hits[$j]{qR});
     push @out, {class => $class, L => $L, R => $R, sco => $sco, seq => $seq, dna => $dna, sign => $sign, len => $R - $L + 1};
     if ($class eq 'circle' and (($$hits[$i]{sR} == $reflens{$dna} and $$hits[$j]{sL} == 1) or ($$hits[$j]{sR} == $reflens{$dna} and $$hits[$i]{sL} == 1))) {$circjunxn ++}
    }
   }
  }
 }
 next unless @out;
 next if $circjunxn; # Circular junctions not all removed by bowtie on permuted index!
 my @others = ();
 @out = sort {$$b{sco} <=> $$a{sco} || $$a{len} <=> $$b{len}} @out;
 for (1..$#out) {
  last unless $out[$_]{len} == $out[0]{len} and $out[$_]{sco} == $out[0]{sco};
  push @others, "$out[$_]{dna}:$out[$_]{L}..$out[$_]{R}$out[$_]{sign}";
 }
 my $other = '.'; if (@others) {$other = join(',', @others)}
 push @{$cts{$out[0]{class}}{$out[0]{dna}}{$out[0]{L}}{$out[0]{R}}}, [$read, $out[0]{sco}, $out[0]{seq}, $out[0]{sign}, $other];
}

for my $class (qw/scar circle/) {
 for my $dna (sort {$a cmp $b} keys %{$cts{$class}}) {
  for my $L (sort {$a <=> $b} keys %{$cts{$class}{$dna}}) {
   for my $R (sort {$a <=> $b} keys %{$cts{$class}{$dna}{$L}}) {
    @{$cts{$class}{$dna}{$L}{$R}} = sort {$$b[1] <=> $$a[1]} @{$cts{$class}{$dna}{$L}{$R}};
    my $len = $R - $L +1;
    next if $len < $min or $len > $max;
    my @others = ();
    for (1..$#{$cts{$class}{$dna}{$L}{$R}}) {
     my $read = $cts{$class}{$dna}{$L}{$R}[$_][0];
     $read .= 'PAL' if $pals{$read};
     push @others, $read;
    }
    $cts{$class}{$dna}{$L}{$R}[0][0] .= 'PAL' if $pals{$cts{$class}{$dna}{$L}{$R}[0][0]}; 
    my $other = '.'; if (@others) {$other = join(',', @others)}
    print join("\t", $class, $len, $dna, $L, $R, scalar(@{$cts{$class}{$dna}{$L}{$R}}), @{$cts{$class}{$dna}{$L}{$R}[0]}, $other), "\n";
    #die "$dna $L $R $cts{$class}{$dna}{$L}{$R}[0][0] $other\n" if $cts{$class}{$dna}{$L}{$R}[0][0] =~ /1.474036_79/;
   }
  }
 }
}

sub Revcomp {my $ret = reverse $_[0]; $ret =~ tr/ACGTacgt/TGCAtgca/; return $ret}
