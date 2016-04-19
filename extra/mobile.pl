#! /usr/bin/perl -w
use strict;

#die "mobiles.txt file already exists\n" if -f "mobiles.txt";
my $juxfile = "juxtas.txt";
my $mobilefile = "../../mobile.bed";

open OUT, ">juxtas.bed";
my (@juxtas, %overlaps);
open IN, $juxfile;
while (<IN>) {
 # scar	1317	chromosome	924198	925514	1	NS500534:23:H57HJBGXX:2:23208:7506:17567_94	181.3	ggCGGCCTGTTTCTCTTTTTCCGCGGCGGCCAGTCGGGTAGCxCAGCGttaagcgctttagcctggtcaagagcggcctgctgcgccttctccgcx	+	.	.
 chomp;
 my ($class, $len, $dna, $L, $R, $ct, $read, $sco, $seq, $sign, $loci, $other) = split "\t";
 push @juxtas, [$read, $_];
 @{$overlaps{$read}} = (10000000, 0);
 $R = $L if $R < $L; # Rare but necessary for bedtools to complete
 print OUT "$dna\t$L\t$R\t$read\n";
}
close IN;
close OUT;

my @hits = `bedtools intersect -wo -a juxtas.bed -b $mobilefile`;
for (@hits) {
 chomp;
 my @f = split "\t";
 my ($read, $mobile, $La, $Ra, $Lm, $Rm) = ($f[3], $f[7], $f[1], $f[2], $f[5], $f[6]);
 my ($dL, $dR, $pct, $lenM, $low) = ($Lm-$La, $Ra-$Rm, 100, $Rm-$Lm, $Lm-$La);
 if ($low > $dR) {$low = $dR}
 if ($dL < 0) {$pct += 100*$dL/$lenM}
 if ($dR < 0) {$pct += 100*$dR/$lenM}
 if ($low < 0) {if ($overlaps{$read}[1] < $pct) {@{$overlaps{$read}} = (0, $pct, $mobile, $Lm, $Rm)}}
 elsif ($low < $overlaps{$read}[0]) {@{$overlaps{$read}} = ($low, $pct, $mobile, $Lm, $Rm)}
}

for (@juxtas) {
 if ($overlaps{$$_[0]}[1]) {print join("\t", $$_[1], @{$overlaps{$$_[0]}}), "\n"}
 else {print "$$_[1]\n"}
}
