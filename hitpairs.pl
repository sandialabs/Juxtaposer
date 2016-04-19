#! /usr/bin/perl -w
use strict;

# Output note: Blast hit 1 in upper case, hit 2 flanked by slashes
my (%dnas, %reads, %seqs, $head, %uids, %tns, %mobs);
my ($delta, $xslen, $endwindow, $uniqlen) = ('^', 20000000, 5, 8);
die "Usage: perl $0 referencePathPrefix maxReadLength\n" unless @ARGV == 2;
my ($refdir, $readmax) = @ARGV;

my $spanlen = $readmax * 2;
open IN, "$refdir.lens" or die "No $refdir.lens\n"; while (<IN>) {chomp; %{$dnas{$1}} = (linear => $_, len => $2) if s/^(\S+)\s+(\d+)\s*//} close IN; # Replicon length, circularity
#for (keys %dnas) {print "$_ $dnas{$_}{len}\n"} exit;
open IN, "nonstandard.fa"; while (<IN>) {chomp; if (/^>(\S+)/) {$head = $1; next} $seqs{$head} .= $_;} close IN; # Read sequences
open IN, "$refdir.mobile.bed"; while (<IN>) {chomp; my @f = split "\t"; %{$mobs{$f[3]}} = (dna => $f[0], L => $f[1]+1, R => $f[2])} close IN; # Mobility genes
if (-f "$refdir.tn.bed") {
 open IN, "$refdir.tn.bed"; # Transposon left and right genomic ends
 while (<IN>) {
  chomp; my @f = split "\t";
  $uids{$f[3]} ++;
  %{$tns{"$f[3].$uids{$f[3]}L"}} = (dna => $f[0], L => $f[1]+1, R => $f[1]+$endwindow, term => $f[1]+1, end => -1);
  %{$tns{"$f[3].$uids{$f[3]}R"}} = (dna => $f[0], L => $f[2]-$endwindow+1, R => $f[2], term => $f[2]  , end =>  1);
 }
 close IN;
}
my %cts = (nonstandard => scalar(keys %seqs), mobs => scalar(keys %mobs), tn => scalar(keys %tns));

# LOAD HITS INTO %reads
open IN, "nonstandard.blast" or die "No nonstandard.blast file\n";
while (<IN>) { # Group Blast hits from each read, by DNA & orientation
 chomp;
 my @f = split "\t";
 my ($read, $dna, $id, $qL, $qR, $sS, $sE, $sco, $btop, $dir, $clos) = (@f[0,1,2,6,7,8,9,11,12], 1, 0);
 my ($sL, $sR) = ($sS, $sE);
 if ($sL > $sR) {($sL, $sR) = ($sE, $sS); $dir = -1}
 if ($dna =~ s/_CLOSURE_(\d+)//) { # Hit to replicon closure entry, only if it spans origin
  next unless $sL <= $1 and $sR > $1; # Ignore closure hits that don't span circular origin
  $clos = 1;
  for ($sE, $sS, $sL, $sR) {$_ = ($_ - $1) % $dnas{$dna}{len}} # Transpose to full replicon coordinates, may yield 0
 }
 my ($seq, $len, $offset) = ($seqs{$read}, $qR-$qL+1, $qL-1);
 $seq = lc(substr($seq, 0, $offset)) . uc(substr($seq, $offset, $len)) . lc(substr($seq, $qR)); # Uppercase hit portion only of read seq
 push @{$reads{$read}}, {dna => $dna, qL => $qL, qR => $qR, sL => $sL, sR => $sR, id => $id, len => $len, seq => $seq, sco => $sco, sS => $sS, sE => $sE,
  dir => $dir, offset => $offset, btop => $btop, clos => $clos};
}
close IN;

# PER READ: FIND SHIFTED HIT PAIRS AND SELECT SINGLE BEST; TEST FOR TRANSPOSITION, CIRCLE/SCAR, PALINDROME
my %juxtas;
for my $read (keys %reads) {
 $cts{genomic} ++;
 my $hitct = @{$reads{$read}};
 next unless $hitct > 1;
 $cts{multi} ++;
 $cts{double} ++ if $hitct == 2;
 my @tophit = (0); # [0] top score of ref-matching hits, [1-many] hits with that top score
 my %closs; # Keys are DNAs, in case multiple closures hit
 for my $i (0 .. $#{$reads{$read}}) { # Collect closure hits, and keep top score
  my $hit = $reads{$read}[$i];
  if ($$hit{clos}) {push @{$closs{$$hit{dna}}}, $$hit{sE}, $$hit{sS}}
  if (not $tophit[0] or $$hit{sco} > $tophit[0]) {@tophit = ($$hit{sco}, $i)} # New top score
  elsif ($$hit{sco} < $tophit[0]) {next} # Less than top score
  else {push @tophit, $i} # Equals previous top score, add to list
 }
 if (%closs) { # RE-INDEX FOR ANY DNA WHOSE CLOSURE WAS HIT
  for my $i (0 .. $#{$reads{$read}}) { # Load all coordinates hit on any closure-hit DNA
   my $hit = $reads{$read}[$i];
   next if not $closs{$$hit{dna}} or $$hit{clos}; # Skip if non-closure-hit DNA, or if already loaded clos hit
   push @{$closs{$$hit{dna}}}, $$hit{sE}, $$hit{sS};
  }
  for my $dna (keys %closs) { # Find temporary new origin for each closure-hit DNA
   my @gap = (0, 0);
   my $prev = 0;
   for (sort {$a <=> $b} @{$closs{$dna}}) { # Find biggest gap among combined L & R hit ends
    if ($_ - $gap[0] <= $gap[1]) {$gap[0] = $_; next}
    @gap = ($_, $_ - $gap[0]);
   }
   @{$closs{$dna}} = ($gap[0] + int(($gap[1]-$gap[0])/2));
  }
  for my $i (0 .. $#{$reads{$read}}) { # Transpose each hit on a reindexed closure-hit DNA
   my $hit = $reads{$read}[$i];
   next unless $closs{$$hit{dna}};
   for (qw/sE sS sL sR/) {$$hit{$_} = ($$hit{$_} + $closs{$$hit{dna}}[0]) % $dnas{$$hit{dna}}{len} }
  }
 }

 my @hitpairs; # EXAMINE ALL READ'S SHIFTED HIT PAIRS; SELECT ONE
 $tophit[0] =~ s/^ +//;
 for my $i (sort {$reads{$read}[$tophit[$a]]{dna} cmp $reads{$read}[$tophit[$b]]{dna} ||
                  $reads{$read}[$tophit[$a]]{sL}  <=> $reads{$read}[$tophit[$b]]{sL}  } 1..$#tophit) {
  for my $j (sort {$reads{$read}[$a]{dna} cmp $reads{$read}[$b]{dna} || $reads{$read}[$a]{sL} <=> $reads{$read}[$b]{sL}} 0..$#{$reads{$read}}) {
   next if $j == $tophit[$i]; # Skip self
   my @gnmOrder = sort {$reads{$read}[$a]{qL} <=> $reads{$read}[$b]{qL}} $tophit[$i], $j; # Left hit on read first
   #my @gnmOrder = sort {$reads{$read}[$a]{dna} cmp $reads{$read}[$b]{dna} || $reads{$read}[$a]{sL} <=> $reads{$read}[$b]{sL}} $tophit[$i], $j;
   my @pair = ($reads{$read}[$gnmOrder[0]], $reads{$read}[$gnmOrder[1]]); # pair[0]=Lower on alphabetically concatenated genome; pair[1]=higher
   my ($p0, $p1) = ($pair[0], $pair[1]);
   next if ${$p0}{qL} > ${$p1}{qL}-$uniqlen or ${$p0}{qR} > ${$p1}{qR}-$uniqlen; # Reject non-shifted configuration
   my $sco = ${$p0}{sco}+${$p1}{sco};  # Sum bitscores
   next if $hitpairs[0] and $sco < $hitpairs[0]; # Reject pairs with lower score than previous
   my @shift;  # [0]=1 if pair[0] on left end of read, else =-1; [1] pair[0] jxn coord; [2] pair[1] jxn coord; [3] overlap; [4] shortest span over origin? [5] shortest span around circular replicon

   my $overlap = ${$p0}{qR} - ${$p1}{qL} + 1; # Hit overlap on read
   my ($goff0, $roff0, $goff1, $roff1) = (0, 0, 0, 0); # Corrections on genome coords and on reads
   if ($overlap > 0) { # Check that positive overlap is perfect otherwise reduce until perfect
    my ($btop, $last) = (${$p0}{btop}, 0);
    while ($btop) {
     if ($btop =~ s/([0-9]+)$//) {$last = $1; last if $last+$roff0 > $overlap}
     next unless $btop =~ s/(.)(.)$//;
     $roff0 += $last+1; $goff0 += $last+1; $last = 0;
     if    ($1 eq '-') {$roff0 --}
     elsif ($2 eq '-') {$goff0 --}
     last if $roff0 > $overlap;
    }
    ${$p0}{sE} -= $goff0 * ${$p0}{dir};
    ${$p0}{qR} -= $roff0;
    ($btop, $last) = (${$p1}{btop}, 0);
    while ($btop) {
     if ($btop =~ s/^([0-9]+)//) {$last = $1; last if $last+$roff1 > $overlap}
     next unless $btop =~ s/^(.)(.)//;
     $roff1 += $last+1; $goff1 += $last+1; $last = 0;
     if    ($1 eq '-') {$roff1 --}
     elsif ($2 eq '-') {$goff1 --}
     last if $roff1 > $overlap;
    }
    ${$p1}{sS} -= $goff1 * ${$p1}{dir};
    ${$p1}{qL} += $roff1;
    $overlap   -= $goff0 + $goff1;
   }
   for my $hit (0,1) {
    my $dna = ${$pair[$hit]}{dna};
    for (${$pair[$hit]}{sE}, ${$pair[$hit]}{sS}, ${$pair[$hit]}{sL}, ${$pair[$hit]}{sR}) {
     $_ -= $closs{$dna}[0] if $closs{$dna};  # Transpose closure hits back to original coordinate system
     #die "$dna\n" unless $dnas{$dna}{len};
     $_ = Register($_, $dnas{$dna}{len});
    }
   }
   @shift = (1, ${$p0}{sE}, ${$p1}{sS}, $overlap, 0, -1);
   my @comp = (${$p0}{dna} cmp ${$p1}{dna}, $shift[1] <=> $shift[2]);
   if ($comp[0] == 1 or ($comp[0] == 0 and $comp[1] == 1)) { # Lexically lower was second; switch so it's first
     @gnmOrder = reverse @gnmOrder;
     @pair = reverse @pair;
     ($p0, $p1) = ($pair[0], $pair[1]);
     @shift = (-1, $shift[2], $shift[1], $overlap, 0, -1);
   }
   if ($comp[0] == 0) {  # Same DNA
    $shift[4] = 1; $shift[5] = $shift[2]-$shift[1];
    my $otherway = $dnas{${$p0}{dna}}{len} - $shift[5]; # Is the shorter distance around the replicon's origin?
    if ($otherway < $shift[5]) {$shift[4] = -1; $shift[5] = $otherway}
   }

   my @config = ('+','+','',''); # [0] hit jxn orient 0; [1] ... 1; [2] concat [0].[1]; [3] circle/scar/pal call;     # Config describes first the lower genomic-position hit, then the higher; + means Left end of hit is distal to joint; - means Right end is distal
   $config[0] = '-' if    ${$p0}{dir}*$shift[0] == -1;
   $config[1] = '-' if -1*${$p1}{dir}*$shift[0] == -1;
   $config[2] = $config[0].$config[1];
   for ($config[0], $config[1]) {$_ = 1 * ($_ . 1)}
   my @tnend; # [0][0] tn1 id, [0][1] tn1 offset, [1][0] tn2 id, [1][1] tn2 offset 
   for my $hit (0, 1) {
    for my $tn (keys %tns) { # Hit1 overlapping known tn ends?
     next if ${$pair[$hit]}{sR} < $tns{$tn}{L} or ${$pair[$hit]}{sL} > $tns{$tn}{R} or $tns{$tn}{dna} ne ${$pair[$hit]}{dna} or
       $config[$hit] * $tns{$tn}{end} == -1; 
     push @{$tnend[$hit]}, [$tn, $config[$hit] * ($shift[$hit+1] - $tns{$tn}{term})];
    }
    if ($tnend[$hit]) {my @flag = sort {abs($$a[1]) <=> abs($$b[1])} @{$tnend[$hit]}; push @{$tnend[$hit]}, $flag[0][0], $flag[0][1]}
    else {push @{$tnend[$hit]}, '.', '.'}
   }
   if ($tnend[1][-1] ne '.' and $tnend[1][-1] > 1000) {die "$tnend[1][-2] $tnend[1][-1] $shift[2] $config[1] $tns{$tnend[1][-2]}{term}\n"} 
   my ($seq, $len, @align) = (${$p0}{seq}, $xslen);
   for ($pair[0], $pair[1]) {
    my @hitseq = split //, substr($$_{seq}, $$_{offset}, $$_{len});
    my ($btop, $offset) = ($$_{btop}, 0);
    #print "$read $btop $offset @shift\n";
    while ($btop) {
     $offset += $1 if $btop =~ s/^([0-9]+)//;
     next unless $btop =~ s/^(.)(.)//;
     if    ($1 eq '-') {splice @hitseq, $offset, 0, '-'}  # BTOP gives query first in pair
     elsif ($2 eq '-') {splice @hitseq, $offset, 1, $delta}
     else {$hitseq[$offset] = lc $2} # Base substitution
     $offset ++;
    }
    #die "$read $btop $offset @shift\n";
    push @align, join('', @hitseq);
   }
   if (${$p0}{dna} eq ${$p1}{dna}) { # Both circleJxn and scar have hits in same orientation, on same DNA
    $config[3] = 'circleJxn' if ($config[2] eq '-+' and $shift[4] ==  1) or ($config[2] eq '+-' and $shift[4] == -1);
    $config[3] = 'scar'   if ($config[2] eq '-+' and $shift[4] == -1) or ($config[2] eq '+-' and $shift[4] ==  1);
    $config[3] = 'palindrome' if ($config[2] eq '--' or $config[2] eq '++') and $shift[5] < 10000;
    $len = $shift[5];
   }
   my $mobct = 0; for (0,1) {if ($tnend[$_][-2] eq '.') {$mobct ++}}
   if (not $config[3] and $mobct == 1) {$config[3] = 'transpose'} # One and only one hit is to a transposon end
   $config[3] = '.' unless $config[3];
   $cts{$config[3]} ++;
   my $second = substr $config[2], 1, 1; $second =~ tr/\-\+/+-/; $config[2] = substr($config[2], 0, 1) . $second; # Revert config to the more literal pair of hit orientations
   my @order = (sort {$a <=> $b} ${$p0}{offset}, ${$p0}{offset}+${$p0}{len}, ${$p1}{offset}, ${$p1}{offset}+${$p1}{len});
   my ($jL, $jR) = ($order[1], $order[2]);
   my $jxn = (substr($seq, $jL - 5, 5) .'/'. substr($seq, $jL, $jR-$jL) .'/'. substr($seq, $jR, 5));
   $seq = substr($seq, 0, ${$p1}{qL}-1) . '-' . substr($seq, ${$p1}{qL}-1, ${$p1}{qR}-${$p1}{qL}+1) . '-' . substr($seq, ${$p1}{qR});
   if ($shift[0] == -1) {for ($seq, $jxn, @align) {$_ = Revcomp($_)}}
   my $ct = 0; $ct = $hitpairs[-1]+1 if $hitpairs[-1] and $sco == $hitpairs[0];
   if (not $hitpairs[0] or $sco > $hitpairs[0] or $len < $hitpairs[1]) {
    @hitpairs = ($sco, $len, $shift[1], $shift[2], $config[2], $seq, $tnend[0][-2], $tnend[0][-1], $tnend[1][-2], $tnend[1][-1], $shift[3],
     $align[0], $align[1], $jxn, $config[3], $shift[4], $gnmOrder[0], $gnmOrder[1], 0);
   } else {$hitpairs[-1] ++}
  }
 }
 next unless @hitpairs;
 $cts{shift} ++;
 $cts{overlap}{$hitpairs[10]} ++;

 $hitpairs[1] = -1 if $hitpairs[1] == $xslen; 
 my ($hit1, $hit2) = ($reads{$read}[$hitpairs[-3]], $reads{$read}[$hitpairs[-2]]);
 my $label = "$$hit1{dna}/$hitpairs[2]/$$hit2{dna}/$hitpairs[3]/$hitpairs[4]"; 
 
 if (not $juxtas{$label} or $juxtas{$label}{sco} < $hitpairs[0]) {
  my $ct = 1; if ($juxtas{$label}) {$ct += $juxtas{$label}{ct}}
  %{$juxtas{$label}} = (read => $read, hit1 => $hitpairs[-3], hit2 => $hitpairs[-2], readseq => $hitpairs[5], ct => $ct,
   sites => $hitpairs[-1], coord => $hitpairs[2], coordp => $hitpairs[3], config => $hitpairs[4], q1R => $hitpairs[7], q2L => $hitpairs[8],
   overlap => $hitpairs[10], sco => $hitpairs[0], mob1 => $hitpairs[6], mob1off => $hitpairs[7], mob2 => $hitpairs[8], mob2off => $hitpairs[9], 
   len => $hitpairs[1], hitseq1 => $hitpairs[11], hitseq2 => $hitpairs[12], jxn => $hitpairs[13], call => $hitpairs[14], origin => $hitpairs[15]);
 } else {$juxtas{$label}{ct} ++}
}

# ASSIGN MOBILOME TO CIRCLES/SCARS
open OUT, ">juxtas.bed";
my (%overlaps);
for my $label (keys %juxtas) {
 next unless $juxtas{$label}{call} =~ /circle|scar/;
 my ($dna, $L, $R) = ($reads{$juxtas{$label}{read}}[$juxtas{$label}{hit1}]{dna}, $juxtas{$label}{coord}-1, $juxtas{$label}{coordp});
 if ($juxtas{$label}{origin} == -1) {print OUT "$dna\t$L\t", $dnas{$dna}{len},"\t$label\n$dna\t0\t$R\t$label\n"} else {
  $R = $L if $R < $L; # Rare but necessary for bedtools to complete
  print OUT "$dna\t$L\t$R\t$label\n";
 }
}
close OUT;

my @hits = `bedtools intersect -wo -a juxtas.bed -b $refdir.mobile.bed`;
unlink 'juxtas.bed';

for (@hits) {
 chomp;
 my @f = split "\t";
 my ($label, $mobile, $Lj, $Rj, $Lm, $Rm) = ($f[3], $f[7], $f[1], $f[2], $f[5], $f[6]);
 my ($dL, $dR, $pct, $lenM, $low) = ($Lm-$Lj, $Rj-$Rm, 100, $Rm-$Lm, $Lm-$Lj);
 if ($low > $dR) {$low = $dR}
 if ($dL < 0) {$pct += 100*$dL/$lenM}
 if ($dR < 0) {$pct += 100*$dR/$lenM}
 if ($low < 0) {
  if (not $overlaps{$label} or $overlaps{$label}[1] < $pct) {@{$overlaps{$label}} = (0, $pct, $mobile)}
 } elsif (not $overlaps{$label} or $low < $overlaps{$label}[0]) {@{$overlaps{$label}} = ($low, $pct, $mobile)}
}

open OUT, ">juxtas.txt";
print OUT '#', join("\t", qw(dnaLo CoordLo dnaHi CoordHi config call overlap origin idLo idHi tnLo tnOffsetLo tnHi tnOffsetHi sampleRead overlapSeq readseq seqLo seqHi
 count sites Mob pctMob shortDistMob)), "\n";
for my $label (sort {$juxtas{$b}{ct} <=> $juxtas{$a}{ct}} keys %juxtas) {
 my $line = $juxtas{$label};
 my ($hit1, $hit2) = ($reads{$$line{read}}[$$line{hit1}], $reads{$$line{read}}[$$line{hit2}]);
 my $mobile = ''; $mobile = join("\t", reverse @{$overlaps{$label}}) if $overlaps{$label};
 print OUT join("\t", $$hit1{dna}, $$line{coord}, $$hit2{dna}, $$line{coordp}, $$line{config}, $$line{call}, $$line{overlap}, $$line{origin}, $$hit1{id}, $$hit2{id}, $$line{mob1},
  $$line{mob1off}, $$line{mob2}, $$line{mob2off}, $$line{read}, $$line{jxn}, $$line{readseq}, $$line{hitseq1}, $$line{hitseq2}, $$line{ct}, $$line{sites}, $mobile) , "\n"; 
}
close OUT;

for (qw/nonstandard genomic multi double shift circleJxn scar palindrome transpose/) {$cts{$_} = 0 unless $cts{$_}; print "$_=$cts{$_} "} print "tn=$cts{tn}\n";
if ($cts{overlap}) {for (sort {$a <=> $b} keys %{$cts{overlap}}) {print "overlap $_=$cts{overlap}{$_}\n"}}

sub Revcomp {my $ret = reverse $_[0]; $ret =~ tr/ACGTacgt/TGCAtgca/; return $ret}
sub Min {return (sort {$a <=> $b} @_)[0]}
sub Register {my ($ret, $len) = @_; $ret = $ret % $len; $ret += $len if $ret < 1; return $ret}
