#! /usr/bin/perl -w
use strict;

#my $oddball = 'Kleb_2.2.2688992_140';
# Output note: Blast hit 1 in upper case, hit 2 flanked by dashes
my (%lens, %reads, %cts, %seqs, $head, %pals, %refs, %uids);
my ($delta, $xslen, $endwindow, $uniqlen) = ('^', 20000000, 5, 8);
my $mobilefile = "/data1/users/kpwilli/mitomycin/mobile.bed";
open IN, $mobilefile;
while (<IN>) {
 chomp; my @f = split "\t";
 $uids{$f[3]} ++;
 %{$refs{"$f[3].$uids{$f[3]}L"}} = (dna => $f[0], L => $f[1]+1, R => $f[1]+$endwindow, term => $f[1]+1, end => -1);
 %{$refs{"$f[3].$uids{$f[3]}R"}} = (dna => $f[0], L => $f[2]-$endwindow+1, R => $f[2], term => $f[2]  , end =>  1);
}
close IN;
$cts{ref} = scalar(keys %refs);

open IN, "nonstandard.fa";
while (<IN>) {
 chomp;
 if (/^>(\S+)/) {$head = $1; $cts{nonstandard} ++; next}
 $seqs{$head} .= $_;
}
close IN;

open IN, "nonstandard.fa.pal";
while (<IN>) {
 chomp;
 $pals{$_} ++;
 $cts{nonstandardpal} ++;
}
close IN;

open IN, "nonstandard.blast";
while (<IN>) { # Group Blast hits from each read, by DNA & orientation
 chomp;
 my @f = split "\t";
 my ($read, $dna, $id, $qL, $qR, $sS, $sE, $sco, $btop, $dir) = (@f[0,1,2,6,7,8,9,11,12], 1);
 #$read =~ /_(\d+)$/; 
 #my $len = $1;
 #$lens{$read} = $len;
 my $seq = $seqs{$read};
 $seq = lc(substr($seq, 0, $qL-1)) . uc(substr($seq, $qL-1, $qR-$qL+1)) . lc(substr($seq, $qR));
 my ($sL, $sR) = ($sS, $sE);
 if ($sL > $sR) {($sL, $sR) = ($sE, $sS); $dir = -1}
 push @{$reads{$read}}, {dna => $dna, qL => $qL, qR => $qR, sL => $sL, sR => $sR, id => $id, len => $qR-$qL+1, #$len,
  seq => $seq, sco => $sco, sS => $sS, sE => $sE, dir => $dir, offset => $qL-1, btop => $btop};
}
close IN;

my %juxtas;
for my $read (keys %reads) {
 $cts{genomic} ++;
 $cts{genomicpal} ++ if $pals{$read};
 my $hitct = @{$reads{$read}};
 next unless $hitct > 1;
 $cts{multi} ++;
 $cts{multipal} ++ if $pals{$read};
 $cts{double} ++ if $hitct == 2;
 $cts{doublepal} ++ if $pals{$read} and $hitct == 2;
 my @tophit = (0); # [0] top score of ref-matching hits, [1-many] hits with that top score
 for my $i (0 .. $#{$reads{$read}}) {
  my $hit = $reads{$read}[$i];
  if (not $tophit[0] or $$hit{sco} > $tophit[0]) {@tophit = ($$hit{sco}, $i)} # New top score
  elsif ($$hit{sco} < $tophit[0]) {next} # Less than top score
  else {push @tophit, $i} # Equals previous top score, add to list
 }
 $tophit[0] =~ s/^ +//;
 my @hitpairs;
 for my $i (sort {$reads{$read}[$tophit[$a]]{dna} cmp $reads{$read}[$tophit[$b]]{dna} || $reads{$read}[$tophit[$a]]{sL} <=> $reads{$read}[$tophit[$b]]{sL}} 1..$#tophit) {
  for my $j (sort {$reads{$read}[$a]{dna} cmp $reads{$read}[$b]{dna} || $reads{$read}[$a]{sL} <=> $reads{$read}[$b]{sL}} 0..$#{$reads{$read}}) {
   next if $j == $tophit[$i];
   my @gnmOrder = sort {$reads{$read}[$a]{dna} cmp $reads{$read}[$b]{dna} || $reads{$read}[$a]{sL} <=> $reads{$read}[$b]{sL}} $tophit[$i], $j;
   my ($hit1, $hit2) = ($reads{$read}[$gnmOrder[0]], $reads{$read}[$gnmOrder[1]]); # hit1=Lower on alphabetically concatenated genome; hit2=higher
   my @shift;  # $shift[0]=1 if hit1 on left end of read, else =-1
   my $readlen = length $$hit1{seq};
   @shift = ( 1, $$hit1{sE}, $$hit2{sS}) if $$hit1{qL} <= $$hit2{qL}-$uniqlen and $$hit2{qR} >= $$hit1{qR}+$uniqlen; # hit1 on left
   @shift = (-1, $$hit1{sS}, $$hit2{sE}) if $$hit2{qL} <= $$hit1{qL}-$uniqlen and $$hit1{qR} >= $$hit2{qR}+$uniqlen;
   next unless @shift;

   push @shift, $$hit1{qR} - $$hit2{qL} + 1; # Hit overlap on read
   $shift[3] = $$hit2{qR} - $$hit1{qL} + 1 if $shift[0] == -1;
   #die "$shift[0] $shift[1] $shift[2] $shift[3] $tophit[$i] $j $$hit1{dna} $$hit2{dna}" if $read eq 'Kleb_1.1.2109411_127';
   my $sco = $$hit1{sco}+$$hit2{sco};  # Sum bitscores
   if ($hitpairs[0] and $sco < $hitpairs[0]) {next}
   my @config = ('+','+'); $config[0] = '-' if $$hit1{dir}*$shift[0] == -1; $config[1] = '-' if -1*$$hit2{dir}*$shift[0] == -1; $config[2] = $config[0].$config[1];
   for ($config[0], $config[1]) {$_ = 1 * ($_ . 1)}
    # config describes first the lower genomic-position hit, then the higher; + means Left end of hit is distal to joint; - means Right end is distal
   my ($mob1, $moboffset1, $mob2, $moboffset2) = ('.', '.', '.', '.');
   my @flag = ();
   for my $ref (keys %refs) { # overlapping known mobiles?
    unless ($$hit1{sR} < $refs{$ref}{L} or $$hit1{sL} > $refs{$ref}{R} or $refs{$ref}{dna} ne $$hit1{dna} or $config[0] * $refs{$ref}{end} == -1) 
    {push @flag, [$ref, $config[0] * ($shift[1] - $refs{$ref}{term})]}
   }
   if (@flag) {@flag = sort {abs($$a[1]) <=> abs($$b[1])} @flag; ($mob1, $moboffset1) = ($flag[0][0], $flag[0][1])}
   @flag = ();
   for my $ref (keys %refs) {
    unless ($$hit2{sR} < $refs{$ref}{L} or $$hit2{sL} > $refs{$ref}{R} or $refs{$ref}{dna} ne $$hit2{dna} or $config[1] * $refs{$ref}{end} == -1) 
    {push @flag, [$ref, $config[1] * ($shift[2] - $refs{$ref}{term})]}
   }
   if (@flag) {@flag = sort {abs($$a[1]) <=> abs($$b[1])} @flag; ($mob2, $moboffset2) = ($flag[0][0], $flag[0][1])}
   if ($mob2 and $mob2 ne '.' and $moboffset2 > 1000) {die "$mob2 $moboffset2 $shift[2] $config[1] $refs{$mob2}{term}\n"} 
   my ($seq, $len, @align) = ($$hit1{seq}, $xslen);
   for ($hit1, $hit2) {
    my @hitseq = split //, substr($$_{seq}, $$_{offset}, $$_{len});
    my ($btop, $offset) = ($$_{btop}, 0);
    while ($btop) {
     #print "$offset $btop\n" if $read eq $oddball;
     $offset += $1 if $btop =~ s/^([0-9]+)//;
     next unless $btop =~ s/^(.)(.)//;
     #print " $1 $2\n" if $read eq $oddball;
     if    ($1 eq '-') {splice @hitseq, $offset, 0, '-'}  # BTOP gives query first in pair
     elsif ($2 eq '-') {splice @hitseq, $offset, 1, $delta}
     else {$hitseq[$offset] = lc $2}
     $offset ++;
     #print join('', @hitseq), "  $offset $btop\n" if $read eq $oddball;
    }
    push @align, join('', @hitseq);
   }
   if ($$hit1{dna} eq $$hit2{dna}) { # both of the following have hits in same orientation, on same DNA
    $config[2] = 'circle' if $config[2] eq '-+';
    $config[2] = 'scar' if $config[2] eq '+-';
    $len = $shift[2]-$shift[1];
   }
   #my ($jL, $jR) = ((sort {$a <=> $b} $$hit1{qL}, $$hit1{qR}, $$hit2{qL}, $$hit2{qR})[1], (sort {$a <=> $b} $$hit1{qL}, $$hit1{qR}, $$hit2{qL}, $$hit2{qR})[2]);
   my $jL = (sort {$a <=> $b} $$hit1{offset}, $$hit1{offset}+$$hit1{len}, $$hit2{offset}, $$hit2{offset}+$$hit2{len})[1];
   my $jR = (sort {$a <=> $b} $$hit1{offset}, $$hit1{offset}+$$hit1{len}, $$hit2{offset}, $$hit2{offset}+$$hit2{len})[2];
   my $jxn = (substr($seq, $jL - 5, 5) .'/'. substr($seq, $jL, $jR-$jL) .'/'. substr($seq, $jR, 5));
   $seq = substr($seq, 0, $$hit2{qL}-1) . '-' . substr($seq, $$hit2{qL}-1, $$hit2{qR}-$$hit2{qL}+1) . '-' . substr($seq, $$hit2{qR});
   if ($shift[0] == -1) {for ($seq, $jxn, @align) {$_ = Revcomp($_)}}
   my $ct = 0; $ct = $hitpairs[-1]+1 if $hitpairs[-1] and $sco == $hitpairs[0];
   if (not $hitpairs[0] or $sco > $hitpairs[0] or $len < $hitpairs[1]) {
    @hitpairs = ($sco, $len, $shift[1], $shift[2], $config[2], $seq, $mob1, $moboffset1, $mob2, $moboffset2, $shift[3], $align[0], $align[1], $jxn, 
     $gnmOrder[0], $gnmOrder[1], 0);
   } else {$hitpairs[-1] ++}
  }
 }
 next unless @hitpairs;
 $cts{shift} ++;
 my $pal = '.'; if ($pals{$read}) {$pal = 'PAL'; $cts{shiftpal} ++}

 $hitpairs[1] = -1 if $hitpairs[1] == $xslen; 
 my ($hit1, $hit2) = ($reads{$read}[$hitpairs[-3]], $reads{$read}[$hitpairs[-2]]);
 my $label = "$$hit1{dna}/$hitpairs[2]/$$hit2{dna}/$hitpairs[3]/$hitpairs[4]"; 
 if (not $juxtas{$label} or $juxtas{$label}{sco} < $hitpairs[0]) {
  my $ct = 1; if ($juxtas{$label}) {$ct += $juxtas{$label}{ct}}
  %{$juxtas{$label}} = (read => $read, hit1 => $hitpairs[-3], hit2 => $hitpairs[-2], pal => $pal, readseq => $hitpairs[5], 
   ct => $ct, sites => $hitpairs[-1], coord => $hitpairs[2], coordp => $hitpairs[3], config => $hitpairs[4], q1R => $hitpairs[7], q2L => $hitpairs[8],
   overlap => $hitpairs[10], sco => $hitpairs[0], mob1 => $hitpairs[6], mob1off => $hitpairs[7], mob2 => $hitpairs[8], mob2off => $hitpairs[9], 
   len => $hitpairs[1], hitseq1 => $hitpairs[11], hitseq2 => $hitpairs[12], jxn => $hitpairs[13]);
 } else {$juxtas{$label}{ct} ++}
}

open OUT, ">juxtas.bed";
my (%overlaps);
for my $label (keys %juxtas) {
 next unless $juxtas{$label}{config} =~ /circle|scar/;
 my ($dna, $L, $R) = ($reads{$juxtas{$label}{read}}[$juxtas{$label}{hit1}]{dna}, $juxtas{$label}{coord}, $juxtas{$label}{coordp});
 $R = $L if $R < $L; # Rare but necessary for bedtools to complete
 print OUT "$dna\t$L\t$R\t$label\n";
}
close OUT;

my @hits = `bedtools intersect -wo -a juxtas.bed -b $mobilefile`;
for (@hits) {
 chomp;
 my @f = split "\t";
 my ($label, $mobile, $Lj, $Rj, $Lm, $Rm) = ($f[3], $f[7], $f[1], $f[2], $f[5], $f[6]);
 my ($dL, $dR, $pct, $lenM, $low) = ($Lm-$Lj, $Rj-$Rm, 100, $Rm-$Lm, $Lm-$Lj);
 if ($low > $dR) {$low = $dR}
 if ($dL < 0) {$pct += 100*$dL/$lenM}
 if ($dR < 0) {$pct += 100*$dR/$lenM}
 if ($low < 0) {
  if (not $overlaps{$label} or $overlaps{$label}[1] < $pct) {@{$overlaps{$label}} = (0, $pct, $mobile, $Lm, $Rm)}
 } elsif (not $overlaps{$label} or $low < $overlaps{$label}[0]) {@{$overlaps{$label}} = ($low, $pct, $mobile, $Lm, $Rm)}
}

open OUT, ">juxtas.txt";
print OUT '#', join("\t", qw(dnaLo jxnCoordLo dnaHi jxnCoordHi config pal overlap idLo idHi mobLo mobOffsetLo mobHi mobOffsetHi sampleRead jxn readseq seqLo seqHi
 count sites shortDistMob pctMob Mob MobL MobR)), "\n";
for my $label (sort {$juxtas{$b}{ct} <=> $juxtas{$a}{ct}} keys %juxtas) {
 my $line = $juxtas{$label};
 my ($hit1, $hit2) = ($reads{$$line{read}}[$$line{hit1}], $reads{$$line{read}}[$$line{hit2}]);
 my $mobile = ''; $mobile = join("\t", @{$overlaps{$label}}) if $overlaps{$label};
 print OUT join("\t", $$hit1{dna}, $$line{coord}, $$hit2{dna}, $$line{coordp}, $$line{config}, $$line{pal}, $$line{overlap}, $$hit1{id}, $$hit2{id}, $$line{mob1},
  $$line{mob1off}, $$line{mob2}, $$line{mob2off}, $$line{read}, $$line{jxn}, $$line{readseq}, $$line{hitseq1}, $$line{hitseq2}, $$line{ct}, $$line{sites}, $mobile), "\n"; 
}
close OUT;

$cts{reads} = `wc -l qf.fq`; $cts{reads} =~ s/ .*//; $cts{reads} /= 4;
if (-f "qf.fq.pal") {$cts{readspal} = `wc -l qf.fq.pal`; chomp $cts{readspal}; $cts{readspal} =~ s/ .*//;}
for my $pal ('', 'pal') {for (qw/reads nonstandard genomic multi double shift/) {$cts{$_.$pal} = 0 unless $cts{$_.$pal}; print "$_$pal=$cts{$_.$pal} "}} print "ref=$cts{ref}\n";

sub Revcomp {my $ret = reverse $_[0]; $ret =~ tr/ACGTacgt/TGCAtgca/; return $ret}
