#! /usr/bin/perl -w
use strict;
#exit;

# Output note: Blast hit 1 in upper case, hit 2 flanked by dashes
my (%lens, %reads, %cts, %seqs, $head);
my $uniqlen = 8;
open IN, "nonstandard.fa";
while (<IN>) {
 chomp;
 if (/^>(\S+)/) {$head = $1; next}
 $seqs{$head} .= $_;
}
close IN;

my (%refs, %uids);
my $endwindow = 5;
open IN, "/data1/users/kpwilli/mitomycin/mobile.bed";
while (<IN>) {
 chomp;
 my @f = split "\t";
 $uids{$f[3]} ++;
 %{$refs{"$f[3].$uids{$f[3]}L"}} = (dna => $f[0], L => $f[1]+1, R => $f[1]+$endwindow, end => $f[1]+1);
 %{$refs{"$f[3].$uids{$f[3]}R"}} = (dna => $f[0], L => $f[2]-$endwindow+1, R => $f[2], end => $f[2]);
}
close IN;
$cts{ref} = scalar(keys %refs);
#die "$cts{ref} refs\n";

open IN, "nonstandard.blast";
while (<IN>) { # Group Blast hits from each read, by DNA & orientation
 chomp;
 my @f = split "\t";
 my ($read, $dna, $id, $qL, $qR, $sS, $sE, $sco, $sign, $dir) = (@f[0,1,2,6,7,8,9,11], '+', 1);
 #if ($id < 94) {next} # Blast run now has perc_identity filter
 $read =~ /_(\d+)$/; 
 my $len = $1;
 $lens{$read} = $len;
 my $seq = $seqs{$read};
 $seq = lc(substr($seq, 0, $qL-1)) . uc(substr($seq, $qL-1, $qR-$qL+1)) . lc(substr($seq, $qR));
 my ($sL, $sR) = ($sS, $sE);
 if ($sL > $sR) {($sL, $sR) = ($sE, $sS); $sign = '-'; $dir = -1}
 push @{$reads{$read}}, {dna => $dna, qL => $qL, qR => $qR, sL => $sL, sR => $sR, id => $id, len => $len, seq => $seq, sco => $sco, sS => $sS, sE => $sE, sign => $sign, dir => $dir};
}
close IN;

for my $read (keys %reads) {
 $cts{read} ++;
 next unless @{$reads{$read}} > 1;
 $cts{multi} ++;
 my @flag = (0); # [0] top score of ref-matching hits, [1-many] hits with that top score
 for my $i (0 .. $#{$reads{$read}}) {
  my $hit = $reads{$read}[$i];
  for my $ref (keys %refs) {
   next if $$hit{sR} < $refs{$ref}{L} or $$hit{sL} > $refs{$ref}{R} or $refs{$ref}{dna} ne $$hit{dna};
   $$hit{ref} = $ref;
   unless ($flag[0]) {@flag = ($$hit{sco}, $i)}
   elsif ($$hit{sco} < $flag[0]) {next} # Less than top score
   elsif ($$hit{sco} > $flag[0]) {@flag = ($$hit{sco}, $i)} # New top score
   else {push @flag, $i} # Equals previous top score, add to list
  }
 }
 next unless $flag[0];
 $flag[0] =~ s/^ +//;
 #print "$read, @flag\n";
 $cts{hit} ++;
 my (@outs, @top);
 for my $i (1..$#flag) {
  my ($hit, $refori) = ($reads{$read}[$flag[$i]], 1);
  #print "$i $flag[$i] $$hit{ref}\n";
  $refori = -1 if $$hit{ref} =~ /L$/;
  for my $j (0 .. $#{$reads{$read}}) {
   next if $j == $flag[$i];
   #next if $seen{$j};
   my $hit2 = $reads{$read}[$j];
   next if $$hit{qR} >= $$hit2{qR}-$uniqlen and $$hit{dir} * $refori == 1;
   next if $$hit{qL} <= $$hit2{qL}+$uniqlen and $$hit{dir} * $refori == -1;
   my $sco = $$hit{sco}+$$hit2{sco};
   if (not $outs[0] or $sco > $outs[0]) {@outs = ($sco, $flag[$i], $j, 0)}
   elsif ($sco == $outs[0]) {$outs[3] ++}
   else {next}
   #print "$i, $flag[$i], $j, $sco, @outs\n";
   #print "$read $$hit{ref} $refori $$hit{qL} $$hit{qR} $$hit{dna} $$hit{sS} $$hit{sE} $$hit{dir}, $$hit2{qL} $$hit2{qR} $$hit2{dna} $$hit2{sS} $$hit2{sE} $$hit2{dir} $$hit2{seq}\n";
  }
 }
 next unless @outs;
 $cts{config} ++;
 my ($hit, $hit2) = ($reads{$read}[$outs[1]], $reads{$read}[$outs[2]]); 
 #$$hit{ref} .= $$hit{qR} - $refs{$$hit{ref}}{end}
 if ($$hit2{ref}) {$$hit2{ref} .= "/$refs{$$hit2{ref}}{end}"} else {$$hit2{ref} = '.'}
 $$hit{seq} = substr($$hit{seq}, 0, $$hit2{qL}-1).'-'.substr($$hit{seq}, $$hit2{qL}-1, $$hit2{qR}-$$hit2{qL}+1).'-'.substr($$hit{seq}, $$hit2{qR});
 print join("\t", $read, "$$hit{ref}/$refs{$$hit{ref}}{end}", "$$hit{qL}..$$hit{qR}/$$hit{id}", "$$hit{dna}/$$hit{sS}-$$hit{sE}",
  $$hit2{ref}, "$$hit2{qL}..$$hit2{qR}/$$hit2{id}", "$$hit2{dna}/$$hit2{sS}-$$hit2{sE}\t$outs[3]", $$hit{seq}), "\n";
 #last;
}

for (sort keys %cts) {warn "$_ $cts{$_}\n"}
