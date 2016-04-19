#! /usr/bin/perl -w
use strict;
my @rawcats = qw/raw meanLen/;
my (%raw, %tposseqs, $ent, $seq);
for (`grep exp sumq`) {next unless /^(\S+)\t(\d+)\t\d+\t([0-9\.]+)/; %{$raw{$1}} = (raw => $2, meanLen => $3)}
#for (sort keys %raw) {print "$_ $raw{$_}{raw} $raw{$_}{meanLen}\n"} exit;
open IN, "tposonFams/all.mark"; 
while (<IN>) {
 chomp;
 if (/>(\S+)/) {$tposseqs{$seq} = $ent if $seq; $seq = ''; $ent = $1; next}
 $_ =~ s/-//g;
 $seq .= $_;
}
close IN;
$tposseqs{$seq} = $ent;
#print "$tposseqs{ACTTGACCACAACAGACTGTTGTGGTCAAAT}\n"; die scalar (keys %tposseqs), " tps\n";
#print "$tposseqs{ATTGAGCCTTGACACATGTTGTAATGGCTCAAT}\n"; die scalar (keys %tposseqs), " tps\n";
my @attcats = `grep -o -P '^>\\S+' atts.fa`;
for (@attcats) {chomp; s/^>//} # print "'$_'\n"}
my @regexcats = qw/IS26 IS4321circle1 IS5075circle1 ISKpn1 ISKpn18 ISKpn21/;
my @tposcats = qw/IS26.1 IS26.2 IS26.3 IS26.4 IS26.6 IS4321.1 IS5075.1 ISKpn1.1 ISKpn1.2 ISKpn1.35 ISKpn1.4 ISKpn18.1 ISKpn18.2 ISKpn21.c ISKpn21.p ISKpn21.x/;
my @ctscats = qw/reads nonstandard genomic multi double shift readspal nonstandardpal genomicpal multipal doublepal shiftpal/;
print "#Expt/Sample"; for (@rawcats, @ctscats, @attcats, @tposcats) {print "\t$_"} print "\n";
for my $e (glob "exp[0-9]") { # Experiments
 for my $s (`ls -l $e`) { # Samples
  next unless $s =~ /^d/;
  chomp $s; $s =~ s/.* /$e\//;
  chdir "$s";
  my %ct;
  for (@rawcats) {$ct{$_} = $raw{$s}{$_} if $raw{$s};}
  my $counts = `cat counts.txt`;
  for (@ctscats) {
   next unless $counts =~ /$_=(\d+)/;
   $ct{$_} = $1;
  }
  $ct{readspal} = `wc -l qf.fq.pal`; chomp $ct{readspal}; $ct{readspal} =~ s/ .*//;
  my @regs = `grep '0=' ct.log`;
  for (@regs) {
   next unless /^# (\S+)\t0=(\d+)/;
   $ct{$1} = $2;
  }
  open IN, "qf.fq.fa.hits";
  my ($n, $ent, $seq);
  while (<IN>) {
   chomp;
   if (/^>\d+\_(\d+)\_\d+\_(\S+)/) {
    #$ct{hitct} ++; $ct{hittot} += $1;
    my ($nnew, $entnew) = ($1, $2);
    if ($seq) {
     #warn "$ent $s $n\n" if $seq eq "ACTTGACCACAACAGACTGTTGTGGTCAAAT";
     #warn "seq $ent $n $seq?\n" unless $tposseqs{$seq};
     if ($ent =~ /circle2|circleS|\//) {($n, $ent, $seq) = ($nnew, $entnew, ''); next}
     warn "$ent $seq\n" unless $tposseqs{$seq};
     $ct{$tposseqs{$seq}} += $n; 
    }
    ($n, $ent, $seq) = ($nnew, $entnew, '');
    next;
   }
   $seq .= $_;
  } 
  close IN;
  $ct{$tposseqs{$seq}} += $n if $seq and $tposseqs{$seq};
  my @atts = `grep PBLR qf.fq.atts`;
  for (@atts) {
   /(\S+)\s+PBLR\t(\d+)\t(\d+)\t(\d+)\t(\d+)/;
   $ct{$1 . '.P'} = $2; $ct{$1 .'.B'} = $3; $ct{$1 .'.L'} = $4; $ct{$1 .'.R'} = $5;
  }
  print $s; 
  #for (@rawcats, @ctscats, @attcats, @tposcats, 'hitct', 'hittot') {if ($ct{$_}) {print "\t$ct{$_}"} else {print "\t0"}}
  for (@rawcats, @ctscats, @attcats, @tposcats) {if ($ct{$_}) {print "\t$ct{$_}"} else {print "\t0"}}
  print "\n";
  chdir "../../";
 }
#last
}
