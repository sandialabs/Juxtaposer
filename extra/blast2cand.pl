#! /usr/bin/perl -w
use strict;

my (%lens, %reads, %cts);
while (<>) { # Group Blast hits from read by DNA & orientation
 chomp;
 my @f = split "\t";
 my ($read, $dna, $id, $qL, $qR, $sL, $sR, $sign) = (@f[0,1,2,6,7,8,9], '+');
 if ($sL > $sR) {($sL, $sR, $sign) = ($sR, $sL, '-')}
 $read =~ /_(\d+)$/; 
 $lens{$read} = $1;
 push @{$reads{$read}{$dna.$sign}}, {dna => $dna, qL =>$qL, qR =>$qR, sL =>$sL, sR =>$sR, sign => $sign, id => $id};
}

my %outs;
for my $read (sort {$lens{$b} <=> $lens{$a}} keys %lens) { # Pairwise 
 for my $set (sort keys %{$reads{$read}}) {
  my $readset = "$read\t$set"; # readID_Length\tDNA.sign
  my $ct = @{$reads{$read}{$set}};
  if ($ct == 1) {$cts{1} ++; next} # Reject if only one hit
  #$cts{$ct} ++;
  my $out = "$read\t$set\t$ct\t";
  my (%tot);
  for (sort {$$a{sL} <=> $$b{sL}} @{$reads{$read}{$set}}) {
   if ($set =~ /-$/) {for ($$_{qL}, $$_{qR}) {$_ = $lens{$read} -$_ +1} ($$_{qL}, $$_{qR}) = ($$_{qR}, $$_{qL})}
   push @{$outs{$ct}{$readset}{hits}}, {sL => $$_{sL}, sR => $$_{sR}, qL => $$_{qL}, qR => $$_{qR}, id => int($$_{id})};
   my $len= $$_{sR} - $$_{sL} + 1;
   $tot{id} += $len * $$_{id};
   $tot{len} += $len;
  }
  my $id = int($tot{id}/$tot{len});
  $outs{$ct}{$readset}{id} = $id;
  $outs{$ct}{$readset}{len} = $lens{$read};
 }
}

for my $ct (sort {$a <=> $b} keys %outs) {
 for my $readset (sort {$outs{$ct}{$b}{id} <=> $outs{$ct}{$a}{id}} keys %{$outs{$ct}}) {
  print "$readset\t$ct\t$outs{$ct}{$readset}{id}\t$outs{$ct}{$readset}{len}\n";
  @{$outs{$ct}{$readset}{hits}} = sort {$$a{sL} <=> $$b{sL}} @{$outs{$ct}{$readset}{hits}};
  for (@{$outs{$ct}{$readset}{hits}}) {print "\t$$_{sL}\t$$_{sR}\t$$_{qL}\t$$_{qR}\t$$_{id}\n"}
  next unless $ct == 2 and $outs{$ct}{$readset}{id} > 94 and $readset =~ /chromosome/;
  my $hits = $outs{$ct}{$readset}{hits};
  if ($$hits[0]{qL} < $$hits[1]{qL} and $$hits[0]{qR} < $$hits[1]{qR}) {$cts{scar}{$$hits[0]{sR} . "\t" . $$hits[1]{sL}} ++}
  if ($$hits[0]{qL} > $$hits[1]{qL} and $$hits[0]{qR} > $$hits[1]{qR}) {$cts{circle}{$$hits[0]{sL} . "\t" . $$hits[1]{sR}} ++}
 }
}
for (sort {$cts{scar}{$b} <=> $cts{scar}{$a}} keys %{$cts{scar}}) {print "$_\tscar\t$cts{scar}{$_}\n"}
for (sort {$cts{circle}{$b} <=> $cts{circle}{$a}} keys %{$cts{circle}}) {print "$_\tcircle\t$cts{circle}{$_}\n"}
