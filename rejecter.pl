my %reads;
while (<>) {
 my @f = split "\t";
 my $sense = abs($f[9]-$f[8])/($f[9]-$f[8]);
 #die "$f[9]-$f[8] $sense\n";
 $reads{$f[0]}{$f[1].$sense} ++;
}

my %ct;
for my $read (keys %reads) {
 my ($flag, $sum);
 for (keys %{$reads{$read}}) {
  $sum += $reads{$read}{$_};
  if ($reads{$read}{$_} > 1) {$flag ++}
 }
 $ct{$sum} ++;
 $ct{flag} ++ if $flag;
}

for (keys %ct) {
 print "$_\t$ct{$_}\n";
}

