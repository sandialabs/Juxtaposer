#! /usr/bin/perl -w
use strict;

my %islein = qw/TandL 583102-615770 49R 1296810-1345944 38RybB 1785570-1823069 29S 1966140-1994683 23SapB 2286456-2309756 42GloC 2345718-2388107 40GuaA 3969748-4010194 37X 4129454-4166020 55F 4603724-4658623 16Fis 4919120-4935114/;
my (%covs, %isles, %distrs);
for (keys %islein) { 
 #print "$_ $islein{$_}\n"; next;
 $islein{$_} =~ /^(\d+)-(\d+)$/;
 #die "$islein{$_} $_\n";
 my $len = $2-$1+1;
 %{$isles{$_}} = (L => $1, R => $2, len => $len);
 my $temp = $len;
 $len = 16000 if $_ eq '42GloC';
 %{$isles{$_."fl"}} = (L => $1-$len, R => $1-1, len => $len);
 $temp = $len;
 $len = 16000 if $_ eq '23SapB';
 %{$isles{$_."fr"}} = (L => $2+1, R => $2+$len, len => $len);
}
#for (sort {$isles{$a}{len} <=> $isles{$b}{len}} keys %isles) {print "$_ $isles{$_}{L} $isles{$_}{R} $isles{$_}{len}\n"} exit;

open IN, $ARGV[0] or die "No infile $ARGV[0]\n";
while (<IN>) {
 #pCuAs	1	115	1
 chomp;
 my @f = split /\s+/;
 next unless $f[3] == 1;
 my $flag = '';
 if ($f[0] eq 'chromosome') {
  for (keys %isles) {
   next unless $f[1] >= $isles{$_}{L};
   $flag = $_ if $f[1] <= $isles{$_}{R};
  }  
 }
 #if (/^chromosome\s+2309757/ or /^chromosome\s+2286455/) {print "$_\n'$f[0]' '$f[1]' '$f[2]' '$f[3]' '$flag'\n"}
 if ($flag) {push @{$covs{$flag}}, $f[2]} #; print "$_\n" if $flag =~ /f[lr]23SapB/}
 unless ($flag and $flag !~ /^f/) {push @{$covs{$f[0]}}, $f[2]}
}
close IN;

my @order; 
for (qw/16Fis 23SapB 29S TandL 37X 38RybB 40GuaA 42GloC 49R 55F/) {push @order, $_."fl", $_, $_."fr"} 
for (qw/pMYS-2014 pHg-85164 pCuAs-117755 pNDMUS-140825 chromosome-5435369/) {/(\S+)-(\S+)/; push @order, $1; $isles{$1}{len} = $2}
for my $isle (@order) {
 print "$isle\t$isles{$isle}{len}"; for (NMedMad(@{$covs{$isle}})) {print "\t$_"} print "\n";
}

sub NMedMad {
 my $n = @_;
 if ($n < 1) {return 0, 0, 0, 0, 0, 0, 0}
 my ($mid, $mean, $med, $sumsquares, $mad, $outhi, $outlo, @devs) = (int($n/2)-1, 0, 0, 0, 0, 0, 0);
 my @values = sort {$a <=> $b} @_;
 if ($n % 2) {$med = $values[$mid]} else {$med = ($values[$mid] + $values[$mid+1])/2}
 for (@values) {push @devs, abs($_ - $med); $mean += $_}
 @devs = sort {$a <=> $b} @devs;
 if ($n	% 2) {$mad = 1.4826 * $devs[$mid]} else {$mad = 1.4826 * ($devs[$mid] + $devs[$mid+1])/2} 
 my ($hi, $lo) = ($med+2.5*$mad, $med-2.5*$mad);
 $mean /= $n;
 for (@values) {
  $outhi ++ if $_ > $hi;
  $outlo ++ if $_ < $lo;
  $sumsquares += ($_ - $mean) ** 2;
 }
 return $n, $mean, sqrt($sumsquares/($n - 1)), $med, $mad, $outhi, $outlo;
}
