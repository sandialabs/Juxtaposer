my @cats = qw/Tn6187 ISKpn1 ISKPn14 ISKpn18 ISKpn21 IS1R IS1F IS1X4 IS26 IS3000 IS6100 IS903B ISEcp1 ISEc22 ISEcl1 IS5075circle1 IS4321circle1/;
my ($ct, %types, $cat, $entry);
while (<>) {
 s/.*-//;
 chomp;
 $types{$cat}{$_} = $entry unless />/;
 s/qf.fq.fa.hits://;
 next unless /(.*)>(\d+\_\d+\_\d+\_(\S+))/;
 ($cat, $entry) = ($3, $1 . $2);
}

for my $cat (@cats) {
 next unless $types{$cat};
 open OUT, ">tposonFams/$cat";
 for my $seq (sort keys %{$types{$cat}}) {
  print OUT ">$types{$cat}{$seq}\n$seq\n";
 }
 close OUT;
 `muscle -in tposonFams/$cat -out tposonFams/$cat.align > /dev/null`;
}
  
