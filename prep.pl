#! /usr/bin/perl -w
use strict;

chdir 'exp5';
for (1..8) {
#for (3) {
 chdir "S$_";
 #system "perl ../../perl/juxtaposer.pl > juxlog";
 system "perl ../../attct/attct.pl qf.fq > attct.log";
 system "perl ../../attct/attct.pl qf.fq > attct.log";
 #system "perl ../../perl/readMersSC.pl qf.fq 21 k";
 #system "perl ../../perl/covmap.pl k.cov ref.fa k.map";
 #system "perl ../collate.pl k.map > collate.out";
 chdir "..";
 next;
 system "ln -s ../../../klebs/genomes/KpnBAA2146.fna ref.fa"; chdir ".."; next;
 system "ln -s ../../attct/atts.fa atts.fa";
 for (qw/ref.nhr ref.nin ref.nsq refperm.1.bt2 refperm.2.bt2 refperm.3.bt2 refperm.4.bt2 refperm.fa refperm.rev.1.bt2 refperm.rev.2.bt2/) {
  system "ln -s ../S1/$_ $_";
 }
 chdir "..";
}

chdir '../exp4';
for (qw/log0mm0h log1mm2h log5mm2h sty1mm2h log1mm1h log5mm1h sty0mm0h sty5mm2h/) {
 chdir $_;
 system "ln -s ../../attct/atts.fa atts.fa";
 #system "perl ../../attct/attct.pl qf.fq > attct.log";
}
