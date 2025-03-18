#! /usr/bin/perl
use strict; use warnings;
use File::Basename;
my $scriptdir = dirname(__FILE__);
die "Usage: perl $0 refPrefix [cpu]\n" unless @ARGV > 0;
my $palflag = 0;
my ($refpfx, $cpu) = @ARGV;
$cpu = 1 unless $cpu;
die "cpu must be an integer\n" if $cpu =~ /[^\d]/;
my ($readfile, $zip) = ('qf.fq', 0);
if (-f 'qf.fq') {} elsif (-f 'qf.fq.gz') {($readfile, $zip) = ('qf.fq.gz', 1)} else {die "requires read qf.fq or qf.fq.gz file\n"}
die "requires multi-replicon $refpfx.fa file\n" unless -f "$refpfx.fa";

# FIND MAXIMAL READLENGTH
my $readlen = 0; # Maximum read length
if ($zip) {open(IN, sprintf("zcat %s |", $readfile))} else {open IN, $readfile} 
for (1..1000) {<IN>; my $seq = <IN>; last unless $seq; chomp $seq; my $len = length $seq; $readlen = $len if $len > $readlen; <IN>; <IN>}
close IN;
warn "Largest read of first 1000 is $readlen\n";

# PROCESS REFERENCE GENOME
doCall("perl $scriptdir/faSizes.pl $refpfx.fa > ${refpfx}.lens") unless -f "${refpfx}.lens";
my $mobile_cpu = $cpu -1;  # Hmmsearch with cpu==0 is safest
doCall("perl $scriptdir/find_mobile/mobilome_hunter.pl --fna $refpfx.fa --cpu $mobile_cpu") unless -f "$refpfx.mobile.bed";
my $clos = "${refpfx}_clos_$readlen";
unless (-f "$clos.1.bt2") {
 doCall("perl $scriptdir/closer.pl $refpfx $readlen") unless -f "$clos.fa";
 doCall("bowtie2-build $clos.fa $clos"); 
}

# PROCESS READS
doCall("bowtie2 -p $cpu --un nonstandard.fq -x $clos -U qf.fq* > /dev/null 2> bowtie.results") unless -f "nonstandard.fq";
doCall("perl $scriptdir/fq2fa.pl nonstandard.fq nonstandard.fa") unless -f "nonstandard.fa";
doCall("makeblastdb -in $clos.fa -out $clos -dbtype nucl") unless -f "${clos}.nin";
doCall("blastn -perc_identity 95 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore btop' -query nonstandard.fa -db $clos > nonstandard.blast") unless -f "nonstandard.blast";
doCall("perl $scriptdir/hitpairs.pl ${refpfx} $readlen $palflag > counts.txt") unless -f "juxtas.txt";
warn "//FINISHED\n";

sub doCall {
 my $call = $_[0];
 warn "$call\n";  # Log the call
 my $code = system $call;
 die "Call failed, exiting: $call\nerror code $code\n" if $code;
}
