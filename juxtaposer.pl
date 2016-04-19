#! /usr/bin/perl -w
use strict;
use File::Basename;
my $scriptdir = dirname(__FILE__);
die "Usage: perl $0 refPrefix readLength\n" unless @ARGV == 2;
my ($refpfx, $readlen) = @ARGV;
die "requires multi-replicon $refpfx.fa file and read qf.fq file\n" unless -f "$refpfx.fa" and -f 'qf.fq';

# FIND MAXIMAL READLENGTH
#my $readlen = 0; # Maximum read length
#open IN, 'qf.fq'; for (1..1000) {<IN>; my $seq = <IN>; chomp $seq; my $len = length $seq; $readlen = $len if $len > $readlen; <IN>; <IN>} close IN;

# PROCESS REFERENCE GENOME
doCall("perl $scriptdir/faSizes.pl $refpfx.fa > ${refpfx}.lens") unless -f "${refpfx}.lens";
doCall("perl $scriptdir/find_mobile/mobilome_hunter.pl --fna $refpfx.fa") unless -f "$refpfx.mobile.bed";
my $clos = "${refpfx}_clos_$readlen";
unless (-f "$clos.1.bt2") {
 doCall("perl $scriptdir/closer.pl $refpfx $readlen") unless -f "$clos.fa";
 doCall("bowtie2-build $clos.fa $clos"); 
}

# PROCESS READS
doCall("bowtie2 -p 16 --un nonstandard.fq -x $clos -U qf.fq > /dev/null 2> bowtie.results") unless -f "nonstandard.fq";
doCall("perl $scriptdir/fq2fa.pl nonstandard.fq nonstandard.fa") unless -f "nonstandard.fa";
doCall("makeblastdb -in $clos.fa -out $clos -dbtype nucl") unless -f "${clos}.nin";
doCall("blastn -perc_identity 95 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore btop' -query nonstandard.fa -db $clos > nonstandard.blast") unless -f "nonstandard.blast";
doCall("perl $scriptdir/hitpairs.pl ${refpfx} $readlen > counts.txt") unless -f "juxtas.txt";

sub doCall {
 my $call = $_[0];
 my $code = system $call;
 die "Call failed, exiting: $call\nerror code $code\n" if $code;
}
