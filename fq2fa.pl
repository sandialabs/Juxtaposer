#!/usr/bin/perl
use strict;

unless ( $ARGV[0] ) { die "usage: $0 inFastq [outFasta]\nunless outFasta given, will name output file: inFastq.fa\n" }
my $infile = $ARGV[0] ;
my $outfile = "$infile.fa" ;
if ( $ARGV[1] ) { $outfile = $ARGV[1] }
open (IN, "<$infile") or die "Cannot access file $infile.\n" ;
open (OUT, ">$outfile") or die "Can't write $outfile\n" ;
#open ( Q , ">$outfile.q" ) ;
my ( $header , $seq , $q , $ct ) ;
while ($header = <IN>) {
  $seq = <IN>;
  <IN>;
  $q = <IN>;
  $ct ++ ;
  chomp $seq;
  my $len = length $seq;
  $header =~ s/^\@(\S+)/>$1_$len/ ;
  print OUT "$header$seq\n";
  #print Q "$header$q";
  #print OUT ">$ct\n$seq";  print Q ">$ct\n$q" ;
}    
close (IN);
close (OUT);
#close Q ;
