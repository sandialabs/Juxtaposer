#!/usr/bin/perl
use strict;
use File::Slurp;
use Parallel::ForkManager;
my $fasta_ref = read_fasta($ARGV[0]);
my %fasta = %$fasta_ref;
my @keys = keys %fasta;
my $pm = Parallel::ForkManager->new($ARGV[1]);
foreach my $key (@keys){
    my $pid = $pm->start and next;
    my $run = "echo '>$fasta{$key}{head}\n$fasta{$key}{body}' | ./pfscan -f - xer.prf";
    my $exit = `$run`;
    if ($exit ne '' && $key ne ''){
	print "$key\t$exit";
    }
    $pm->finish;
}

sub read_fasta {
    my $file = shift;
    my @lines = read_file($file);
    my %hash;
    my $prot_name;
    foreach my $line (@lines){
	chomp $line;
	if ($line =~ m/^>/){
	    $prot_name = $line;
	    $prot_name =~ s/>//;
	    $hash{$prot_name}{name} = $prot_name;
	} else {
	    $hash{$prot_name}{body} .= $line;
	}
    }
    return(\%hash);
}
