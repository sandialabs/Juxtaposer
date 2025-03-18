#!/usr/bin/perl
use strict;
use Getopt::Long;
use Pod::Usage;
use Attribute::Handlers;
use Cwd;
use File::Basename;
############################################
#                Globals                   #
############################################
our $verbose;
our $scriptname = $0;
our $dir = dirname($0);
our $VERSION = '0.1';
my ($cpu, $fna) = (1, "");
############################################
#                Options                   #
############################################
if (@ARGV < 1){
    print "\n Try '$scriptname --man' for full info\n\n";
    exit(0);
}
else{
    GetOptions('help' => sub {pod2usage(1);},
	       'version' => sub {print STDOUT "\n $scriptname version $VERSION\n"; exit()},
	       'man' => sub {pod2usage(-exitstatus => 0, -verbose => 2);},
	       'verbose' => \$verbose,
	       'fna=s' => \$fna,
               'cpu=i' => \$cpu,
	) ;
}
unless(-e $fna){
	print "Need fna file to run mobilome_hunter.pl\n";
	exit(0);
}

my ($faa, $gff) = run_prodigal($fna);
my ($phage_integrase, $integron_integrase) = run_integrase_finder($faa, $gff);
my ($mobile_tbl) = run_tpn_finder($faa);
collect_mobilome($phage_integrase, $integron_integrase, $mobile_tbl, $gff);

sub collect_mobilome{
	my $phage_integrase = shift;
	my $integron_integrase = shift;
	my $mobile_tbl = shift;
	my $gff = shift;
	my %mobile_hash;
	my $output = "";
	my ($filename, $dirs, $suffix) = fileparse($gff, qr/\.[^.]*/);
	my $output_gff = $dirs.$filename.".mobile.gff";
	my $output_bed = $dirs.$filename.".mobile.bed";
	open(FILE, "<", $phage_integrase);
	while(<FILE>){
		my $line = $_;
		chomp $line;
		my $value =  $line."MobileType=PhageIntegrase;\n";
		$output .= $value;
	}
	close FILE;
	open(FILE, "<", $integron_integrase);
	while(<FILE>){
		my $line = $_;
		chomp $line;
               	my $value =  $line."MobileType=IntegronIntegrase;\n";
		$output .= $value;
	}
	close FILE;
	my $stdout = `sort $mobile_tbl -k1,1 -k3g | awk '(\$5 < 1e-04){print}' | awk '((\$1!=h)&&(\$1 !~ /^#/)) {print} {h=\$1}'`;
	my @tpn_out = split(/\n/, $stdout);
	my %tpn_hash;
	my $id;
	foreach (@tpn_out){
		my $line = $_;
		chomp $line;
		my @line_array = split(/\s+/, $line);
		$id = $line_array[-1];
		$id =~ s/.*ID=(.*?);.*/$1/;
		$tpn_hash{$id} = "MobileType=$line_array[2];Score=$line_array[4];"
	}
	open(FILE, "<", $gff);
	while(<FILE>){
		my $line = $_;
		chomp $line;
	        my @l_array = split(/\t/, $line);
        	($id = $line) =~ s/.*ID=(.*?);.*/$1/g;
        	if (exists $tpn_hash{$id}){
	               	my $value =  $line."$tpn_hash{$id};\n";
			$output .= $value;
		}
	}
	close FILE;
	my @outputs = split(/\n/, $output);
	my @sorted =
		map {$_->[0]}
		sort {$a->[1] cmp $b->[1]
			or $a->[2] <=> $b->[2]
			or $a->[0] cmp $b->[0]
		}
		map { /(.*?)\t(.*?)\t(.*?)\t(.*?)\t.*/; [$_, $1, $4]}
		@outputs;
	open(GFF, ">", $output_gff);
	open(BED, ">", $output_bed);
	my %ids;
	foreach my $sorted_line (@sorted){
		my @bed_array = split(/\t/, $sorted_line);
		$sorted_line =~ s/(.*)\t.*/\1/g;
		print BED $bed_array[0]."\t".$bed_array[3]."\t".$bed_array[4]."\t";
		$bed_array[-1] =~ s/.*MobileType=(.*?);.*/\1/;
		$bed_array[-1] =~ s/curated//g;
		$bed_array[-1] =~ s/combined//g;
		$bed_array[-1] =~ s/original//g;
		$bed_array[-1] =~ s/[_|\.]*$//g;
		if (exists $ids{$bed_array[-1]}){
			$ids{$bed_array[-1]}++;
		} else {
			$ids{$bed_array[-1]} = 1;
		}
		print BED $bed_array[-1].".$ids{$bed_array[-1]}\t.\t".$bed_array[6]."\n";
		print GFF $sorted_line."\tID=$bed_array[-1].$ids{$bed_array[-1]};\n";
	}
	close BED;
	close GFF;
}

sub run_tpn_finder{
	my $faa = shift;
	my ($filename, $dirs, $suffix) = fileparse($faa, qr/\.[^.]*/);
	my $tbl = $dirs.$filename.".mobile.tbl";
	my $domtbl = $dirs.$filename.".mobile.domtbl";
	`$dir/bin/hmmsearch --tblout $tbl --domtblout $domtbl --noali --cut_tc --cpu $cpu $dir/db/TnpPred_HMM_Profiles3.hmm $faa`;
	return($tbl);
}

sub run_integrase_finder{
	my $faa = shift;
	my $gff = shift;
        my ($filename, $dirs, $suffix) = fileparse($fna, qr/\.[^.]*/);
	my $phage_integrase = $dirs.$filename.".phage";
	my $integron_integrase = $dirs.$filename.".integrons";
	my $output = `perl $dir/bin/integrase_finder.pl --faa $faa --gff $gff --cpu $cpu --fam $dir/db/PF00589_seed.hmm --int $dir/db/famint9.hmm --xer $dir/db/xers.prf --out`;
	print $output;
	return($phage_integrase, $integron_integrase);
}

sub run_prodigal{
	my $fna = shift;
	my ($filename, $dirs, $suffix) = fileparse($fna, qr/\.[^.]*/);
	my $faa = $dirs.$filename.".faa";
	my $gff = $dirs.$filename.".gff";
	`$dir/bin/prodigal -a $faa -f gff -i $fna -o $gff`;
	return($faa, $gff);
}

## POD Documentation
=pod

=head1 NAME

mobileome_hunter.pl

=head1 VERSION

Documentation for mobilome_hunter.pl version 0.1

=head1 SYNOPSIS

mobilome_hunter.pl PATH_TO_FNA

=head1 DESCRIPTION

This program reads in an .fna file and produces a list of mobile elements.

=head1 USAGE

    Examples:
    
    perl mobilome_hunter.pl PATH_TO_FNA

=head1 AUTHOR
   
Written by Corey M. Hudson and Kelly P. Williams.

=head1 DEPENDENCIES

uses Getopt::Long, Pod::Usage, List::Util qw(max min), Attribute::Handlers
These should be standard in most implementation of perl. More important than this
are the programs necessary to determine tRNAs and integrases.

=head1 LICENSE AND COPYRIGHT

Copyright (2016) Sandia Corporation. Under the terms of
Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
retains certain rights in this software.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

=cut
