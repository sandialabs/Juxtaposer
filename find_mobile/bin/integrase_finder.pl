#!/usr/bin/perl
use strict;
use Getopt::Long;
use Pod::Usage;
#use Parallel::ForkManager;
use Data::Dumper;
use File::Basename;
use Cwd;
############################################
#                Globals                   #
############################################
our $verbose;
our $cpu = 1;
my $faa;
my $gff;
my $pfam;
my $xer;
my $intfam;
our $out;
my $scriptname = $0;
our $dir = dirname($0);
my $VERSION = '0.1';
our $directory;
our $prodigal;
our $seq2;

############################################
#                Options                   #
############################################
if (@ARGV < 1){
    print "\n Try '$scriptname --man' for full info\n\n";
    exit(0);
}
else {
    GetOptions( 'help' => sub {pod2usage(1);},
                'version' => sub {print STDOUT "\n $scriptname version $VERSION\n"; exit()},
                'man' => sub {pod2usage(-exitstatus => 0, -verbose =>2);},
                'faa=s' => \$faa,
                'cpu=i' => \$cpu,
                'fam=s' => \$pfam,
                'int=s' => \$intfam,
                'gff=s' => \$gff,
		        'xer=s' => \$xer,
		        'out' => \$out,
                'directory=s' => \$directory,
                'prodigal' => \$prodigal,
                'sequence' => \$seq2,
        );
}
print "Checking recombinases\n";
my $recombinase_list = check_integrases($faa, $pfam);
my $recombinases = read_fasta($faa, $recombinase_list);
my $integrases = remove_xers($recombinases, $xer);
my $phage_ints = remove_integrons($faa, $integrases, $intfam);
print_output($faa, $gff, $recombinases, $integrases, $phage_ints);
if ($seq2){
    print_sequence($faa);
}

sub print_sequence {
    my $faa = shift;
    my $base = $faa;
    $base =~ s/\.faa$//g;
    my $xer = $base.".xer";
    my $phg = $base.".phage";
    my $integrons = $base.".integrons";
    open(FILE, "<", $faa);
    my %seq_hash;
    my $header = "";
    while(<FILE>){
        my $line = $_;
        chomp $line;
        if ($line =~ m/^>/){
            my @arr = split(/\s+/, $line);
            $arr[0] =~ s/>//;
            $header = $arr[0];
        }
        else {
            $seq_hash{$header} .= $line;
        }
    }
    close FILE;
    my $outfile = $base.".int.sequence.faa";
    open(FILE2, ">", $outfile);
    foreach my $filename ($xer, $phg, $integrons){
        open(FILE, "<", $filename);
        my @file = <FILE>;
        close FILE;
        foreach my $line (@file){
            chomp $line;
            my @line_array = split(/\t/, $line);
            my $id = $line_array[8];
            $id =~ s/ID=\d+(.*?);.*/$line_array[0]$1/;
            print FILE2 ">".$id." ".$filename."\n";
            print FILE2 $seq_hash{$id}."\n";
        }
    }
    close FILE2;
}

sub print_output {
    my $faa = shift;
    my $base = $faa;
    $base =~ s/\.faa$//g;
    my $int = $base.".int";
    my $xer = $base.".xer";
    my $phg = $base.".phage";
    my $integrons = $base.".integrons";
    if ($directory){
        $int = $directory.".user.int";
        $xer = $directory.".user.xer";
        $phg = $directory.".user.phage";
        $integrons = $directory.".user.integrons";
    }
    #print "$faa\t";
    open(INT, ">", $int) or die "Cannot open $int $!";
    open(XER, ">", $xer) or die "Cannot open $xer $!";
    open(PHAGE, ">", $phg) or die "Cannot open $phg $!";
    open(INTEGRON, ">", $integrons) or die "Cannot open $integrons $!";
    my $gff = shift;
    my $recombinases = shift;
    my $integrases = shift;
    my $phages = shift;
    my %r_hash = %$recombinases;
    my %i_hash = %$integrases;
    my %p_hash = %$phages;
    my @lines = read_file($gff);
    my $id;
    my %xers;
    foreach (keys %r_hash) { $xers{$_} = $r_hash{$_} if not exists $i_hash{$_};}
    foreach my $line (@lines){
        chomp $line;
        my @l_array = split(/\t/, $line);
        ($id = $line) =~ s/.*ID=\d+(.*?);.*/$l_array[0]$1/g;
        if (exists $r_hash{$id}){
            if (!exists $xers{$id}){
                if (!exists $p_hash{$id}){
                    print INTEGRON $line."\n";
                }
            }
        }
        if (exists $xers{$id}){
	    if ($out){
		print ">".$base."|".$r_hash{$id}{head}."\n";
		print $r_hash{$id}{body}."\n";
	    }
            print XER $line."\n";
        }
        if (exists $i_hash{$id}){
            print INT $line."\n";
        }
        if (exists $p_hash{$id}){
            print PHAGE $line."\n";
        }
    }
    close INT;
    close XER;
    close PHAGE;
    close INTEGRON;
}

sub remove_integrons {
    my $faa = shift;
    my $integrases = shift;
    my %phages = %$integrases;
    my $intfam = shift;
    my $base = $faa;
    $base =~ s/\.faa$//g;
    my $hmm_tbl = $base.".integron.tbl";
    my $hmm_dom = $base.".integron.domtbl";
    my $hmm_log = $base.".integron.hmm.txt";
    if ($directory){
        my $hmm_tbl = $directory.".user.integron.tbl";
        my $hmm_dom = $directory.".user.integron.domtbl";
        my $hmm_log = $directory.".user.integron.hmm.txt";
    }
    my $runhmmer = "$dir/hmmsearch --acc -o $hmm_log --tblout $hmm_tbl --domtblout $hmm_dom --cpu $cpu --cut_tc $intfam $faa";
    print $runhmmer."\n";
    system($runhmmer);
    my @lines = read_file($hmm_tbl);
    my $id;
    foreach my $line (@lines){
        if ($line !~ m/^#/){
            my @line_array = split(/\s+/, $line);
            ($id = $line_array[0]) =~ s/.*id\|([0-9]*).*/$1/g;
            if (exists $phages{$id}){
               if (sprintf('%G', $line_array[4]) < sprintf('%G', '1.2e-24')){
                 delete $phages{$id};
                }
           } 
        }
    }
    return \%phages;
}

sub check_integrases {
    my $faa = shift;
    my $pfam = shift;
    my $base = $faa;
    $base =~ s/\.faa$//g;
    my $hmm_log = $base.".recomb.hmm.txt";
    my $hmm_tbl = $base.".recomb.tbl";
    my $hmm_dom = $base.".recomb.domtbl";
    my $runhmmer = "$dir/hmmsearch --acc -o $hmm_log --tblout $hmm_tbl --domtblout $hmm_dom --cpu $cpu --cut_tc $pfam $faa";
    system($runhmmer);
    my @lines = read_file($hmm_tbl);
    my %hash;
    my $id;
    foreach my $line (@lines){
        if ($line !~ m/^#/){
            my @line_array = split(/\s+/, $line);
            ($id = $line_array[0]) =~ s/.*id\|([0-9]*).*/$1/g;
            $hash{$id} = 1;
        }
    } 
    return \%hash;
}

sub remove_xers {
    my $fasta_ref = shift;
    my $xer = shift;
    my %fasta = %$fasta_ref;
    my @keys = keys %fasta;
    my %ints;
    my @xers;
    foreach my $key (@keys){
        my $run = "echo '>$fasta{$key}{head}\n$fasta{$key}{body}' | pfscan -f - $xer"; #ET edit Nov14,2024 remove path to binary
        my $exit = `$run`;
        if ($exit eq '') {
            $ints{$key} = $fasta{$key};
        }
        else {
            push @xers, $key;
        }
    }
    #if ($#xers < 1){
    #    my %ints2;
    #    foreach my $key (@keys){          
    #        my $run = "echo '>$fasta{$key}{head}\n$fasta{$key}{body}' | pfscan -C -1 -f - $xer";
    #        my $exit = `$run`;
    #        if ($exit eq '') {
    #            $ints2{$key} = $fasta{$key};
    #        }
    #    }
    #    return \%ints2;
    #} else {
        return \%ints;
    #}
}

sub read_fasta {
    my $file = shift;
    my $included_ref = shift;
    my %included = %$included_ref;
    my @lines = read_file($file);
    my %hash;
    my $prot_name;
    my $bool;
    foreach my $line (@lines){
    	chomp $line;
    	if ($line =~ m/^>/){
            $bool = 0;
    	    $prot_name = $line;
    	    $prot_name =~ s/>//;
            $prot_name =~ s/\s+.*//g;
            $prot_name =~ s/.*id\|([0-9]*).*/$1/g;
            if (exists $included{$prot_name}){
        	    $hash{$prot_name}{head} = $prot_name;
                $bool = 1;
            }
    	} else {
            if ($bool == 1){
        	    $hash{$prot_name}{body} .= $line;
            }
    	}
    }
    return(\%hash);
}

sub read_file {
    my $filename = shift;
    open(FILE, "<", $filename);
    my @file = <FILE>;
    close FILE;
    return @file;
}
