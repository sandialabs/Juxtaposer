#!/usr/bin/perl
use strict;
#use warnings;
use List::Util qw(sum max min);
use File::Basename;
use Getopt::Long;
use POSIX qw(ceil) ;
#use Bio::Tools::SeqPattern ;

# VARIABLES , DEFAULTS
my $qscore_threshold = 30 ;
my $endscore_threshold = 25 ;
my $qtrim_method = 'bwa' ;
my $nametype = 'serial' ;
my $qual_type = 3 ;
my $guess_qual = 0 ;
my $debarcode = 0 ;
my $paired = '' ;
my $primersets_available ;
my $partmin = 5 ;
my $minLen = 30 ;
my $primerlist ;
my $NLen = 1 ;
my $Nright = 0;
my $dustoff ;
my $debug ;
my $version = "1.9";
my $threads = 1 ;

#other global variables:
my ( $infile , $prefix , $verbose , $outnames , $outfastq , %ct ) ;
my $t0 = time ;
my @time = localtime(time) ; #sec,min,hour,mday,mon,year,wday,yday,isdst
my $timestamp = sprintf "%4d-%02d-%02d_%02d:%02d:%02d" , $time[5] + 1900 , $time[4] + 1 , @time[3,2,1,0] ;
my $call = "on $timestamp, perl $0 @ARGV" ;
my $RAPTOR_PATH = dirname(__FILE__);
$RAPTOR_PATH =~ s/scripts$//; 
my $primerfile = "$RAPTOR_PATH/library/RaptorPrimers.txt" ;
my %dusthash ;

my %ambig = (M => "[AC]" ,
             R => "[AG]" ,
             W => "[AT]" ,
             S => "[CG]" ,
             Y => "[CT]" ,
             K => "[GT]" ,
             V => "[ACG]" ,
             H => "[ACT]" ,
             D => "[AGT]" ,
             B => "[CGT]" ,
             N => "[GATC]" ) ;

# MESSAGES , ARGUMENTS
if ( @ARGV < 1 ) {
 print STDERR <<EOF;
FastQ pre-alignment (quality) filter. Rejects any-N, too short after min-Q tail trim. Low-complexity filtering. Scores and sorts on avg-Q. 
Removes primers in forward or reverse orientation using 14-mer parts with up to 1 base-substitution, and nibbles smaller 3' primer fragments at termini.
Options:
-qscore-threshold <int> avg quality score threshold for filtering. Default $qscore_threshold
-endscore-threshold <int> trim in from end all bases at or below this quality.  Default $endscore_threshold
-qtrim-method <simple|bwa> algorithm for deciding where to trim low quality bases. Default $qtrim_method
-nametype <serial|original> sequence names: serial numbers or original. Default $nametype
-guess_qual Guess between phred+33 and phred+64 for -infile. Reads 100 qual lines, guesses, reports with analysis, then exits.
-qual-type <int> 0=sanger qualities (phred+33), 1=illumina qualities pipeline>=1.3(phred+64), 2=illumina qualities pipeline<1.3. 3=guess between 0 and 1. Default $qual_type.
-minLen <int> Minimum nt length post-trim. Default $minLen.
-debarcode <int> Number of bases to trim in-sequence barcode from 5' end. Default no-debarcode.
-primers <pre-named set or custom comma-separated string of primer sequences> For list of pre-named sets: $0 -primersets_available
-primersets_available Lists primer sets and sequences available, then exits.
-partmin <integer> For terminal primer-matching, minimal length of 3' end. Default = $partmin.
-Nlength <num> number of consecutive ambiguous bases (N) causing rejection. Default = $NLen.
-Nright trims everything righward from first N. Overrides -Nlength. Default = off.
-dustoff turns off DUST filter.
-debug 
-infile <input file> fastq file. Required. 
-prefix <prefix used for all output files in a single run> Required.
-threads <integer> number of threads to use for dust-masker.
EOF
 exit(1);
}

GetOptions ( 
 'qscore-threshold=i' => \$qscore_threshold,
 'endscore-threshold=i' => \$endscore_threshold,
 'qtrim-method=s' => \$qtrim_method,
 'nametype=s' => \$nametype,
 'qual-type=i' => \$qual_type,
 'guess_qual' => \$guess_qual ,
 'debarcode=i' => \$debarcode ,
 'primers=s' => \$primerlist ,
 'primersets_available' => \$primersets_available ,
 'partmin=i' => \$partmin ,
 'minLen=i' => \$minLen ,
 'Nlength=i' => \$NLen ,
 'Nright' => \$Nright ,
 'dustoff' => \$dustoff ,
 'debug' => \$debug ,
 'infile=s' => \$infile,
 'prefix=s' => \$prefix,
 'verbose' => \$verbose,
 'threads=i' => \$threads
) ;

if ( !defined $guess_qual and !defined $primersets_available and ( !defined $infile || !defined $prefix ) ) {
 print STDERR "Error: Must specify infile and prefix\n";
 exit(1);
}
if ( $infile eq "$prefix.fq" ) { die "Prefix choice '$prefix' will cause overwrite of infile '$infile'. Choose another prefix.\n" }
if ($Nright) {$NLen = 0}
my ( $header , $seq , $qual , %scorecounts , %seqs , %seqload , %outtypes , %dupes , $seqsperfile , $guess, $score ) ;
my $ct = 0 ;
my $serial = 0 ;
my $infilecat ;

my $qtrim ;
if    ( $qtrim_method eq 'bwa'    ) { $qtrim = \&qtrim_bwa ; }
elsif ( $qtrim_method eq 'simple' ) { $qtrim = \&qtrim_simple ; }
else  { die "Unknown quality trim method.\n"; }

# PROCESS PRIMER FILE
my ( $type , $collect , %primers , @primerparts ) ;
if ( open IN , $primerfile ) {
    while ( <IN> ) {
	chomp ;
	if ( /^(PRIMERS(ETS)*)/ ) { $type = $1 }
	elsif ( /^\t*$/ ) { next }
	elsif ( /^UniqueName/ ) { next }
	elsif ( /EXPERIMENTS/ ) { last }
	elsif ( /^#(\S+)/ ) { 
	    $collect = $1 ; 
	    if ( $primersets_available ) { print "\nPRIMER SET $1\n" }
	} 
	elsif ( $type eq 'PRIMERS' ) { /^(\S+)\s+(\S+)/ ; $primers{seqs}{$1} = $2 }
	elsif ( $type eq 'PRIMERSETS' ) { 
	    /^(\S+)/ ; 
	    my $primer = $1 ;
	    if ( $primersets_available ) { print "$primer" }
	    if ( $primer =~ s/#$/1/ and $primersets_available ) { print ", e.g. (uppercase barcode varies)" }
	    if ( $primersets_available ) { print "\t$primers{seqs}{$primer}\n" }
	    push @{$primers{sets}{$collect}} , $primer ;
	} 
    }
    close IN ;
}
if ( $primersets_available ) { 
    unless ( scalar keys %primers and scalar keys %{$primers{sets}} ) { 
	print "None available: " ;
	if ( -f $primerfile ) { print "check parseability of RaptorPrimers.txt file\n" }
	else { print "$primerfile file is missing\n" }
    } 
    print "\n" ;
    exit ( 1 ) 
}


if ( -f $infile) {
    if ($infile =~ m/\.gz$/) { $infilecat = "zcat $infile" }
    else                     { $infilecat = "cat $infile" }
}
else { die "Infile $infile not found." }

# GUESS QUALITY TYPE
if ( ( $guess_qual or $qual_type == 3 ) and not $primersets_available ) {
    $qual_type = guess_qual_type($infilecat);
}


# PROCESS PRIMERS (MAKE PRIMER PARTS)
if ($primerlist eq 'empty') {$primerlist = ''}
my ( $partmax ) = ( 14 ) ;
if ( $verbose ) { print time-$t0, " ---> reading primers and making primerparts with max length=$partmax\n" }
if ( $primerlist ) {
    my ( @primers , %uniqs ) ;
    if ( $primers{sets} and $primers{sets}{$primerlist} ) { 
	for ( @{$primers{sets}{$primerlist}} ) { 
            push @primers , $primers{seqs}{$_};
            if ( /(L|R)(EFT|IGHT)ONLY/ ) { $primers[-1] .= ">$1" }
        }
    } else {
	if ( $primerlist !~ /,/ ) {
	    my ( $testdna , $len ) = ( $primerlist , length $primerlist ) ;
	    $testdna =~ s/[acgt]//gi ;
	    if ( ( length $testdna ) / $len > 0.2 ) { die "$primerlist doesn't look like DNA, and is not an available primer set\n" }
	}
	for ( split "," , $primerlist ) { push @primers , $_ }
    }
    for my $primer ( @primers ) {
	$primer = uc $primer ;
	print "primer $primer\n" if $debug ;
	my ( $seq , $m , $side ) = ( '' , 0 , '' ) ;
	if ( $primer =~ s/>(.)$// ) { $side = $1 } # Special symbol for side pref, added above not in primer file
	while ( $primer =~ s/[^ACGT]$// ) { $m ++ } #trim 3' ambiguous characters
	for ( qw/A C G T/ ) { if ( $primer =~ s/($_{5})($_+)$/$1/ ) { $m += length $2 } } #trim 3' homopolymers until 5 nt
	my $part = substr $primer , -1 * $partmax , $partmax ;

	if ( not defined $uniqs{$part} and not $side eq 'R' ) { # 3' 14-mer
	    $uniqs{$part} = scalar @primerparts ; $ct{part} ++ ;
	    push @primerparts , {part => $part, 'm' => $m, len => length($part), serial => $ct{part}, dir => 'f', threeprime => 1};
	} elsif ( not $side eq 'R' and $m > $primerparts[$uniqs{$part}]{'m'} ) { $primerparts[$uniqs{$part}]{'m'} = $m }
	$part = Revcomp ( $part ) ;
	if ( not defined $uniqs{$part} and not $side eq 'L' ) { # 3' 14-mer backwards
	    $uniqs{$part} = scalar @primerparts ; $ct{part} ++ ;
	    push @primerparts , {part => $part, 'm' => $m, len => length($part), serial => $ct{part}, dir => 'r', threeprime => 1};
	} elsif ( not $side eq 'L' and $m > $primerparts[$uniqs{$part}]{'m'} ) { $primerparts[$uniqs{$part}]{'m'} = $m }
	$part = substr $primer , 0 , $partmax , '' ;
	$m += length $primer ;
	if ( not defined $uniqs{$part} and not $side eq 'R' ) { # 5' 14-mer
	    $uniqs{$part} = scalar @primerparts ; $ct{part} ++ ;
	    push @primerparts , {part => $part, 'm' => $m, len => length($part), serial => $ct{part}, dir => 'f', threeprime => 0};
	} elsif ( not $side eq 'R' and $m > $primerparts[$uniqs{$part}]{'m'} ) { $primerparts[$uniqs{$part}]{'m'} = $m }
	$part = Revcomp ( $part ) ;
	if ( not defined $uniqs{$part} and not $side eq 'L' ) { # 5' 14-mer backwards
	    $uniqs{$part} = scalar @primerparts ; $ct{part} ++ ;
	    push @primerparts , {part => $part, 'm' => $m, len => length($part), serial => $ct{part}, dir => 'r', threeprime => 0};
	} elsif ( not $side eq 'L' and $m > $primerparts[$uniqs{$part}]{'m'} ) { $primerparts[$uniqs{$part}]{'m'} = $m }
    }
}
if ( $verbose or $debug ) { for ( @primerparts ) { print "Primerpart $$_{serial}: $$_{part} m=$$_{m} dir=$$_{dir} 3'=$$_{threeprime}\n" } }

# preset conversion subroutine to avoid per-position choice
my ( $sum_qual_score, $mean_qual_score, $offset ) ;
if ( $qual_type == 0 ) { $sum_qual_score = \&sum_qual_score_0 ; $mean_qual_score = \&mean_qual_score_0 ; $offset = 33 ; }
if ( $qual_type == 1 ) { $sum_qual_score = \&sum_qual_score_1 ; $mean_qual_score = \&mean_qual_score_1 ; $offset = 64 ; }
if ( $qual_type == 2 ) { $sum_qual_score = \&sum_qual_score_2 ; $mean_qual_score = \&mean_qual_score_2 ; $offset = 64 ; }
my $min_symb = join ('', ("[", map ( {chr($_ + $offset)} (0..$endscore_threshold) ), "]" ) ); 

my $qualregex = qr/$min_symb+$/ ;

## PRECOMPILE REGULAR EXPRESSIONS
if ( $verbose ) { print time-$t0, " ---> pre-compiling regular expressions\n" ; }
for my $part ( sort { $$b{'m'} <=> $$a{'m'} } @primerparts ) { 
    my @re ;  # for collecting pre-compiled maxpart regexes
    if ( $debug ) { print "Primerpart $$part{serial}: $$part{part} m=$$part{'m'} dir=$$part{dir}\n" }
    for my $j ( 0 .. $$part{len} - 1 ) { #make each single-substitution variant regex
	my $var = $$part{part} ; 
	substr $var , $j , 1 , '.' ;
	if ( $j == 0 or $j == $$part{len} - 1 ) { $var =~ s/\.// } #remove terminal wildcard 
	if ( $$part{dir} eq 'f' ) { # remove primer and left flank
	    if ( $$part{'m'} and $j != $$part{len} - 1 ) { $var .= ".{0,$$part{m}}" }
	    if ( $j == $$part{len} - 1 ) { $var .= ".{0," . ( $$part{'m'} + 1 ) . "}" } #last variant removes 3' nt, raising m
	    foreach my $amb (keys %ambig) { my $repl = $ambig{$amb}; $var =~ s/$amb/$repl/g ; }  # to deal with ambiguity codes
	    push @re , qr/.*$var/ ; 
	} else { # reverse part, to remove reverse primer and right flank
	    my $pre = '' ;
	    if ( $$part{'m'} and $j != 0 ) { $pre = ".{0,$$part{m}}" }
	    if ( $j == 0 ) { $pre = ".{0," . ( $$part{'m'} + 1 ) . "}" } #first variant removes 3' nt, raising m
	    foreach my $amb (keys %ambig) { my $repl = $ambig{$amb}; $var =~ s/$amb/$repl/g ; }  # to deal with ambiguity codes
	    push @re , qr/$pre$var.*/ ; 
	}
    }
    $$part{regex1} = [@re]  ;
    if ($debug) { foreach ( @{$$part{regex1}} ) { print "$_\n"} }
    my %re2 ; # for collecting pre-compiled short-mer regexes
    if ( $$part{threeprime} ) {    # skip non-3' primer parts
	for my $n ( 2 .. $partmax - $partmin ) {     # start with $n=2 because 13 mers have already been checked
	    my ( $var , $plus ) = ( $$part{part} , '' ) ; 
	    if ( $$part{dir} eq 'f' ) { $var =~ s/^.{$n}// } else { $var =~ s/.{$n}$// } 
	    if ( $$part{'m'} ) { $plus = ".{0,$$part{m}}" }
	    foreach my $amb (keys %ambig) { my $repl = $ambig{$amb}; $var =~ s/$amb/$repl/g ; }  # to deal with ambiguity codes
	    if ( $$part{dir} eq 'f' ) { 
		$re2{$n} = qr/^$var$plus/  ; 
	    } else { 
		$re2{$n} = qr/$plus$var$/  ; 
	    }
	}
	$$part{regex2} = {%re2} ;
       	if ($debug) { foreach ( sort {$a <=> $b} keys %{$$part{regex2}} ) {print "$$part{regex2}{$_}\tn=$_\n" } }
    }
}

# MAIN LOOP
unless ($dustoff) {
    my $dustheader ;
    my $dustcmd ;
    
    if ($threads > 1) { $dustcmd = "$infilecat | sed -n '1~4s/^@\\(\\S\\+\\).*/>\\1/p;2~4p' | parallel -j $threads --block 10M --recstart '>' --pipe '$RAPTOR_PATH/bin/dustmasker | grep -B 1 ^[0-9]' > $prefix.dust"}
    else { $dustcmd = "$infilecat | sed -n '1~4s/^@\\(\\S\\+\\).*/>\\1/p;2~4p' | $RAPTOR_PATH/bin/dustmasker | grep -B 1 ^[0-9] > $prefix.dust"}

    if ($verbose) { print time-$t0, " ---> running dustmasker to file $infile using $threads threads\n" }
    if ($debug) { print "Dust command: $dustcmd\n" }
    system($dustcmd);
    open DUSTIN , "< $prefix.dust" ;
    while ( <DUSTIN> ) {
	if    ( $_ =~ m/^>(\S+).*/ )          { $dustheader = $1; @{$dusthash{$dustheader}} = () }
	elsif ( $_ =~ m/([0-9]+) - ([0-9]+)/ ) { push @{$dusthash{$dustheader}} , $1, $2 }
	else                                   { next }
    }
    close DUSTIN ;
}

if ($verbose) { print time-$t0, " ---> applying quality filter to file $infile using $threads threads\n" }


open ( IN , "-|", $infilecat ) or die "Cannot access file $infile\n";
open ( FASTQ,    ">$prefix.fq" ) ;
open ( INDEX,    ">$prefix.index" ) ;
print INDEX "# Call: $call\n# Quality Type: $guess\n# Fields: readID, header, prelength, length, prescore, score, ltrim, rtrim, proc, dupes\n" ;

$serial = 0 ;
while ( $header = <IN> ) {

    my @procs;
    my $ltrim = 0;   ## how much is being removed from left hand side (cumulative)
    my $rtrim = 0;   ## how much is being removed from right hand side (cumulative)
    my $left ;       ## temporary holder to make proc reporting easier
    my $right ;      ## temporary holder to make proc reporting easier
    my $dustsum = 0;

    $header =~ s/^@(\S+).*$/$1/ ;
    $seq = uc <IN>;
    <IN>; # ignore header2
    $qual = <IN>;
    chomp ( $header, $seq , $qual ) ;

    $serial ++ ;

    my $seqlength = length $seq ;   ## original seq length -- does not change
    my $currseqlength ;
    my $origscore = &$mean_qual_score(($qual, $seqlength)) ;
    my $score = 0;

    #process B: debarcode
    if ( $debarcode ) {
	$seq  =~ s/^.{$debarcode}// ;                   ######## sequence changing length here ########
	$qual =~ s/^.{$debarcode}// ;
	push @procs, "B[1-$debarcode]" ;
	$ltrim = $debarcode ;
    }

    #process T3 and T5: low-quality tails
######################
    # $left is the left-most base of the right trim
    $left = &$qtrim ( $qual )  ;
#    print "$serial RIGHT #############$left\n$qual\n";
    if ($left > 0) {
	$rtrim += length($qual) - $left + 1;
	$seq = substr $seq , 0 , $left - 1;         ######## sequence changing length here ########
	$qual = substr $qual , 0 , $left - 1 ;
#	print "$qual\n";
 	$left += $ltrim ;   # for correct reporting
	push @procs, "T3[$left-$seqlength]" ;
    }

    # $right is the right-most base of the left trim
    $right = length($qual) - &$qtrim ( scalar (reverse ($qual))) + 1;
    if ($right > length($qual)) {$right = 0}
#    print "$serial LEFT #############$right\n$qual\n";
    if ($right > 0) {
	$seq = substr $seq, $right, length $seq ;        ######## sequence changing length here ########
	$qual = substr $qual, $right, length $qual ;
	$left = $ltrim + 1;  # for correct reporting
	$ltrim += $right;   #adjust ltrim
#	print "$qual\n";
	push @procs, "T5[$left-$ltrim]" ;
    }

    $currseqlength = length $seq ;
    if ( $currseqlength < $minLen ) {
        push @procs, 'R' ;
    }

    #process N (unless Nright on): reject sequences with >= NLen number of N's
    #Nrightward if on
    if ($procs[-1] ne 'R' and $Nright and $seq =~ s/(N.*)//) {
     $rtrim += length $1;
     $qual = substr $qual, 0, length($seq);  
     push @procs, "N";
      $currseqlength = length $seq ;
      if ( $currseqlength < $minLen ) {
          push @procs, 'R' ;
      }
    }

    if ($procs[-1] ne 'R' and $NLen and $seq =~ /N{$NLen}/)  { 
	push @procs , 'N' ;
	push @procs , 'R' ;
    }

    #process P: trim primer sequences
    if ( $procs[-1] ne 'R' ) { 
	my $lprimertrim = 0;
	my $rprimertrim = 0;
	# phase I: max-mer parts and their single-substitution variants
	PART: foreach my $part ( @primerparts ) { 
	    foreach my $re ( @{$$part{regex1}} ) {
		if ($seq =~ $re ) {
		    $left = $ltrim + $-[0] + 1;
		    $right = $ltrim + $+[0] ;
		    push @procs, "P$$part{serial}\[$left-$right]";
		    if ( $debug ) {print "$serial\t$seq\tP$$part{serial}\[$left-$right]\t$re\t$-[0]\t$+[0]\n" } 
		    if    ( $-[0] == 0 )                   { $lprimertrim = max ( $+[0] , $lprimertrim ) }
		    elsif ( $+[0] == $seqlength-$ltrim-$rtrim ) { $rprimertrim = max ( $+[0] - $-[0] , $rprimertrim ) }
		    else  { die "Problem: $seq $re $seqlength $ltrim $rtrim @- @+ fuzzy primer regex didn't match an end\n" }                    ## sanity check
		    next PART ;
		}
	    }
	}

	# phase II: terminal nibbling
	foreach my $n ( 2 .. $partmax - $partmin ) {
	    foreach my $part ( @primerparts ) {
		if ( $$part{threeprime} ) {    # skip non-3' primer parts
		    if ( ($$part{dir} eq 'f' and $lprimertrim == 0) or ($$part{dir} eq 'r' and $rprimertrim == 0)) { #skip if end already primer-trimmed
			my $re = ${$$part{regex2}}{$n} ;
			if ($seq =~ $re ) {
			    $left = $ltrim + $-[0] + 1;
			    $right = $ltrim + $+[0] ;
			    push @procs, "p$$part{serial}\[$left-$right]";
			    if ( $debug ) {print "$serial\t$seq\tP$$part{serial}\[$left-$right]\t$re\t$-[0]\t$+[0]\t$lprimertrim\t$rprimertrim\n" } 
			    if      ( $-[0] == 0 )                         { $lprimertrim = max ( $+[0] , $lprimertrim ) if $lprimertrim == 0 }   ## left trim
			    elsif   ( $+[0] == $seqlength-$ltrim-$rtrim )  { $rprimertrim = max ( $+[0] - $-[0] , $rprimertrim ) if $rprimertrim == 0 }   ## right trim
			    else    { die "Problem: $seq $re $seqlength $ltrim $rtrim @- @+ nmer primer regex didn't match an end\n" }
			}
		    }
		}
	    }
	}

	if ( $lprimertrim != 0 or $rprimertrim != 0 ) {
	    $ltrim += $lprimertrim;
	    $rtrim += $rprimertrim;
	    if ($seqlength - $ltrim - $rtrim <= 0) { $seq = ""; $qual = ""; $currseqlength = 0; }
	    else {
		$seq  = substr $seq  , $lprimertrim , $seqlength - $ltrim - $rtrim ;         ######## sequence changing length here ########
		$qual = substr $qual , $lprimertrim , $seqlength - $ltrim - $rtrim ;
		$currseqlength = length $seq ;
	    }
	    if ( $currseqlength < $minLen ) {
		push @procs, 'R' ;
	    }
	}
    }

    #process D: check to see if DUSTing results in too short a sequence
    if ( $procs[-1] ne 'R' and not $dustoff ) {
	if (exists $dusthash{$header}) {
	    my $i;
	    for ($i=0; $i < (scalar @{$dusthash{$header}} - 1); ) {  ## looping through the dust-hash 2 array elements at a time (because dust coords always comes in pairs)
		$left = ${$dusthash{$header}}[$i] + 1;
		$right= ${$dusthash{$header}}[$i+1] ;
		push @procs, "D[$left-$right]";
		$dustsum += max( 0, min(${$dusthash{$header}}[$i+1],$seqlength-$rtrim) - max(${$dusthash{$header}}[$i],$ltrim)) ;
		$i += 2;
	    }
	    if ( $currseqlength - $dustsum < $minLen ) {
		push @procs, 'R' ;
	    }
	}
    }

    #process S: scoring overall quality
    $score = &$mean_qual_score(($qual,$currseqlength)) ;
    if ( $procs[-1] ne 'R' and $currseqlength != 0 and $score  < $qscore_threshold) { 
	push @procs , 'S' ;
	push @procs , 'R' ;
    }

#    $currseqlength = length $seq ;
#    if ($currseqlength != length $qual ) { print "$serial $header unequal seq $seq and qual $qual @procs!\n" }

    ## OUTPUT FASTQ and INDEX FILES HERE
    if ( $procs[-1] ne 'R' ) {
	if ( $nametype eq 'serial' ) { print FASTQ "\@$serial\n$seq\n+\n$qual\n" }
	else { print FASTQ "\@$header\n$seq\n+\n$qual\n" }
    }
    printf INDEX "$serial\t$header\t%u\t%u\t%.1f\t%.1f\t$ltrim\t$rtrim\t%s\n", ($seqlength, $currseqlength, $origscore, $score, join(";", @procs) );
}
close IN ;
close INDEX ;
close FASTQ ;

    
open INDEX, "<$prefix.index" ;
$ct = 0;
my $unrejected = 0;
my $length = 0;
my @line;
my @proc;
my %length;
my %qual ;
my %procs;
my %rejects ;
my %series ;
while (<INDEX>) {
    if ( $_ =~ m/^#/ ) {next}
    chomp $_ ;
    @line = split "\t", $_ ;
    $qual{int($line[4]+.5)}{before} ++ ;
    $length{$line[2]}{before} ++ ;
    $ct ++ ;
    @proc = split ";", $line[8];
    @proc = map { local $_ = $_; s/(\w*).*/$1/; $_ } @proc ;
    $series{join(",",@proc)} ++ ;
    foreach (@proc) {$procs{$_} ++ }
    if ($proc[-1] ne 'R') {
	$unrejected ++ ;
        $length += $line[3];
	$qual{int($line[5]+.5)}{after} ++ ;
	$length{$line[3]}{after} ++ ;
    }
    else { $rejects{$proc[-2]} ++ ; }
    #else { foreach (@proc) {$rejects{$_} ++ } }
}
close INDEX;

open ( SCORESUM, ">$prefix.scoresum" ) ;
print SCORESUM "# Call: $call\n" ;
print SCORESUM "Reads: $ct\n" ;
print SCORESUM "Unrejected: $unrejected\n" ;

if ($unrejected > 0) {print SCORESUM "Average Length of Unrejected: " , sprintf ( "%.1f" , $length / $unrejected ) , "\n" }
else                 {print SCORESUM "Average Length of Unrejected: 0\n" }

print SCORESUM "Quality Type: $guess\n" ;

print SCORESUM "\nProcesses per read\n" ;
foreach ( sort { $procs{$b} <=> $procs{$a} } keys %procs ) { 
    print SCORESUM "$_\t$procs{$_}\n" ; 
}

print SCORESUM "\nProcesses per reject\n" ;
foreach ( sort { $rejects{$b} <=> $rejects{$a} } keys %rejects ) { 
    print SCORESUM "$_\t$rejects{$_}\n" ;
}

print SCORESUM "\nProcess Series per read\n" ;
foreach ( sort { $series{$b} <=> $series{$a} } keys %series ) {
    print SCORESUM "$_\t$series{$_}\n" ;
}

print SCORESUM "\nAveQ\tBefore\tAfter\n" ;
foreach ( sort { $a <=> $b } keys %qual ) {
    print SCORESUM "$_\t$qual{$_}{before}\t$qual{$_}{after}\n" ;
}

print SCORESUM "\nLength\tBefore\tAfter\n" ;
foreach ( sort { $a <=> $b } keys %length ) {
    print SCORESUM "$_\t$length{$_}{before}\t$length{$_}{after}\n" ;
}

close SCORESUM;

   



# SUBROUTINES

sub mean_qual_score_0 { # obtain mean score of a sequence using sanger
    if ($_[1] == 0) {return 0}
    else {return (  sum(map(ord($_) - 33, split(//,$_[0] ) ) ) / $_[1] ) }
}
sub mean_qual_score_1 { # obtain mean score of a sequence using illumina>=1.3
    if ($_[1] == 0) {return 0}
    else {return (  sum(map(ord($_) - 64, split(//,$_[0] ) ) ) / $_[1] ) }
}
sub mean_qual_score_2 { # obtain mean score of a sequence using illumina<1.3
    if ($_[1] == 0) {return 0}
    else {return ( 10 * log(1 + 10 ** ((ord($_[0]) - 64) / 10.0)) / log(10)) }
} 

sub sum_qual_score_0 { # obtain mean score of a sequence using sanger
    if ($_[1] == 0) {return 0}
    else {return (  sum(map(ord($_) - 33, split(//,$_[0] ) ) ) ) }
}
sub sum_qual_score_1 { # obtain mean score of a sequence using illumina>=1.3
    if ($_[1] == 0) {return 0}
    else {return (  sum(map(ord($_) - 64, split(//,$_[0] ) ) ) ) }
}
sub sum_qual_score_2 { # obtain mean score of a sequence using illumina<1.3
    if ($_[1] == 0) {return 0}
    else {return ( 10 * log(1 + 10 ** ((ord($_[0]) - 64) / 10.0)) / log(10)) }
} 


sub Revcomp { # Compute the reverse complement of DNA sequence, even if degenerate
 my $dna = $_[0];
 my $revcomp = reverse $dna;
 $revcomp =~ tr/ACGTRYMKVBHDacgtrymkvbhd/TGCAYRKMBVDHtgcayrkmbvdh/; #N,S,W unchanged
 return $revcomp;
}

sub guess_qual_type {
 my $guessinfile=$_[0];
 if ( $verbose ) { print time-$t0, " ---> guessing quality (phred) type for file $guessinfile" }
 open ( IN , "-|" , $guessinfile ) or die "Cannot access file $guessinfile\n";
 while ( <IN> ) {
  $serial ++ ;
  <IN>; <IN>; # ignore all but qual
  $qual .= <IN>;
  chomp ( $qual ) ;
  if ( $serial == 100 ) { last }
 }
 my $ct = length $qual ;
 my ( %qual_chars , $avg ) ;
 for ( unpack ( '(a)*' , $qual ) ) { $qual_chars{$_} ++ }
 my $min = ord ( ( sort keys %qual_chars )[0] ) ;
 my $mode = ord ( ( sort { $qual_chars{$b} <=> $qual_chars{a} } keys %qual_chars )[0] ) ;
 for ( sort keys %qual_chars ) { 
  $avg += ord ( $_ ) * $qual_chars{$_} ;
 } 
 $avg = sprintf "%.1f" , $avg / $ct ;
 unless ( $qual_chars{'#'} ) { $qual_chars{'#'} = 0 }
 unless ( $qual_chars{'B'} ) { $qual_chars{'B'} = 0 }
 if ( $min < 64 ) { $guess = 'Phred+33 (-qual-type 0)' ; $qual_type = 0 }
 elsif ( $avg > 73 or $qual_chars{'B'} / $ct > 0.4 ) { $guess = 'Phred+64 (-qual-type 1)' ; $qual_type = 1 }
 else { $guess = 'Undetermined' }
 if ( $guess_qual or $verbose ) {
     $min .= ' (' . chr ( $min ) . ')' ;
     $mode .= ' (' . chr ( $mode ) . ')' ;
     $avg .= ' (' . chr ( int ( $avg ) ) . ')' ;
     if ($verbose) { print ", guess is: $guess\n" }
     if ($debug) { 
	 print "Mean ASCII= $avg\nMode ASCII= $mode\nMin ASCII= $min\nTotal Chars Counted = $ct\n" .
	 "Count for usual minimum score ASCII 35 (#) in phred+33 = $qual_chars{'#'}\n" .
	 "Count for usual minimum score ASCII 66 (B) in phred+64 = $qual_chars{'B'}\nGuess=$guess\n" ;
	 for ( sort keys %qual_chars ) { print ord ( $_ ) , "\t$_\t$qual_chars{$_}\n" }
     }
 }
 if ( $guess_qual ) { exit }
 if ( $qual_type == 3 ) { die "Can't guess qual_type. For fuller analysis run with -guess_qual option.\n" }
 close IN;
 return $qual_type ;
}


sub qtrim_simple {
    my $quality = shift;
    if ( $quality =~ m/$qualregex/ )  { return $-[0] + 1 }
    else                              { return 0         }
}

####################
sub qtrim_bwa {
    my $quality = shift;
    my $s = 0;
    my $max_qual = 0;
    my $max_pos = 0;
    my @qualityscores = map(ord($_) - $offset, split(//,$quality ) ) ;
    my $pos = length($quality) - 1;
    while ($pos >= 0 and $s >= 0) {
	$s += ( $endscore_threshold - $qualityscores[$pos] ) ;
	if ($s > $max_qual) { $max_qual = $s; $max_pos = $pos; }
	$pos--;
    }
    return $max_pos ;
}

## from http://wiki.bioinformatics.ucdavis.edu/index.php/TrimBWAstyle.pl
#my $h1;  my $s;  my $h2;  my $q;
#my $pos;  my $maxPos;  my $area;  my $maxArea;
#while ($h1 = <>) {  # read first header
#    $s = <>;  # read sequence
#    $q = <>;  # read quality scores
#    $pos = length($q);
#    $maxPos = $pos;
#    $area = 0;
#    $maxArea = 0;
#    while ($pos>0 and $area>=0) {
#	$area += $opt_q - (ord(substr($q,$pos-1,1))-33);
#	if ($area > $maxArea) {
#	    $maxArea = $area;
#	    $maxPos = $pos;
#	}
#	$pos--;
#    }
#    if ($pos==0) { $s = "N\n";  $q = "#\n" }  # scanned whole read and didn't integrate to zero?  replace with "empty" read ...
#    else {  # integrated to zero?  trim before position where area reached a maximum (~where string of qualities were still below 20 ...)
#	$s = substr($s,0,$maxPos)."\n";
#	$q = substr($q,0,$maxPos)."\n";
#    }
#    print $h1.$s.$h2.$q;
#}

