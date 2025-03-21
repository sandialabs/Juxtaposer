# Juxtaposer

## INTRO
Juxtaposer looks for NGS reads that may have resulted from genomic recombination, and tests them for the
possibility of being: circles/scars from mobile elements, transpositions of transposons, or a palindrome artifact.

## AUTHORS
Kelly Williams, kpwilli@sandia.gov
Corey Hudson

Repo maintained by: Ellis Torrance, eltorra@sandia.gov

## CITATIONS
If using Juxtaposer in your research please cite:

[Schoeniger JS, Hudson CM, Bent ZW, Sinha A, Williams KP. Experimental single-strain mobilomics reveals events that shape pathogen emergence. Nucleic Acids Res. 2016 Aug 19;44(14):6830-9. doi: 10.1093/nar/gkw601. Epub 2016 Jul 4. PMID: 27378783; PMCID: PMC5001619.](https://academic.oup.com/nar/article/44/14/6830/2468219)

Please also cite the associate publications of the Juxtaposer dependencies.

## SOFTWARE DEPENDENCIES
 - bowtie2 and bowtie2-build in path
 - blastn and makeblastdb in path
 - hmmscan in path
 - prodigal: binary is included, but if the binary fails on your machine, replace with fresh install

Instutions for installing these dependencies in a conda environment is available in INSTALLATION below. 

## INSTALLATION

1. Download the github repo from https://github.com/sandialabs/Juxtaposer or clone it with git:
```
git clone https://github.com/sandialabs/Juxtaposer.git
```
2. Make Juxtaposer files callable system wide by adding the following to your bash profile. Then, restart the terminal session or use source to reload your bash profile 
```
export PATH=<absolute path to juxtaposer folder>:$PATH
```
3. Install Juxtaposer dependencies. See below for instructions on installation of these dependencies using conda.

### (OPTIONAL) Creating a Juxtaposer [Conda](https://www.anaconda.com/docs/getting-started/miniconda/install) environment
1. Create the Conda environment
```
conda create --name Juxtaposer -y
```
2. activate the Conda environment (Note, you will need to activate the conda environment every time before you run Juxtaposer)
```
conda activate Juxtaposer
```
3. install juxtaposer dependencies
```
conda config --add channels bioconda && conda install -y bowtie2 && conda install -y blast && conda install -y hmmer && conda install -y pftools && conda install -y bedtools
```

## INPUT DIRECTORIES AND FILES

Juxtaposer requires an initial directory containing a single fastq file named <ins>qf.fq</ins> (or qf.fq.gz if gunzipped) and the corresponding assembled genome in fasta format named <ins>genome.fa</ins> where all output files will be written. For the purpose of the readme, this folder will be called <ins>"ngsdir"</ins>.
We recommend using a different ngsdir directory for each experimental sample.

### PREPARING READS

Juxtaposer requires an input file qf.fq, that contains NGS data from which low quality terminal regions and library
primer sequences have been trimmed. Primer trimming allows the standard-read filter's global search to correctly
identify standard genomic reads. We supply the quality filter qfilter.pl that we use. To use qfilter.pl properly,
you must supply it the primer sequences used for library preparation. Qfilter.pl reads a file called primers.txt,
that already contains some commonly used primer sets. To see those:
$ qfilter.pl -primersets_available
If your library used a primer set already in primers.txt simply supply that primerset name with the -primers
option. If your library used different primers, you can write the set into the primers.txt file, or specify them
directly with the -primers option.
We have been using SE reads, but if you have PE reads, we recommend the following: merge pairs using a tool like
bbduk or PEAR, aiming to remove primer sequences and low quality regions, and combine all merged read pairs and
unmerged reads into the qf.fq input file.

## RUNNING JUXTAPOSER
Juxtaposer must be run from the ngsdir. We also recomend sending print statements to a log file you can send to us in case troubleshooting is needed.
```
cd <path to ngsdir>; perl <path to juxtaposer>juxtaposer.pl genome &> log.txt
```

## OUTPUTS
An example of the files that are created by Juxtaposer are in the directory Test. You can also test run juxtaposer is installed correctly with all necessary dependencies by running it on this folder.

Here is a description of files and folders created by juxtaposer:
refdir: Directory containing information about the reference genome, to which Juxtaposer will write some files.
refPfx: refDir/prefix, where prefix is the part of the reference genome file preceding the .fa suffix
refPfx.fa: Multi-fasta file with all replicons (dnas) in the genome
[refPfx.tn.bed]: Recommended but optional transposon annotation file. Precise coordinates of genome's
 transposons/ISs, used to annotate reads from potential transposition junctions. An initial run of Juxtaposer may
 help discover circularizing transposons to add to this list. The user can also add other mobile elements, such as
 genomic islands, self-splicing introns, integron cassettes, etc., that may generate circular DNAs or deletion scars.
 Format, tab-delimited columns: [1] replicon; [2] left coordinate minus one; [3] right coordinate; [4] element
  name; [5] 1 if transposable, otherwise 0; [6] orientation (+ or -)

refPfx.lens: Two tab-delimited columns: dnaName, Length (basepairs).
refPfx_clos_n.fa: The input ref.fa to which has been added closure sequences for each circular replicon, based
 on n, the length of the longest read in qf.fq
refPfx: bowtie2 index and blast library for refPfx_clos_n.fa, many mobility gene outputs
refPfx.mobile.bed: List of mobility enzyme genes (integrases, transposases, etc.)
juxtas.txt: see below for field descriptions

JUXTAS.TXT FIELDS
dnaLo: replicon containing the lower-position hit
CoordLo: junction coordinate for the lower-position hit
dnaHi: replicon containing the higher-position hit
CoordHi: junction coordinate for the higher-position hit
config: orientations of Lo and Hi hits
call: possibility of circularJunction, scar, palindrome or transposition end
overlap: length of overlap (negative if unmatched spacer)
origin: 0 if different DNAs, 1 if shortest distance around DNA not spanning origin, -1 if spanning origin
idLo: percent identity of Lo hit
idHi: percent identity of Hi hit
tnLo: transposon end in Lo hit
tnOffsetLo: distance from transposon end to Lo hit end
tnHi: transposon end in Lo hit
tnOffsetHi: distance from transposon end to Lo hit end
sampleRead: read ID
overlapSeq: sequence of slash-marked overlap, with 5 bp flanks. Uppercase, hit Lo (on left)
readseq: sequence of read. Uppercase, hit Lo (on left). Slash-marked, hit Hi
seqLo: Mismatch-marked hit Lo
seqHi: Mismatch-marked hit Hi
count: Number of similar reads
sites: Number of additional identical-scoring hit pairs for same read
Mob: mobility gene spanned by hitpair
pctMob: percentage of the mobility gene spanned
shortDistMob: closest distance between mobility gene and hitpair span end
