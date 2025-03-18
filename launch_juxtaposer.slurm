#!/bin/bash
#SBATCH --job-name=juxtaposer_job      # Job name
#SBATCH --time=01:00:00                  # Time limit hrs:min:sec
#SBATCH --ntasks=1                       # Number of tasks
#SBATCH --mem=4G                         # Memory limit. YOU MAY NEED TO CHANGE THIS

source activate Juxtaposer

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -reads) reads="$2"; shift ;;
        -reference) reference="$2"; shift ;;
        -outdir) outdir="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done
#-readlen) readlen="$2"; shift ;;

# Check if all required arguments are provided
if [[ -z "$reads" || -z "$reference" || -z "$outdir" ]]; then
    echo "Usage: script.slurm -reads <path to reads> -reference <path to reference> -outdir <path to outdir>"
    exit 1
fi

# Create output directory
mkdir -p "$outdir"

# Create aliases for reads and reference
ln -s "$reads" "$outdir/qf.fq"
ln -s "$reference" "$outdir/genome.fa"

# Navigate to the output directory and run the program
cd "$outdir" || exit
#juxtaposer.pl genome &> log
juxtaposer.pl genome &> log

# End of script

