use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(name = "rna_sieve")]
#[command(about = "Tool to filter out rRNA and host reads using FM-index and q-gram filtering.", long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    #[command(about = "Build an FM-index from a reference sequences")]
    BuildIndex(BuildArgs),
    #[command(about = "Filter reads that match a reference sequence")]
    Filter(FilterArgs),
}

#[derive(clap::Args)]
pub struct BuildArgs {
    pub fasta: String,
    pub index: String,
    #[arg(long, short='n', default_value_t = 50, help="Length of the spacer between sequences in the index.")]
    pub spacer_len: usize,
    #[arg(long, short='f', default_value_t = 64, help="BWT occurance sampling rate. If sample interval is k, every k-th entry will be kept.")]
    pub sample_rate: usize,
}

/// Arguments for filtering reads

#[derive(clap::Args)]
pub struct FilterArgs {
    pub index: String,
    pub outdir: String,
    /// Required R1 reads (fastq, fastq.gz, or fasta).
    #[arg(long, short = '1')]
    pub read1: String,
    /// Optional R2 reads (fastq, fastq.gz, or fasta).
    #[arg(long, short = '2')]
    pub read2: Option<String>,
    /// Size of the seeds to extract from the reads.
    #[arg(long, short='s', default_value_t = 20)]
    pub seed_size: usize,
    /// Interval/spacing between seeds in the reads.
    #[arg(long, short='i', default_value_t = 5)]
    pub seed_interval: usize,
    /// Minimum number of seeds that must match the index for a read to be considered a hit.
    #[arg(long, short='c', default_value_t = 1)]
    pub cutoff: usize,
    /// Output format for the reads. Options: 'fastq', 'fasta', 'fastq.gz' [default: matches input format].
    /// Dummy quality scores added for fastq if not present in the input.
    #[arg(long)]
    pub output_format: Option<String>,
}

