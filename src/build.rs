use bio::alphabets::dna;
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::fmindex::{FMIndex};
use bio::data_structures::suffix_array::suffix_array;
use bio::io::fasta;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;
use anyhow::Result;

use crate::cli::BuildArgs;

pub fn build_index(args: &BuildArgs) -> Result<()> {
    println!("Reading and concatenating sequences from {}", args.fasta);
    let concat = read_and_concat_fasta_with_spacer(&args.fasta, args.spacer_len)
        .expect("Failed to process FASTA");

    let alphabet = dna::n_alphabet();
    // 1. Create suffix array
    let sa = suffix_array(&concat);

    // 2. Build BWT from the text + suffix array
    let bwt = bwt(&concat, &sa);

    let less = less(&bwt, &alphabet);
    let occ = Occ::new(&bwt, args.sample_rate.try_into().unwrap(), &alphabet);

    let fm_index = FMIndex::new(&bwt, &less, &occ);

    // Save the FM-index using bincode (or keep it in memory if you prefer)
    let output_path = format!("{}.fm", args.index);
    let file = File::create(&output_path).expect("Failed to create index file");
    bincode::serialize_into(BufWriter::new(file), &fm_index)
        .expect("Failed to write FM-index");

    println!("FM-index saved to {}", output_path);
    Ok(())
}

/// Helper to read FASTA and join sequences with spacer Ns
fn read_and_concat_fasta_with_spacer<P: AsRef<Path>>(
    path: P,
    spacer_len: usize,
) -> std::io::Result<Vec<u8>> {
    let reader = fasta::Reader::new(BufReader::new(File::open(path)?));
    let spacer = vec![b'N'; spacer_len];

    let mut concatenated = Vec::new();
    let mut first = true;

    for result in reader.records() {
        let record = result?;
        if !first {
            concatenated.extend_from_slice(&spacer);
        }
        concatenated.extend_from_slice(record.seq());
        first = false;
    }
    // append a lexicographically small sentinel symbol $ to the end of the concatenated sequence.
    concatenated.push(b'$');
    println!("Concatenated sequence length: {}", concatenated.len());
    // println!("Concatenated sequence: {:?}", std::str::from_utf8(&concatenated).unwrap());
    Ok(concatenated)
}
