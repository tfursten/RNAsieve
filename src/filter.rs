use bio::data_structures::fmindex::{FMIndex, FMIndexable, BackwardSearchResult, Interval};
use bio::data_structures::bwt::{Occ};
use needletail::{parse_fastx_file, Sequence};
use std::fs::File;
use std::path::{Path, PathBuf};
use std::io::{BufReader, Write, BufWriter};
use anyhow::Result;
use bincode;

use flate2::write::GzEncoder;
use flate2::Compression;

use crate::cli::FilterArgs;

fn ensure_leading_dot(ext: &str) -> String {
    if ext.starts_with('.') {
        ext.to_string()
    } else {
        format!(".{}", ext)
    }
}



fn build_output_path(input_path: &str, prefix: &str, output_format: &str) -> PathBuf {
    let path = Path::new(input_path);
    let parent = path.parent().unwrap_or_else(|| Path::new(""));
    let stem = path.file_stem().unwrap().to_string_lossy();
    let filename = path.file_name().unwrap().to_string_lossy();

    let output_filename = if output_format.is_empty() {
        format!("{}{}", prefix, filename)
    } else {
        format!("{}{}{}", prefix, stem, ensure_leading_dot(output_format))
    };

    parent.join(output_filename)
}

fn infer_format_from_filename(path: &str) -> &str {
    let path = path.to_ascii_lowercase();
    if path.ends_with(".fastq") || path.ends_with(".fastq.gz") || path.ends_with(".fq") || path.ends_with(".fq.gz") {
        "fastq"
    } else if path.ends_with(".fasta") || path.ends_with(".fasta.gz") || path.ends_with(".fa") || path.ends_with(".fa.gz") {
        "fasta"
    } else {
        "fastq" // fallback
    }
}


pub fn filter(args: &FilterArgs) -> Result<()> {
    let index_file = File::open(&args.index)?;
    let reader = BufReader::new(index_file);

    // You'll need to know your FMIndex types – here's a common default setup:
    let fm_index: FMIndex<Vec<u8>, Vec<usize>, Occ> = bincode::deserialize_from(reader)?;


    let paired = args.read2.is_some();

    let mut reader1 = parse_fastx_file(&args.read1)?;
    let mut reader2 = if paired {
        Some(parse_fastx_file(args.read2.as_ref().unwrap())?)
    } else {
        None
    };
    
    
    let format = match &args.output_format {
        Some(fmt) => fmt.as_str(), // use user-supplied format
        None => {
            // Infer from read1 file extension
            infer_format_from_filename(&args.read1)
        }
    };
    
    let matched_r1_path = build_output_path(&args.read1, &args.matched_prefix, format);
    let matched_r2_path = args.read2.as_ref().map(|r2| {
        build_output_path(r2, &args.matched_prefix, format)
    });
    
    let filtered_r1_path = build_output_path(&args.read1, &args.filtered_prefix, format);
    let filtered_r2_path = args.read2.as_ref().map(|r2| {
        build_output_path(r2, &args.filtered_prefix, format)
    });


    let mut matched_r1 = create_writer(matched_r1_path.to_str().unwrap())?;
    let mut matched_r2 = if paired {
        Some(create_writer(matched_r2_path.unwrap().to_str().unwrap())?)
    } else {
        None
    };
    
    let mut filtered_r1 = create_writer(filtered_r1_path.to_str().unwrap())?;
    let mut filtered_r2 = if paired {
        Some(create_writer(filtered_r2_path.unwrap().to_str().unwrap())?)
    } else {
        None
    };


    

    while let Some(record1) = reader1.next() {
        let rec1 = record1?;
        let id = String::from_utf8_lossy(rec1.id());
        let seq = rec1.seq();
        let qual = rec1.qual();
    
        let mut r2_seq = None;
        let mut r2_id = None;
        let mut r2_qual = None;
    
        let r1_hits = seed_hits(&rec1.seq(), args, &fm_index)?;
    
        let r2_hits = if let Some(reader2) = reader2.as_mut() {
            if let Some(record2) = reader2.next() {
                let rec2 = record2?;
                r2_seq = Some(rec2.seq().to_vec());
                r2_id = Some(String::from_utf8_lossy(rec2.id()).to_string());
                r2_qual = rec2.qual().map(|q| q.to_vec());

                seed_hits(&rec2.seq(), args, &fm_index)?
            } else {
                0
            }
        } else {
            0
        };
    
        if r1_hits < args.cutoff && r2_hits < args.cutoff {
            write_record(&id, &seq, qual.as_deref(), &mut filtered_r1, format)?;
            if let (Some(r2_id), Some(r2_seq)) = (&r2_id, &r2_seq) {
                write_record(
                    r2_id,
                    r2_seq,
                    r2_qual.as_deref(),
                    filtered_r2.as_mut().unwrap(),
                    format,
                )?;
            }
        } else {
            write_record(&id, &seq, qual.as_deref(), &mut matched_r1, format)?;
            if let (Some(r2_id), Some(r2_seq)) = (&r2_id, &r2_seq) {
                write_record(
                    r2_id,
                    r2_seq,
                    r2_qual.as_deref(),
                    matched_r2.as_mut().unwrap(),
                    format,
                )?;
            }
        }
    }

    Ok(())
}




/// Extracts seeds from a sequence and reverse complement using the specified length and interval.
fn get_seeds_both_strands(seq: &[u8], seed_len: usize, interval: usize) -> Vec<Vec<u8>> {
    let mut seeds = Vec::new();

    // Forward strand
    let mut pos = 0;
    while pos + seed_len <= seq.len() {
        seeds.push(seq[pos..pos + seed_len].to_vec());
        pos += interval;
    }

    // Reverse complement strand
    let rc = seq.reverse_complement(); // calls the Sequence trait method
    pos = 0;
    while pos + seed_len <= rc.len() {
        seeds.push(rc[pos..pos + seed_len].to_vec());
        pos += interval;
    }

    seeds
}


fn seed_hits(
    seq: &[u8],
    args: &FilterArgs,
    fm_index: &FMIndex<Vec<u8>, Vec<usize>, Occ>,
) -> Result<usize> {
    let seeds = get_seeds_both_strands(seq, args.seed_size, args.seed_interval);
    let mut count = 0;

    for seed in &seeds {
        let interval = fm_index.backward_search(seed.iter());

        let mut interval_upper = 0;
        let mut interval_lower = 0;
        match interval {
            BackwardSearchResult::Complete(sai) => {
                interval_upper = sai.upper;
                interval_lower = sai.lower;
                sai
            }
            BackwardSearchResult::Partial(sai, _l) => { 
                sai
            }
            BackwardSearchResult::Absent => {
                Interval {
                    upper: 0,
                    lower: 0
                }
            }
        };
        // If no interval is returned no seed hits were found                 
        if (interval_upper == 0) && (interval_lower == 0) {
            continue;
        }
        let n_hits = interval_upper - interval_lower;
        if n_hits > 0 {
            count += 1;
        }
        // Optional: if you want positions from suffix array:
        // let positions = result.occ(&suffix_array);  <-- Only valid if you pass it in
    }

    Ok(count)
}

fn create_writer(path: &str) -> Result<Box<dyn Write>, std::io::Error> {
    let file = File::create(path)?;
    let writer: Box<dyn Write> = if path.ends_with(".gz") {
        Box::new(GzEncoder::new(file, Compression::default()))
    } else {
        Box::new(BufWriter::new(file))
    };
    Ok(writer)
}




fn write_record(
    id: &str,
    seq: &[u8],
    qual: Option<&[u8]>,
    writer: &mut dyn Write,
    output_format: &str,
) -> Result<()> {
    let seq_str = String::from_utf8_lossy(seq);

    match output_format {
        "fastq" | "fastq.gz" => {
            let qual_str = qual
                .map(|q| String::from_utf8_lossy(q).to_string())
                .unwrap_or_else(|| "I".repeat(seq.len()));
            writeln!(writer, "@{}\n{}\n+\n{}", id, seq_str, qual_str)?;
        }
        "fasta" | "fasta.gz" => {
            writeln!(writer, ">{}\n{}", id, seq_str)?;
        }
        _ => return Err(anyhow::anyhow!("Unsupported output format: {}", output_format)),
    }

    Ok(())
}



