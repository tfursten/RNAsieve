# RNAsieve
**RNAsieve** is a Rust command-line tool for filtering RNA sequencing reads based on seed-level matches against a prebuilt FM-index. It supports FASTA, FASTQ, and gzipped inputs, and can handle both single-end and paired-end reads.

## Features

- Seed-based filtering using an FM-index
- Paired-end and single-end read support
- Input formats: FASTA, FASTQ, `.gz` compressed files
- Outputs matched and unmatched reads to separate files
- Checks forward and reverse complement for each read
- Built for speed using Rust and `needletail`/`rust-bio`

---

## 📦 Installation

Clone the repository and build with cargo:

```bash
git clone https://github.com/tfursten/rnasieve.git
cd RNAsieve
cargo build --release
```

The resulting binary will be in `target/release/RNAsieve`.

---

## Usage

### Build an index

```bash
RNAsieve build-index \
  --fasta reference_amplicons.fasta \
  --output index.fm
```

### Filter reads

```bash
RNAsieve filter \
  index.fm results \
  --read1 reads_R1.fastq.gz \
  --read2 reads_R2.fastq.gz \
  --cutoff 5 \
  --seed-size 20 \
  --seed-interval 5
```

The output files will be prefixed with `matched_` and `filtered_`:
- `matched_reads_R1.fastq.gz`, `matched_reads_R2.fastq.gz` For reads that matched index 
- `filtered_reads_R1.fastq.gz`, `filtered_reads_R2.fastq.gz` For reads that did not match index

You can override output formats with:
- `--output-format` (`fastq`, `fastq.gz`, `fasta`, `fasta.gz`)
- Dummy quality scores will be added if converting from fasta to fastq

---

## ⚙️ Options

| Option               | Description                                         |
|----------------------|-----------------------------------------------------|
| `index`              | Path to FM-index file                               |
| `outdir`             | Output directory                                    |
| `--read1`            | Required: path to R1 reads                          |
| `--read2`            | Optional: path to R2 reads                          |
| `--cutoff`           | Minimum number of seed matches to consider a "hit"  |
| `--seed-size`        | Length of each seed                                 |
| `--seed-interval`    | Spacing between seeds in the read                   |
| `--output-format`    | Output format: `fastq`, `fasta`, `fastq.gz`, etc.   |



---

## 📚 Dependencies

- [rust-bio](https://crates.io/crates/bio)
- [needletail](https://crates.io/crates/needletail)
- [flate2](https://crates.io/crates/flate2)
- [anyhow](https://crates.io/crates/anyhow)
- [clap](https://crates.io/crates/clap)

---

## 📄 License

MIT License. See `LICENSE` file for details.

---

## 👩‍💻 Author

Maintained by Tara Furstenau.


