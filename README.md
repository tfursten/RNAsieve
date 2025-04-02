# RNAsieve
**RNAsieve** is a Rust command-line tool for filtering RNA sequencing reads based on seed-level matches against a prebuilt FM-index. It supports FASTA, FASTQ, and gzipped inputs, and can handle both single-end and paired-end reads.

## Features

- Seed-based filtering using an FM-index
- Paired-end and single-end read support
- Input formats: FASTA, FASTQ, `.gz` compressed files
- Outputs matched and unmatched reads to separate files
- Automatically detects reverse complement seed matches
- Built for speed using Rust and `needletail`/`rust-bio`

---

## 📦 Installation

Clone the repository and build with cargo:

```bash
git clone https://github.com/yourusername/rnasieve.git
cd rnasieve
cargo build --release
```

The resulting binary will be in `target/release/rnasieve`.

---

## Usage

### Build an index

```bash
rnasieve build-index \
  --fasta reference_amplicons.fasta \
  --output index.fm
```

### Filter reads

```bash
rnasieve filter \
  --index index.fm \
  --read1 reads_R1.fastq.gz \
  --read2 reads_R2.fastq.gz \
  --cutoff 5 \
  --seed-size 20 \
  --seed-interval 5
```

By default, the output files will be:
- `matched_reads_R1.fastq.gz`, `matched_reads_R2.fastq.gz` For reads that matched index 
- `filtered_reads_R1.fastq.gz`, `filtered_reads_R2.fastq.gz` For reads that did not match index

You can override output prefixes or formats with:
- `--matched-prefix`
- `--filtered-prefix`
- `--output-format` (`fastq`, `fastq.gz`, `fasta`, `fasta.gz`)

---

## ⚙️ Options

| Option               | Description                                         |
|----------------------|-----------------------------------------------------|
| `--index`            | Path to FM-index file                               |
| `--read1`            | Required: path to R1 reads                          |
| `--read2`            | Optional: path to R2 reads                          |
| `--cutoff`           | Minimum number of seed matches to consider a "hit"  |
| `--seed-size`        | Length of each seed                                 |
| `--seed-interval`    | Spacing between seeds in the read                   |
| `--output-format`    | Output format: `fastq`, `fasta`, `fastq.gz`, etc.   |
| `--matched-prefix`   | Custom prefix for matched reads output              |
| `--filtered-prefix`  | Custom prefix for filtered reads output             |


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
```

---

