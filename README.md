# fcount-rs

[![CI](https://github.com/rnabioco/fcount-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/rnabioco/fcount-rs/actions/workflows/ci.yml)
[![Docs](https://img.shields.io/badge/docs-mkdocs-blue)](https://rnabioco.github.io/fcount-rs)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**🚧 fcount-rs is under active development. Caveat emptor. 🚧**

Ultrafast feature counting for RNA-seq data. A high-performance Rust alternative to featureCounts.

## Performance

fcount-rs is significantly faster than featureCounts, especially with multiple threads:

| Command | Mean [s] | Relative |
|:---|---:|---:|
| featureCounts (1 thread) | 2.41 ± 0.01 | 7.67× slower |
| featureCounts (8 threads) | 0.84 ± 0.01 | 2.68× slower |
| fcount-rs (1 thread) | 1.16 ± 0.02 | 3.70× slower |
| **fcount-rs (8 threads)** | **0.31 ± 0.01** | **1.00** |

*Benchmark: chr22 subset (~1.4M paired-end reads), measured with [hyperfine](https://github.com/sharkdp/hyperfine). featureCounts v2.0.6, Subread.*

## Installation

```bash
cargo install fcount-rs
```

Or build from source:

```bash
git clone https://github.com/rnabioco/fcount-rs
cd fcount-rs
cargo build --release
```

## Usage

```bash
# Basic usage
fcount -a annotation.gtf -o counts.txt sample.bam

# Paired-end mode
fcount -a annotation.gtf -o counts.txt -p sample.bam

# Multiple samples with custom names
fcount -a annotation.gtf -o counts.txt -p control=ctrl.bam treated=treat.bam

# With 8 threads, strand-specific, primary alignments only
fcount -a annotation.gtf -o counts.txt -p -t 8 -s 2 --primary sample.bam
```

### Key Options

```
Input/Output:
  -a, --annotation <GTF>    GTF/GFF annotation file
  -o, --output <FILE>       Output count matrix file

Paired-End:
  -p, --paired              Count fragments instead of reads
  -B, --both-aligned        Require both ends aligned
  -C, --no-chimeric         Exclude chimeric fragments

Filtering:
  -Q, --min-mapq <INT>      Minimum mapping quality [default: 0]
  --primary                 Count primary alignments only
  --ignore-dup              Ignore duplicate reads

Strandedness:
  -s, --strand <MODE>       0=unstranded, 1=stranded, 2=reversely stranded [default: 0]

Multi-mapping:
  -M, --multi-mapping       Count multi-mapping reads
  --fraction                Fractional counting (1/NH) for multi-mappers

Overlap:
  -O, --multi-overlap       Allow reads to overlap multiple features
  --min-overlap <INT>       Minimum overlapping bases [default: 1]
  --largest-overlap         Assign to feature with largest overlap

Performance:
  -t, --threads <INT>       Number of threads (0 = auto) [default: 0]
```

See the [full documentation](https://rnabioco.github.io/fcount-rs) for all options and examples.

## License

MIT
