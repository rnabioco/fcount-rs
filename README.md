# fcount-rs

[![CI](https://github.com/rnabioco/fcount-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/rnabioco/fcount-rs/actions/workflows/ci.yml)
[![Docs](https://img.shields.io/badge/docs-mkdocs-blue)](https://rnabioco.github.io/fcount-rs)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Ultrafast feature counting for RNA-seq data. A high-performance Rust alternative to featureCounts.

## Performance

fcount-rs is significantly faster than featureCounts, especially with multiple threads:

| Command | Mean [s] | Relative |
|:---|---:|---:|
| featureCounts (1 thread) | 3.34 ± 0.25 | 4.69× slower |
| featureCounts (8 threads) | 1.45 ± 0.03 | 2.03× slower |
| fcount-rs (1 thread) | 2.26 ± 0.34 | 3.18× slower |
| **fcount-rs (8 threads)** | **0.71 ± 0.02** | **1.00** |

*Benchmark: chr22 subset (~1.4M paired-end reads), measured with [hyperfine](https://github.com/sharkdp/hyperfine). featureCounts v2.1.1.*

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
fcount-rs -a annotation.gtf -o counts.txt sample.bam
```

For paired-end data:

```bash
fcount-rs -p -a annotation.gtf -o counts.txt sample.bam
```

See the [full documentation](https://rnabioco.github.io/fcount-rs) for all options and examples.

## License

MIT
