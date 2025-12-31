# fcount-rs

[![CI](https://github.com/rnabioco/fcount-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/rnabioco/fcount-rs/actions/workflows/ci.yml)
[![Docs](https://img.shields.io/badge/docs-mkdocs-blue)](https://rnabioco.github.io/fcount-rs)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Ultrafast feature counting for RNA-seq data. A high-performance Rust alternative to featureCounts.

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
