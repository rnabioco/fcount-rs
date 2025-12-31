# fcount-rs

Ultrafast feature counting for RNA-seq data. A high-performance Rust alternative to [featureCounts](http://subread.sourceforge.net/).

## Features

- **Fast**: Optimized for speed with parallel processing, custom GTF parser, and efficient data structures
- **Pure Rust**: No C dependencies, uses the [noodles](https://github.com/zaeleus/noodles) library for BAM/SAM parsing
- **Compatible**: Drop-in replacement for featureCounts with similar CLI interface
- **Flexible**: Supports single-end and paired-end data, stranded protocols, and multi-mapping reads

## Quick Start

```bash
# Install
cargo install fcount-rs

# Count reads
fcount-rs -a annotation.gtf -o counts.txt sample.bam
```

## Getting Help

```bash
fcount-rs --help
```
