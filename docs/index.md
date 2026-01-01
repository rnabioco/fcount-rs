# fcount-rs

Ultrafast feature counting for RNA-seq data. A high-performance Rust alternative to [featureCounts](http://subread.sourceforge.net/).

## Performance

fcount-rs is **2× faster** than featureCounts with 8 threads on typical RNA-seq data:

| Tool | Threads | Time (s) | Speedup |
|:---|:---:|---:|:---:|
| featureCounts | 1 | 3.34 | — |
| featureCounts | 8 | 1.45 | — |
| fcount-rs | 1 | 2.26 | 1.5× |
| fcount-rs | 8 | 0.71 | **2.0×** |

*Benchmark: chr22 subset (~1.4M paired-end reads). featureCounts v2.1.1.*

## Features

- **Fast**: 2× faster than featureCounts with parallel processing, custom GTF parser, and efficient data structures
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
