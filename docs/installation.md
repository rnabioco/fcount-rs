# Installation

## From crates.io

```bash
cargo install fcount-rs
```

## From source

```bash
git clone https://github.com/rnabioco/fcount-rs
cd fcount-rs
cargo build --release
```

The binary will be at `target/release/fcount-rs`.

## Requirements

- Rust 1.70+ (uses edition 2021)
- Linux or macOS (uses jemalloc on non-MSVC platforms)
