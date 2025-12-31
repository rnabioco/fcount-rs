#!/bin/bash
# Benchmark script for fcount-rs vs featureCounts
#
# Usage: ./benchmark.sh [--regenerate] [--build]
#
# Options:
#   --regenerate  Regenerate chr22 test files from full dataset
#   --build       Force rebuild of release binary before benchmarking

set -e

# Paths
FCOUNT_RS="./target/release/fcount-rs"
TEST_GTF="/tmp/test_chr22.gtf"
TEST_BAM="/tmp/test_chr22.bam"
TEST_OUTPUT="/tmp/fcount_benchmark.txt"

# Full dataset paths (for regenerating test files)
FULL_GTF="/beevol/home/jhessel/devel/rnabioco/nmd-exons/results/assembly/merged/merged.gtf"
FULL_BAM="/beevol/home/jhessel/devel/rnabioco/nmd-exons/results/align/bam/HELA-smg1i-3.bam"

# Check if we need to regenerate test files
regenerate_tests() {
    echo "Regenerating test files from chr22..."
    grep "^chr22" "$FULL_GTF" > "$TEST_GTF"
    samtools view -b "$FULL_BAM" chr22 > "$TEST_BAM"
    samtools index "$TEST_BAM"
    echo "Created:"
    echo "  GTF: $TEST_GTF ($(wc -l < $TEST_GTF) lines)"
    echo "  BAM: $TEST_BAM ($(samtools view -c $TEST_BAM) reads)"
}

# Build release binary
build() {
    echo "Building release binary..."
    cargo build --release 2>&1 | grep -E "(Compiling|Finished|error)"
}

# Check binary exists
check_binary() {
    if [[ ! -x "$FCOUNT_RS" ]]; then
        echo "Binary not found at $FCOUNT_RS"
        echo "Run with --build flag or: cargo build --release"
        exit 1
    fi
}

# Run benchmarks
run_benchmarks() {
    local gtf=$1
    local bam=$2

    echo ""
    echo "=== BENCHMARK: chr22 subset (~1.4M reads) ==="
    echo "GTF: $gtf"
    echo "BAM: $bam"
    echo ""

    # Check if hyperfine is available
    if ! command -v hyperfine &> /dev/null; then
        echo "hyperfine not found, install with: cargo install hyperfine"
        exit 1
    fi

    # Use pixi run for featureCounts if available
    local FC_CMD="featureCounts"
    if ! command -v featureCounts &> /dev/null; then
        if [[ -f "pixi.toml" ]]; then
            FC_CMD="pixi run featureCounts"
        else
            echo "featureCounts not found"
            echo "Running fcount-rs only..."
            hyperfine --warmup 1 --runs 5 \
                -n "fcount-rs-1t" "$FCOUNT_RS -a $gtf -o $TEST_OUTPUT -p $bam -T 1" \
                -n "fcount-rs-8t" "$FCOUNT_RS -a $gtf -o $TEST_OUTPUT -p $bam -T 8"
            return
        fi
    fi

    echo "Comparing fcount-rs vs featureCounts..."

    hyperfine --warmup 1 --runs 5 \
        -n "featureCounts-1t" "$FC_CMD -a $gtf -o /tmp/fc_out.txt -p -T 1 $bam 2>/dev/null" \
        -n "featureCounts-8t" "$FC_CMD -a $gtf -o /tmp/fc_out.txt -p -T 8 $bam 2>/dev/null" \
        -n "fcount-rs-1t" "$FCOUNT_RS -a $gtf -o $TEST_OUTPUT -p $bam -T 1" \
        -n "fcount-rs-8t" "$FCOUNT_RS -a $gtf -o $TEST_OUTPUT -p $bam -T 8" \
        --export-markdown /tmp/benchmark_results.md

    echo ""
    cat /tmp/benchmark_results.md
}

# Main
main() {
    local do_build=false

    for arg in "$@"; do
        case $arg in
            --regenerate)
                regenerate_tests
                exit 0
                ;;
            --build)
                do_build=true
                ;;
        esac
    done

    # Check test files exist
    if [[ ! -f "$TEST_GTF" || ! -f "$TEST_BAM" ]]; then
        echo "Test files not found, regenerating..."
        regenerate_tests
    fi

    # Build only if requested
    if $do_build; then
        build
    else
        check_binary
    fi

    run_benchmarks "$TEST_GTF" "$TEST_BAM"
}

main "$@"
