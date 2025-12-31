# Performance Regression Investigation

## Current State (BROKEN)
- **fcount-rs**: 71 seconds on full dataset (52M reads)
- **featureCounts**: 23 seconds on same dataset
- **Ratio**: 3x slower than featureCounts

## What We Know
1. GTF parsing: 14s (vs ~10s for featureCounts)
2. BAM counting: 55s (vs ~13s for featureCounts) 
3. Parallelism achieved: 1.55x (110s user / 71s real) with 8 threads

## Changes Made (Potentially Problematic)
1. Added per-chromosome parallelism with IndexedBamReader
2. Switched between coitrees and rust-lapper (no difference)
3. Added streaming vs Vec collection (no difference)
4. Modified thread allocation for BGZF
5. Added buffered overlap lookups

## Key Questions to Answer
1. What was the ACTUAL original baseline performance?
2. Is the MultithreadedReader being used correctly?
3. Where is time actually spent? (Need proper profiling)
4. Why is parallelism so low (1.55x with 8 threads)?

## First Principles Analysis Needed
1. **Benchmark raw noodles BAM reading** - just iterate records, no processing
2. **Benchmark interval queries alone** - with synthetic data
3. **Profile with flamegraph** - find actual hotspots
4. **Compare with a minimal implementation** - strip down to essentials

## Hypothesis
The sequential MultithreadedReader approach should work well because:
- 8 threads decompress BGZF blocks in parallel
- Main thread processes records (fast, simple operations)
- Counting/overlap queries are O(log n)

But we're seeing only 1.55x parallelism, suggesting:
- Decompression is NOT the bottleneck
- Something in record processing is slow
- Or there's contention/synchronization overhead

## Next Steps
1. Create minimal benchmark to isolate noodles BAM reading speed
2. Profile with `perf` or `flamegraph` to find real hotspots
3. Compare record-by-record processing time vs featureCounts
4. Consider if htslib bindings are necessary for competitive performance
