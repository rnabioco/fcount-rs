# Examples

## Single-End RNA-seq

```bash
fcount-rs -a annotation.gtf -o counts.txt sample.bam
```

## Paired-End RNA-seq

```bash
fcount-rs -p -a annotation.gtf -o counts.txt sample.bam
```

## Stranded Library (e.g., dUTP method)

```bash
fcount-rs -p -s 2 -a annotation.gtf -o counts.txt sample.bam
```

## Multi-Threaded Processing

```bash
fcount-rs -T 8 -a annotation.gtf -o counts.txt sample.bam
```

## Multiple BAM Files

```bash
fcount-rs -a annotation.gtf -o counts.txt sample1.bam sample2.bam sample3.bam
```

## Count Multi-Mapping Reads

```bash
fcount-rs -M -a annotation.gtf -o counts.txt sample.bam
```

With fractional counting:

```bash
fcount-rs -M --fraction -a annotation.gtf -o counts.txt sample.bam
```

## Feature-Level Counting

Count at exon level instead of gene level:

```bash
fcount-rs -f -a annotation.gtf -o counts.txt sample.bam
```

## Quality Filtering

Require minimum MAPQ of 10 and ignore duplicates:

```bash
fcount-rs -Q 10 --ignore-dup -a annotation.gtf -o counts.txt sample.bam
```
