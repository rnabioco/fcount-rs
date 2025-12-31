# Usage

## Basic Command

```bash
fcount-rs -a <annotation.gtf> -o <output.txt> <bam_files>...
```

## Options

### Required

| Option | Description |
|--------|-------------|
| `-a, --annotation <PATH>` | GTF/GFF annotation file |
| `-o, --output <PATH>` | Output count matrix file |
| `<BAM_FILES>` | One or more BAM/SAM input files |

### Counting

| Option | Description | Default |
|--------|-------------|---------|
| `-t, --type <TYPE>` | Feature type to count | `exon` |
| `-g, --gene-id <ATTR>` | GTF attribute for gene ID | `gene_id` |
| `-f, --feature-level` | Count at feature level instead of gene level | off |

### Paired-End

| Option | Description | Default |
|--------|-------------|---------|
| `-p, --paired` | Count fragments instead of reads | off |
| `-B, --both-aligned` | Require both mates aligned | off |
| `-C, --count-chimeric` | Count chimeric fragments | off |
| `--min-frag-len <INT>` | Minimum fragment length | 50 |
| `--max-frag-len <INT>` | Maximum fragment length | 600 |

### Strandedness

| Option | Description | Default |
|--------|-------------|---------|
| `-s, --strand <MODE>` | Strand mode: 0=unstranded, 1=stranded, 2=reversely stranded | 0 |

### Multi-Mapping

| Option | Description | Default |
|--------|-------------|---------|
| `-M, --multi-mapping` | Count multi-mapping reads | off |
| `--fraction` | Use fractional counting | off |
| `--primary` | Count primary alignments only | off |

### Overlap

| Option | Description | Default |
|--------|-------------|---------|
| `-O, --multi-overlap` | Allow reads assigned to multiple features | off |
| `--min-overlap <INT>` | Minimum overlapping bases | 1 |
| `--frac-overlap <FLOAT>` | Minimum fraction of read overlapping feature | 0.0 |
| `--frac-overlap-feature <FLOAT>` | Minimum fraction of feature overlapped | 0.0 |
| `--largest-overlap` | Assign to feature with largest overlap only | off |

### Filtering

| Option | Description | Default |
|--------|-------------|---------|
| `-Q, --min-mapq <INT>` | Minimum mapping quality | 0 |
| `--ignore-dup` | Ignore duplicate reads | off |

### Performance

| Option | Description | Default |
|--------|-------------|---------|
| `-T, --threads <INT>` | Number of threads | 1 |

### Output

| Option | Description |
|--------|-------------|
| `-R, --details <PATH>` | Output detailed read assignments |
| `-q, --quiet` | Suppress progress output |
