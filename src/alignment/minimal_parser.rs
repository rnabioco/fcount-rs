//! Minimal BAM record parser for high-performance feature counting.
//!
//! This module provides zero-allocation parsing of BAM records, extracting
//! only the fields needed for feature counting. This mimics the approach
//! used by featureCounts for maximum performance.

use super::cigar::Interval;
use rustc_hash::FxHasher;
use smallvec::SmallVec;
use std::hash::{Hash, Hasher};

/// CIGAR operation types (from BAM specification)
const CIGAR_M: u8 = 0; // Match or mismatch
const CIGAR_I: u8 = 1; // Insertion
const CIGAR_D: u8 = 2; // Deletion
const CIGAR_N: u8 = 3; // Skip (intron)
const CIGAR_S: u8 = 4; // Soft clip
const CIGAR_H: u8 = 5; // Hard clip
const CIGAR_P: u8 = 6; // Padding
const CIGAR_EQ: u8 = 7; // Sequence match
const CIGAR_X: u8 = 8; // Sequence mismatch

/// SAM flags
pub const FLAG_UNMAPPED: u16 = 0x4;
pub const FLAG_MATE_UNMAPPED: u16 = 0x8;
pub const FLAG_REVERSE: u16 = 0x10;
pub const FLAG_SECONDARY: u16 = 0x100;
pub const FLAG_DUPLICATE: u16 = 0x400;
pub const FLAG_SUPPLEMENTARY: u16 = 0x800;

/// Minimal BAM record with only fields needed for feature counting
#[derive(Debug, Clone)]
pub struct MinimalRecord {
    /// SAM flags
    pub flags: u16,
    /// Reference sequence ID (-1 if unmapped)
    pub ref_id: i32,
    /// 0-based leftmost position
    pub pos: i32,
    /// Mapping quality
    pub mapq: u8,
    /// Genomic intervals from CIGAR
    pub intervals: SmallVec<[Interval; 8]>,
    /// NH tag value (number of alignments, 1 if not present)
    pub nh: u16,
    /// Mate reference ID (for paired-end chimeric check)
    pub mate_ref_id: i32,
    /// FxHash of the read name (for paired-end mate matching).
    /// Set only when `need_read_name` is true during parse; 0 otherwise.
    /// We hash directly from the BAM bytes to avoid copying the name into
    /// a SmallVec just to hash it.
    pub read_name_hash: u64,
    // Note: mate_pos and tlen are present in the BAM header but never
    // consumed by counting logic — we don't parse them.
}

impl MinimalRecord {
    /// Create a new empty record (for reuse)
    pub fn new() -> Self {
        MinimalRecord {
            flags: 0,
            ref_id: -1,
            pos: 0,
            mapq: 0,
            intervals: SmallVec::new(),
            nh: 1,
            mate_ref_id: -1,
            read_name_hash: 0,
        }
    }

    /// Clear for reuse
    pub fn clear(&mut self) {
        self.flags = 0;
        self.ref_id = -1;
        self.pos = 0;
        self.mapq = 0;
        self.intervals.clear();
        self.nh = 1;
        self.mate_ref_id = -1;
        self.read_name_hash = 0;
    }

    #[inline]
    pub fn is_unmapped(&self) -> bool {
        self.flags & FLAG_UNMAPPED != 0
    }

    #[inline]
    pub fn is_mate_unmapped(&self) -> bool {
        self.flags & FLAG_MATE_UNMAPPED != 0
    }

    #[inline]
    pub fn is_reverse(&self) -> bool {
        self.flags & FLAG_REVERSE != 0
    }

    #[inline]
    pub fn is_secondary(&self) -> bool {
        self.flags & FLAG_SECONDARY != 0
    }

    #[inline]
    pub fn is_duplicate(&self) -> bool {
        self.flags & FLAG_DUPLICATE != 0
    }

    #[inline]
    pub fn is_supplementary(&self) -> bool {
        self.flags & FLAG_SUPPLEMENTARY != 0
    }
}

impl Default for MinimalRecord {
    fn default() -> Self {
        Self::new()
    }
}

/// Parse a BAM record from raw bytes into a MinimalRecord
///
/// # Arguments
/// * `data` - Raw BAM record bytes (NOT including the 4-byte block_size prefix)
/// * `record` - Output record (will be cleared and filled)
/// * `need_read_name` - Whether to parse read name (for paired-end)
///
/// # Returns
/// * `Ok(())` if parsing succeeded
/// * `Err(msg)` if parsing failed
#[inline]
pub fn parse_bam_record(
    data: &[u8],
    record: &mut MinimalRecord,
    need_read_name: bool,
    need_nh_tag: bool,
) -> Result<(), &'static str> {
    // Minimum BAM record size: 32 bytes fixed header + 1 byte read name (null)
    if data.len() < 33 {
        return Err("Record too short");
    }

    // Parse fixed-size fields (BAM spec: all little-endian). Offsets within
    // record data (after block_size):
    //   0-3 refID (i32)     | 4-7 pos (i32)     | 8 l_read_name (u8) | 9 mapq (u8)
    //   10-11 bin (ignored) | 12-13 n_cigar_op (u16) | 14-15 flag (u16)
    //   16-19 l_seq (i32)   | 20-23 next_refID (i32)
    //   24-31 next_pos + tlen — intentionally skipped (see struct comment)
    //
    // Read the header as a single 32-byte chunk so the compiler sees one
    // bounds check feeding all field extractions, instead of 10 indexed reads.

    record.clear();

    let hdr: &[u8; 32] = data[..32]
        .try_into()
        .expect("data.len() ≥ 33 checked above");

    record.ref_id = i32::from_le_bytes(hdr[0..4].try_into().unwrap());
    record.pos = i32::from_le_bytes(hdr[4..8].try_into().unwrap());
    let l_read_name = hdr[8] as usize;
    record.mapq = hdr[9];
    let n_cigar_op = u16::from_le_bytes(hdr[12..14].try_into().unwrap()) as usize;
    record.flags = u16::from_le_bytes(hdr[14..16].try_into().unwrap());
    let l_seq = i32::from_le_bytes(hdr[16..20].try_into().unwrap()) as usize;
    record.mate_ref_id = i32::from_le_bytes(hdr[20..24].try_into().unwrap());

    // Variable data starts at offset 32
    let var_start = 32;

    // Hash the read name directly from BAM bytes if needed (for paired-end
    // mate matching). Avoids copying ~50-100 bytes per record into a SmallVec
    // when all we ever do with the name is hash it.
    if need_read_name && l_read_name > 0 {
        let name_end = var_start + l_read_name - 1; // Exclude null terminator
        if name_end <= data.len() {
            let mut hasher = FxHasher::default();
            data[var_start..name_end].hash(&mut hasher);
            record.read_name_hash = hasher.finish();
        }
    }

    // CIGAR starts after read name
    let cigar_start = var_start + l_read_name;
    let cigar_end = cigar_start + n_cigar_op * 4;

    if cigar_end > data.len() {
        return Err("CIGAR extends past record");
    }

    // Parse CIGAR into intervals
    parse_cigar_ops(
        &data[cigar_start..cigar_end],
        record.pos,
        &mut record.intervals,
    );

    // Parse NH tag from aux data only if needed (for multi-mapper filtering)
    if need_nh_tag {
        let seq_len = l_seq.div_ceil(2);
        let aux_start = cigar_end + seq_len + l_seq;
        if aux_start < data.len() {
            record.nh = find_nh_tag(&data[aux_start..]).unwrap_or(1);
        }
    }
    // Otherwise nh stays at default value of 1 (set in clear())

    Ok(())
}

/// Parse raw CIGAR ops into genomic intervals
#[inline]
fn parse_cigar_ops(cigar_data: &[u8], start_pos: i32, out: &mut SmallVec<[Interval; 8]>) {
    out.clear();

    if start_pos < 0 {
        return;
    }

    let mut ref_pos = (start_pos + 1) as u32; // Convert to 1-based
    let mut current_interval: Option<Interval> = None;

    // Each CIGAR op is 4 bytes (u32 little-endian)
    // Lower 4 bits: op type, upper 28 bits: length
    for chunk in cigar_data.chunks_exact(4) {
        let op_raw = u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
        let op_type = (op_raw & 0xF) as u8;
        let op_len = op_raw >> 4;

        match op_type {
            // Operations that consume reference and are part of alignment
            CIGAR_M | CIGAR_EQ | CIGAR_X => {
                let interval_end = ref_pos + op_len - 1;
                match current_interval.as_mut() {
                    Some(interval) => {
                        interval.end = interval_end;
                    }
                    None => {
                        current_interval = Some(Interval {
                            start: ref_pos,
                            end: interval_end,
                        });
                    }
                }
                ref_pos += op_len;
            }

            // Deletion: consumes reference, extends interval
            CIGAR_D => {
                let interval_end = ref_pos + op_len - 1;
                match current_interval.as_mut() {
                    Some(interval) => {
                        interval.end = interval_end;
                    }
                    None => {
                        current_interval = Some(Interval {
                            start: ref_pos,
                            end: interval_end,
                        });
                    }
                }
                ref_pos += op_len;
            }

            // Skip (intron): ends current interval
            CIGAR_N => {
                if let Some(interval) = current_interval.take()
                    && interval.end >= interval.start
                {
                    out.push(interval);
                }
                ref_pos += op_len;
            }

            // Insertion, soft clip: don't advance reference
            CIGAR_I | CIGAR_S => {}

            // Hard clip, padding: no effect
            CIGAR_H | CIGAR_P => {}

            _ => {}
        }
    }

    // Save final interval
    if let Some(interval) = current_interval
        && interval.end >= interval.start
    {
        out.push(interval);
    }
}

/// Find NH tag in auxiliary data
///
/// Aux data format: TAG (2 bytes) + TYPE (1 byte) + VALUE
/// NH tag is 'N' 'H' with integer type
#[inline]
fn find_nh_tag(aux_data: &[u8]) -> Option<u16> {
    let mut offset = 0;

    while offset + 3 <= aux_data.len() {
        let tag1 = aux_data[offset];
        let tag2 = aux_data[offset + 1];
        let val_type = aux_data[offset + 2];

        // Check for NH tag
        if tag1 == b'N' && tag2 == b'H' {
            return parse_int_tag(&aux_data[offset + 3..], val_type);
        }

        // Skip to next tag based on type
        offset += 3;
        offset += match val_type {
            b'A' => 1,        // char
            b'c' | b'C' => 1, // int8
            b's' | b'S' => 2, // int16
            b'i' | b'I' => 4, // int32
            b'f' => 4,        // float
            b'Z' | b'H' => {
                // NUL-terminated string or hex. Scan with a *local* cursor: this
                // arm's value is added to `offset` by the caller, so advancing
                // `offset` in here too would skip the value twice and land the
                // walk mid-tag. That silently lost every NH tag written after a
                // Z-type tag (e.g. HISAT2's MD:Z, or an added RG:Z).
                let mut end = offset;
                while end < aux_data.len() && aux_data[end] != 0 {
                    end += 1;
                }
                end - offset + 1 // Include null terminator
            }
            b'B' => {
                // array
                if offset + 5 > aux_data.len() {
                    return None;
                }
                let arr_type = aux_data[offset];
                let arr_len = u32::from_le_bytes([
                    aux_data[offset + 1],
                    aux_data[offset + 2],
                    aux_data[offset + 3],
                    aux_data[offset + 4],
                ]) as usize;
                let elem_size = match arr_type {
                    b'c' | b'C' => 1,
                    b's' | b'S' => 2,
                    b'i' | b'I' | b'f' => 4,
                    _ => return None,
                };
                5 + arr_len * elem_size
            }
            _ => return None, // Unknown type
        };
    }

    None
}

/// Parse integer value from aux tag.
///
/// Widened through `i64` and clamped rather than cast straight to the return
/// type: NH is a count, and truncating it corrupts results silently. `NH:i:300`
/// used to wrap to 44 (giving the read ~7× its correct fractional weight), and
/// `NH:i:256` wrapped to 0, which then divided by zero in `calculate_count`.
/// Clamping the low end at 1 also absorbs a malformed `NH:i:0`.
#[inline]
fn parse_int_tag(data: &[u8], val_type: u8) -> Option<u16> {
    let raw: i64 = match val_type {
        b'c' if !data.is_empty() => data[0] as i8 as i64,
        b'C' if !data.is_empty() => data[0] as i64,
        b's' if data.len() >= 2 => i16::from_le_bytes([data[0], data[1]]) as i64,
        b'S' if data.len() >= 2 => u16::from_le_bytes([data[0], data[1]]) as i64,
        b'i' if data.len() >= 4 => i32::from_le_bytes([data[0], data[1], data[2], data[3]]) as i64,
        b'I' if data.len() >= 4 => u32::from_le_bytes([data[0], data[1], data[2], data[3]]) as i64,
        _ => return None,
    };
    Some(raw.clamp(1, u16::MAX as i64) as u16)
}

/// Get the size of a BAM record from its first 4 bytes
#[inline]
pub fn get_record_size(block_size_bytes: &[u8]) -> usize {
    if block_size_bytes.len() < 4 {
        return 0;
    }
    u32::from_le_bytes([
        block_size_bytes[0],
        block_size_bytes[1],
        block_size_bytes[2],
        block_size_bytes[3],
    ]) as usize
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build an aux block from `(tag, type, value_bytes)` triples.
    fn aux(tags: &[(&[u8; 2], u8, Vec<u8>)]) -> Vec<u8> {
        let mut out = Vec::new();
        for (tag, ty, val) in tags {
            out.extend_from_slice(*tag);
            out.push(*ty);
            out.extend_from_slice(val);
        }
        out
    }

    fn nh_i(n: i32) -> (&'static [u8; 2], u8, Vec<u8>) {
        (b"NH", b'i', n.to_le_bytes().to_vec())
    }

    fn z_tag(tag: &'static [u8; 2], s: &str) -> (&'static [u8; 2], u8, Vec<u8>) {
        let mut v = s.as_bytes().to_vec();
        v.push(0); // NUL terminator
        (tag, b'Z', v)
    }

    #[test]
    fn test_nh_tag_alone() {
        assert_eq!(find_nh_tag(&aux(&[nh_i(5)])), Some(5));
    }

    /// Regression: the Z/H arm used to advance the cursor *and* return the
    /// length, so the walk double-skipped and every NH after a string tag was
    /// lost — silently disabling multi-mapper filtering on HISAT2 output
    /// (MD:Z precedes NH) and on anything with a read group added.
    #[test]
    fn test_nh_tag_after_string_tags() {
        // MD:Z before NH — the HISAT2 layout.
        assert_eq!(find_nh_tag(&aux(&[z_tag(b"MD", "100"), nh_i(5)])), Some(5));
        // A shorter Z value, to catch an off-by-one that only shows at one length.
        assert_eq!(find_nh_tag(&aux(&[z_tag(b"RG", "s1"), nh_i(3)])), Some(3));
        // Empty Z value: just the terminator.
        assert_eq!(find_nh_tag(&aux(&[z_tag(b"RG", ""), nh_i(7)])), Some(7));
        // Several string tags in a row.
        assert_eq!(
            find_nh_tag(&aux(&[
                z_tag(b"MD", "76"),
                z_tag(b"RG", "sample_a"),
                z_tag(b"SA", "chr1,100,+,50M,60,0;"),
                nh_i(11),
            ])),
            Some(11)
        );
    }

    /// The other aux types must still be skipped correctly.
    #[test]
    fn test_nh_tag_after_mixed_tags() {
        let block = aux(&[
            (b"AS", b'i', 42i32.to_le_bytes().to_vec()),
            (b"XN", b'c', vec![3]),
            (b"XS", b's', 7i16.to_le_bytes().to_vec()),
            (b"XF", b'f', 1.5f32.to_le_bytes().to_vec()),
            (b"XA", b'A', vec![b'+']),
            // B array: 3 x int32
            (b"ZB", b'B', {
                let mut v = vec![b'i'];
                v.extend_from_slice(&3u32.to_le_bytes());
                v.extend_from_slice(&1i32.to_le_bytes());
                v.extend_from_slice(&2i32.to_le_bytes());
                v.extend_from_slice(&3i32.to_le_bytes());
                v
            }),
            z_tag(b"MD", "76"),
            nh_i(9),
        ]);
        assert_eq!(find_nh_tag(&block), Some(9));
    }

    #[test]
    fn test_nh_tag_absent() {
        assert_eq!(find_nh_tag(&aux(&[z_tag(b"MD", "100")])), None);
        assert_eq!(find_nh_tag(&[]), None);
    }

    /// Regression: NH used to be parsed into a `u8`. NH:i:256 wrapped to 0 and
    /// then divided by zero in `calculate_count`; NH:i:300 wrapped to 44.
    #[test]
    fn test_nh_tag_large_values_do_not_wrap() {
        assert_eq!(find_nh_tag(&aux(&[nh_i(255)])), Some(255));
        assert_eq!(find_nh_tag(&aux(&[nh_i(256)])), Some(256));
        assert_eq!(find_nh_tag(&aux(&[nh_i(300)])), Some(300));
        assert_eq!(find_nh_tag(&aux(&[nh_i(70_000)])), Some(u16::MAX));
        // Malformed non-positive NH clamps to 1 rather than producing a zero
        // divisor downstream.
        assert_eq!(find_nh_tag(&aux(&[nh_i(0)])), Some(1));
        assert_eq!(find_nh_tag(&aux(&[nh_i(-4)])), Some(1));
    }

    #[test]
    fn test_cigar_parsing() {
        // Test simple match: 100M
        // CIGAR op: (100 << 4) | 0 = 1600
        let cigar_data: [u8; 4] = 1600u32.to_le_bytes();
        let mut intervals = SmallVec::new();
        parse_cigar_ops(&cigar_data, 999, &mut intervals); // 0-based pos 999 -> 1-based 1000

        assert_eq!(intervals.len(), 1);
        assert_eq!(intervals[0].start, 1000);
        assert_eq!(intervals[0].end, 1099);
    }

    /// Encode one CIGAR op the way BAM stores it: length in the upper 28 bits,
    /// opcode in the low 4.
    fn cigar_op(len: u32, op: u8) -> [u8; 4] {
        ((len << 4) | op as u32).to_le_bytes()
    }

    #[test]
    fn test_cigar_with_intron() {
        // Test 50M100N50M
        let mut cigar_data = Vec::new();
        cigar_data.extend_from_slice(&cigar_op(50, CIGAR_M));
        cigar_data.extend_from_slice(&cigar_op(100, CIGAR_N));
        cigar_data.extend_from_slice(&cigar_op(50, CIGAR_M));

        let mut intervals = SmallVec::new();
        parse_cigar_ops(&cigar_data, 999, &mut intervals);

        assert_eq!(intervals.len(), 2);
        assert_eq!(intervals[0].start, 1000);
        assert_eq!(intervals[0].end, 1049);
        assert_eq!(intervals[1].start, 1150);
        assert_eq!(intervals[1].end, 1199);
    }
}
