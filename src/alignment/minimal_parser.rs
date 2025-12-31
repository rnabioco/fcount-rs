//! Minimal BAM record parser for high-performance feature counting.
//!
//! This module provides zero-allocation parsing of BAM records, extracting
//! only the fields needed for feature counting. This mimics the approach
//! used by featureCounts for maximum performance.

use super::cigar::Interval;
use smallvec::SmallVec;

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
pub const FLAG_PAIRED: u16 = 0x1;
pub const FLAG_PROPER_PAIR: u16 = 0x2;
pub const FLAG_UNMAPPED: u16 = 0x4;
pub const FLAG_MATE_UNMAPPED: u16 = 0x8;
pub const FLAG_REVERSE: u16 = 0x10;
pub const FLAG_MATE_REVERSE: u16 = 0x20;
pub const FLAG_FIRST: u16 = 0x40;
pub const FLAG_SECOND: u16 = 0x80;
pub const FLAG_SECONDARY: u16 = 0x100;
pub const FLAG_FILTERED: u16 = 0x200;
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
    pub intervals: SmallVec<[Interval; 4]>,
    /// NH tag value (number of alignments, 1 if not present)
    pub nh: u8,
    /// Mate reference ID (for paired-end)
    pub mate_ref_id: i32,
    /// Mate position (for paired-end)
    pub mate_pos: i32,
    /// Template length (for paired-end)
    pub tlen: i32,
    /// Read name bytes (for paired-end mate matching)
    pub read_name: SmallVec<[u8; 64]>,
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
            mate_pos: 0,
            tlen: 0,
            read_name: SmallVec::new(),
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
        self.mate_pos = 0;
        self.tlen = 0;
        self.read_name.clear();
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
    pub fn is_mate_reverse(&self) -> bool {
        self.flags & FLAG_MATE_REVERSE != 0
    }

    #[inline]
    pub fn is_first_in_pair(&self) -> bool {
        self.flags & FLAG_FIRST != 0
    }

    #[inline]
    pub fn is_second_in_pair(&self) -> bool {
        self.flags & FLAG_SECOND != 0
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

    #[inline]
    pub fn is_paired(&self) -> bool {
        self.flags & FLAG_PAIRED != 0
    }

    #[inline]
    pub fn is_proper_pair(&self) -> bool {
        self.flags & FLAG_PROPER_PAIR != 0
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
) -> Result<(), &'static str> {
    // Minimum BAM record size: 32 bytes fixed header + 1 byte read name (null)
    if data.len() < 33 {
        return Err("Record too short");
    }

    // Parse fixed-size fields (BAM spec: all little-endian)
    // Offsets within record data (after block_size):
    // 0-3: refID (i32)
    // 4-7: pos (i32)
    // 8: l_read_name (u8)
    // 9: mapq (u8)
    // 10-11: bin (u16) - ignored
    // 12-13: n_cigar_op (u16)
    // 14-15: flag (u16)
    // 16-19: l_seq (i32)
    // 20-23: next_refID (i32)
    // 24-27: next_pos (i32)
    // 28-31: tlen (i32)

    record.clear();

    record.ref_id = i32::from_le_bytes([data[0], data[1], data[2], data[3]]);
    record.pos = i32::from_le_bytes([data[4], data[5], data[6], data[7]]);
    let l_read_name = data[8] as usize;
    record.mapq = data[9];
    let n_cigar_op = u16::from_le_bytes([data[12], data[13]]) as usize;
    record.flags = u16::from_le_bytes([data[14], data[15]]);
    let l_seq = i32::from_le_bytes([data[16], data[17], data[18], data[19]]) as usize;
    record.mate_ref_id = i32::from_le_bytes([data[20], data[21], data[22], data[23]]);
    record.mate_pos = i32::from_le_bytes([data[24], data[25], data[26], data[27]]);
    record.tlen = i32::from_le_bytes([data[28], data[29], data[30], data[31]]);

    // Variable data starts at offset 32
    let var_start = 32;

    // Parse read name if needed (for paired-end mate matching)
    if need_read_name && l_read_name > 0 {
        let name_end = var_start + l_read_name - 1; // Exclude null terminator
        if name_end <= data.len() {
            record
                .read_name
                .extend_from_slice(&data[var_start..name_end]);
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

    // Aux data starts after seq and qual
    let seq_len = (l_seq + 1) / 2;
    let aux_start = cigar_end + seq_len + l_seq;

    // Parse NH tag from aux data
    if aux_start < data.len() {
        record.nh = find_nh_tag(&data[aux_start..]).unwrap_or(1);
    }

    Ok(())
}

/// Parse raw CIGAR ops into genomic intervals
#[inline]
fn parse_cigar_ops(cigar_data: &[u8], start_pos: i32, out: &mut SmallVec<[Interval; 4]>) {
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
        let op_len = (op_raw >> 4) as u32;

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
                if let Some(interval) = current_interval.take() {
                    if interval.end >= interval.start {
                        out.push(interval);
                    }
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
    if let Some(interval) = current_interval {
        if interval.end >= interval.start {
            out.push(interval);
        }
    }
}

/// Find NH tag in auxiliary data
///
/// Aux data format: TAG (2 bytes) + TYPE (1 byte) + VALUE
/// NH tag is 'N' 'H' with integer type
#[inline]
fn find_nh_tag(aux_data: &[u8]) -> Option<u8> {
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
                // string or hex
                let start = offset;
                while offset < aux_data.len() && aux_data[offset] != 0 {
                    offset += 1;
                }
                offset - start + 1 // Include null terminator
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

/// Parse integer value from aux tag
#[inline]
fn parse_int_tag(data: &[u8], val_type: u8) -> Option<u8> {
    match val_type {
        b'c' if !data.is_empty() => Some(data[0] as i8 as u8),
        b'C' if !data.is_empty() => Some(data[0]),
        b's' if data.len() >= 2 => Some(i16::from_le_bytes([data[0], data[1]]) as u8),
        b'S' if data.len() >= 2 => Some(u16::from_le_bytes([data[0], data[1]]) as u8),
        b'i' if data.len() >= 4 => {
            Some(i32::from_le_bytes([data[0], data[1], data[2], data[3]]) as u8)
        }
        b'I' if data.len() >= 4 => {
            Some(u32::from_le_bytes([data[0], data[1], data[2], data[3]]) as u8)
        }
        _ => None,
    }
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

    #[test]
    fn test_cigar_with_intron() {
        // Test 50M100N50M
        let mut cigar_data = Vec::new();
        cigar_data.extend_from_slice(&((50u32 << 4) | 0).to_le_bytes()); // 50M
        cigar_data.extend_from_slice(&((100u32 << 4) | 3).to_le_bytes()); // 100N
        cigar_data.extend_from_slice(&((50u32 << 4) | 0).to_le_bytes()); // 50M

        let mut intervals = SmallVec::new();
        parse_cigar_ops(&cigar_data, 999, &mut intervals);

        assert_eq!(intervals.len(), 2);
        assert_eq!(intervals[0].start, 1000);
        assert_eq!(intervals[0].end, 1049);
        assert_eq!(intervals[1].start, 1150);
        assert_eq!(intervals[1].end, 1199);
    }
}
