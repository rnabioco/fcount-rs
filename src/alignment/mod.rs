pub mod block_reader;
mod cigar;
pub mod minimal_parser;
mod paired;
mod reader;

pub use cigar::{parse_cigar_intervals, total_overlap, Interval};
pub use paired::{MateTracker, PendingMate};
pub use reader::{AlignmentReader, AlignmentRecord};
