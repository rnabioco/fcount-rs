pub mod block_reader;
mod cigar;
pub mod minimal_parser;

pub use cigar::{Interval, total_overlap};
