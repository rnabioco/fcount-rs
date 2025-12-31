/// QC statistics for read counting
#[derive(Debug, Default, Clone)]
pub struct ReadCounters {
    /// Successfully assigned reads/fragments
    pub assigned: u64,
    /// Unmapped reads (FLAG 0x4)
    pub unassigned_unmapped: u64,
    /// Only one mate mapped (paired-end with -B)
    pub unassigned_singleton: u64,
    /// Below minimum mapping quality
    pub unassigned_mapping_quality: u64,
    /// Mates on different chromosomes
    pub unassigned_chimeric: u64,
    /// Fragment length outside range
    pub unassigned_fragment_length: u64,
    /// Duplicate reads (FLAG 0x400 with --ignore-dup)
    pub unassigned_duplicate: u64,
    /// Multi-mapping without -M
    pub unassigned_multimapping: u64,
    /// Secondary/supplementary alignments
    pub unassigned_secondary: u64,
    /// No overlapping features
    pub unassigned_no_features: u64,
    /// Below minimum overlap threshold
    pub unassigned_overlap_length: u64,
    /// Multiple features without -O
    pub unassigned_ambiguous: u64,
}

impl ReadCounters {
    /// Total unassigned reads
    pub fn total_unassigned(&self) -> u64 {
        self.unassigned_unmapped
            + self.unassigned_singleton
            + self.unassigned_mapping_quality
            + self.unassigned_chimeric
            + self.unassigned_fragment_length
            + self.unassigned_duplicate
            + self.unassigned_multimapping
            + self.unassigned_secondary
            + self.unassigned_no_features
            + self.unassigned_overlap_length
            + self.unassigned_ambiguous
    }

    /// Total reads processed
    pub fn total(&self) -> u64 {
        self.assigned + self.total_unassigned()
    }

    /// Assignment rate as a fraction
    pub fn assignment_rate(&self) -> f64 {
        let total = self.total();
        if total == 0 {
            0.0
        } else {
            self.assigned as f64 / total as f64
        }
    }

    /// Merge another counter into this one
    pub fn merge(&mut self, other: &ReadCounters) {
        self.assigned += other.assigned;
        self.unassigned_unmapped += other.unassigned_unmapped;
        self.unassigned_singleton += other.unassigned_singleton;
        self.unassigned_mapping_quality += other.unassigned_mapping_quality;
        self.unassigned_chimeric += other.unassigned_chimeric;
        self.unassigned_fragment_length += other.unassigned_fragment_length;
        self.unassigned_duplicate += other.unassigned_duplicate;
        self.unassigned_multimapping += other.unassigned_multimapping;
        self.unassigned_secondary += other.unassigned_secondary;
        self.unassigned_no_features += other.unassigned_no_features;
        self.unassigned_overlap_length += other.unassigned_overlap_length;
        self.unassigned_ambiguous += other.unassigned_ambiguous;
    }
}

impl std::fmt::Display for ReadCounters {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Assigned:                    {}", self.assigned)?;
        writeln!(f, "Unassigned_Unmapped:         {}", self.unassigned_unmapped)?;
        writeln!(f, "Unassigned_Singleton:        {}", self.unassigned_singleton)?;
        writeln!(f, "Unassigned_MappingQuality:   {}", self.unassigned_mapping_quality)?;
        writeln!(f, "Unassigned_Chimera:          {}", self.unassigned_chimeric)?;
        writeln!(f, "Unassigned_FragmentLength:   {}", self.unassigned_fragment_length)?;
        writeln!(f, "Unassigned_Duplicate:        {}", self.unassigned_duplicate)?;
        writeln!(f, "Unassigned_MultiMapping:     {}", self.unassigned_multimapping)?;
        writeln!(f, "Unassigned_Secondary:        {}", self.unassigned_secondary)?;
        writeln!(f, "Unassigned_NoFeatures:       {}", self.unassigned_no_features)?;
        writeln!(f, "Unassigned_OverlapLength:    {}", self.unassigned_overlap_length)?;
        writeln!(f, "Unassigned_Ambiguity:        {}", self.unassigned_ambiguous)?;
        Ok(())
    }
}
