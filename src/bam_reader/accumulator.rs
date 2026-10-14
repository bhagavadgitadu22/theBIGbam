//! Contig accumulation independent of BAM fetching and worker scheduling.
use super::{prepare_cigar, BamPhaseTimings, MAX_REASONABLE_INSERT_SIZE, NEAR_CONTIG_END_DISTANCE};
use crate::features::{
    process_read, recycle_feature_arrays, take_feature_arrays, FeatureArrays, ModuleFlags,
};
use crate::types::SequencingType;
use rust_htslib::bam;

/// One contig's state. A reader starts it, feeds records, then finishes it
/// before advancing to the next reference. Scratch buffers belong to the reader.
pub(crate) struct ContigAccumulator {
    arrays: FeatureArrays,
    ref_length: usize,
    seq_type: SequencingType,
    flags: ModuleFlags,
    circular: bool,
    min_clipping_length: u32,
    has_reads: bool,
    primary_count: u64,
    process_read_ns: u64,
    time_on: bool,
}
impl ContigAccumulator {
    pub(crate) fn new(
        ref_length: usize,
        seq_type: SequencingType,
        flags: ModuleFlags,
        circular: bool,
        min_clipping_length: u32,
        time_on: bool,
    ) -> Self {
        Self {
            arrays: take_feature_arrays(ref_length, flags),
            ref_length,
            seq_type,
            flags,
            circular,
            min_clipping_length,
            has_reads: false,
            primary_count: 0,
            process_read_ns: 0,
            time_on,
        }
    }
    pub(crate) fn primary_count(&self) -> u64 { self.primary_count }

    pub(crate) fn consume_record(
        &mut self,
        record: &bam::Record,
        cigar_buf: &mut Vec<(u32, u32)>,
        seq_buf: &mut Vec<u8>,
    ) {
        let ref_length = self.ref_length;
        let seq_type = self.seq_type;
        let flags = self.flags;
        let circular = self.circular;
        let min_clipping_length = self.min_clipping_length;
        let time_on = self.time_on;
        let need_md = flags.needs_md();
        let need_seq = flags.mapping_metrics;
        let arrays = &mut self.arrays;
        if record.is_unmapped() {
            return;
        }

        self.has_reads = true;

        // Count primary alignments (not secondary, not supplementary)
        if !record.is_secondary() && !record.is_supplementary() {
            self.primary_count += 1;
        }

        // Parse CIGAR string and MD tag (if needed)
        let alignment_end = prepare_cigar(record.raw_cigar(), record.pos(), cigar_buf);

        let md_tag: Option<&[u8]> = if need_md {
            record.aux(b"MD").ok().and_then(|aux| match aux {
                rust_htslib::bam::record::Aux::String(s) => Some(s.as_bytes()),
                _ => None,
            })
        } else {
            None
        };

        // Extract SEQ bytes for sequence tracking (reuses buffer to avoid per-read allocation)
        if need_seq {
            seq_buf.clear();
            let seq = record.seq();
            let seq_len = seq.len();
            seq_buf.reserve(seq_len);
            for i in 0..seq_len {
                seq_buf.push(seq[i]);
            }
        }

        // Compute template length and proper_pair
        // BAM proper pair flag is unreliable (mapper-dependent), so recompute for all short-paired reads
        let (template_length, is_proper_pair, is_non_inward) = if seq_type.is_short_paired() {
            let pos1 = record.pos() as usize;
            let pos2 = record.mpos() as usize;

            let (corrected_tlen, corrected_proper) = if circular {
                // Circular: shortest path around the genome using actual (non-duplicated) length
                let p1 = pos1 % ref_length;
                let p2 = pos2 % ref_length;
                let direct = (p1 as i32 - p2 as i32).abs();
                let wrapped = ref_length as i32 - direct;
                let tlen = direct.min(wrapped);

                let same_ref = record.tid() == record.mtid();
                let opposite_strands = record.is_reverse() != record.is_mate_reverse();
                let reasonable = tlen > 0 && tlen < MAX_REASONABLE_INSERT_SIZE;
                (tlen, same_ref && opposite_strands && reasonable)
            } else {
                // Linear: use absolute insert size, recompute proper pair
                let tlen = record.insert_size().abs() as i32;
                let same_ref = record.tid() == record.mtid();
                let opposite_strands = record.is_reverse() != record.is_mate_reverse();
                let reasonable = tlen > 0 && tlen < MAX_REASONABLE_INSERT_SIZE;
                (tlen, same_ref && opposite_strands && reasonable)
            };

            let non_inward =
                record.tid() == record.mtid() && !record.is_mate_unmapped() && !corrected_proper;

            (corrected_tlen, corrected_proper, non_inward)
        } else {
            (
                record.insert_size().abs() as i32,
                record.is_proper_pair(),
                false,
            )
        };

        // Track circularising reads (primary alignments only, circular mode only)
        // Origin-crossing reads: alignment end exceeds contig length
        let _is_circularising_read =
            if circular && !record.is_secondary() && !record.is_supplementary() {
                let raw_start = record.pos() as usize;
                let raw_end = alignment_end as usize;

                if raw_end > ref_length {
                    arrays.circularising_reads_count += 1;
                    // Gate: confirm circularity if ≥20bp mapped on both sides
                    let left_overlap = ref_length - raw_start;
                    let right_overlap = raw_end - ref_length;
                    if left_overlap >= 20 && right_overlap >= 20 {
                        arrays.circularising_confirmed = true;
                    }
                    arrays
                        .circularising_min_overlaps
                        .push(left_overlap.min(right_overlap));
                    true
                } else {
                    false
                }
            } else {
                false
            };

        // Track circularising inserts & contig-end anomalies (primary paired-end only)
        // Skip reads already counted as circularising reads to avoid double-counting
        if seq_type.is_short_paired() && !record.is_secondary() && !record.is_supplementary() {
            // Detect reads near contig ends for circularisation analysis
            let pos = record.pos() as usize;
            let near_left = pos < NEAR_CONTIG_END_DISTANCE;
            let near_right = pos >= ref_length.saturating_sub(NEAR_CONTIG_END_DISTANCE);

            if near_left || near_right {
                // Circularising inserts: read1 near one end, mate near opposite end
                if record.is_first_in_template() {
                    if circular {
                        // In SAM-spec circular mode, both mates have POS < LN.
                        // Junction-spanning pairs have one mate near each end.
                        let mpos = record.mpos() as usize;
                        let mate_near_left = mpos < NEAR_CONTIG_END_DISTANCE;
                        let mate_near_right =
                            mpos >= ref_length.saturating_sub(NEAR_CONTIG_END_DISTANCE);
                        if (near_left && mate_near_right) || (near_right && mate_near_left) {
                            let bam_proper = record.is_proper_pair();
                            let bam_non_inward = !record.is_proper_pair()
                                && record.tid() == record.mtid()
                                && !record.is_mate_unmapped();
                            if bam_proper || bam_non_inward {
                                arrays.circularising_inserts_count += 1;
                                arrays.circularising_insert_sizes.push(template_length);
                            }
                        }
                    } else if is_non_inward {
                        // Linear mode: non-inward pairs with mate on opposite end
                        let mpos = record.mpos() as usize;
                        let mate_near_left = mpos < NEAR_CONTIG_END_DISTANCE;
                        let mate_near_right =
                            mpos >= ref_length.saturating_sub(NEAR_CONTIG_END_DISTANCE);
                        if (near_left && mate_near_right) || (near_right && mate_near_left) {
                            arrays.circularising_inserts_count += 1;
                            // Use wrapped distance: these pairs span the contig boundary,
                            // so the true insert size is ref_length - raw_tlen
                            let wrapped_tlen = ref_length as i32 - template_length;
                            arrays.circularising_insert_sizes.push(wrapped_tlen.abs());
                        }
                    }
                }

                // Unmapped mate with orientation suggesting mate beyond contig end
                // Forward on right end → mate expected further right → beyond right end
                // Reverse on left end  → mate expected further left  → beyond left end
                if record.is_mate_unmapped() {
                    if (near_right && !record.is_reverse()) || (near_left && record.is_reverse()) {
                        arrays.contig_end_unmapped_mates += 1;
                    }
                }

                // Mate on another contig near contig ends
                if record.tid() != record.mtid() && !record.is_mate_unmapped() {
                    arrays.contig_end_mates_mapped_on_another_contig += 1;
                }
            }

            // Collect all proper-pair insert sizes for overall mean/median baseline
            if record.is_first_in_template() && template_length > 0 && is_proper_pair {
                arrays.all_proper_insert_sizes.push(template_length);
            }
        }

        let t_pr = if time_on {
            Some(std::time::Instant::now())
        } else {
            None
        };
        process_read(
            arrays,
            record.pos(),
            alignment_end,
            record.seq_len() as i32,
            template_length,
            record.is_first_in_template(),
            is_proper_pair,
            record.is_reverse(),
            record.is_secondary(),
            record.is_supplementary(),
            record.is_mate_unmapped(),
            record.tid() != record.mtid(),
            cigar_buf,
            md_tag,
            if need_seq { seq_buf.as_slice() } else { &[] },
            record.mapq(),
            seq_type,
            flags,
            circular,
            min_clipping_length,
        );
        if let Some(t) = t_pr {
            self.process_read_ns += t.elapsed().as_nanos() as u64;
        }
    }
    /// Include reader decode/error overhead in records_ns, excluding process_read.
    pub(crate) fn record_loop_timing(&self, elapsed_ns: u64, timing: &mut Option<BamPhaseTimings>) {
        if let Some(pt) = timing.as_mut() {
            pt.records_ns = elapsed_ns.saturating_sub(self.process_read_ns);
            pt.process_read_ns = self.process_read_ns;
        }
    }
    pub(crate) fn finish(
        self,
        min_aligned_fraction: f64,
        phase_timings: &mut Option<BamPhaseTimings>,
    ) -> Option<(FeatureArrays, f64, u64)> {
        let Self {
            mut arrays,
            has_reads,
            primary_count,
            flags,
            seq_type,
            time_on,
            ..
        } = self;
        if !has_reads {
            recycle_feature_arrays(arrays);
            return None;
        }

        // Check coverage percentage before post-processing and feature calculation
        // This saves time on low-coverage contigs by skipping:
        // 1. finalize_strands() (phagetermini-specific, expensive)
        // 2. All feature compression in the caller (very expensive)
        // 3. Database writes
        let t_fin = if time_on {
            Some(std::time::Instant::now())
        } else {
            None
        };
        let coverage_pct = arrays.coverage_percentage();
        if coverage_pct < min_aligned_fraction {
            if let (Some(pt), Some(t)) = (phase_timings.as_mut(), t_fin) {
                pt.finalize_ns = t.elapsed().as_nanos() as u64;
            }
            recycle_feature_arrays(arrays);
            return None;
        }

        // Finalize strand-specific reads_starts/reads_ends for phagetermini
        // This merges forward/reverse strand tracking into final arrays
        // Only needed for phagetermini; other features are already finalized
        if flags.phagetermini {
            arrays.finalize_strands(seq_type);
        }
        if let (Some(pt), Some(t)) = (phase_timings.as_mut(), t_fin) {
            pt.finalize_ns = t.elapsed().as_nanos() as u64;
        }

        Some((arrays, coverage_pct, primary_count))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::header::HeaderRecord;

    fn record(flag: u16) -> bam::Record {
        let mut header = bam::Header::new();
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "a")
                .push_tag(b"LN", 100),
        );
        let view = bam::HeaderView::from_header(&header);
        let line = format!("read\t{flag}\ta\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\t*\tMD:Z:10");
        bam::Record::from_sam(&view, line.as_bytes()).unwrap()
    }

    #[test]
    fn finish_preserves_primary_counts_coverage_and_filter_boundary() {
        let flags = ModuleFlags::from_modules(&["Coverage".to_owned()]);
        for threshold in [0.0, 10.0, 10.1] {
            let mut accumulator =
                ContigAccumulator::new(100, SequencingType::ShortSingle, flags, false, 5, true);
            let (mut cigar, mut sequence) = (Vec::new(), Vec::new());
            for flag in [0, 256, 2048] {
                accumulator.consume_record(&record(flag), &mut cigar, &mut sequence);
            }
            let mut timing = Some(BamPhaseTimings::default());
            accumulator.record_loop_timing(u64::MAX, &mut timing);
            assert_eq!(
                timing.as_ref().unwrap().records_ns + timing.as_ref().unwrap().process_read_ns,
                u64::MAX
            );
            let result = accumulator.finish(threshold, &mut timing);
            if threshold > 10.0 {
                assert!(result.is_none());
            } else {
                let (arrays, coverage, primary) = result.unwrap();
                assert_eq!(primary, 1);
                assert_eq!(coverage, 10.0);
                assert_eq!(&arrays.primary_reads[..10], &[1; 10]);
                assert!(arrays.primary_reads[10..].iter().all(|&x| x == 0));
                recycle_feature_arrays(arrays);
            }
        }
    }

    #[test]
    fn empty_contig_is_absent_even_without_coverage_filter() {
        let flags = ModuleFlags::from_modules(&["Coverage".to_owned()]);
        let accumulator =
            ContigAccumulator::new(100, SequencingType::ShortSingle, flags, false, 5, false);
        assert!(accumulator.finish(0.0, &mut None).is_none());
    }
}
