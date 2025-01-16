// Workflow for the preparation of wgs samples
// Import required processes
include { wes_trimming               }   from '../modules/wes_trimming.nf'
include { wes_trimming_s             }   from '../modules/wes_trimming_s.nf'
include { alignment                  }   from '../modules/alignment.nf'
include { alignment_s                }   from '../modules/alignment_s.nf'
include { wes_bam_prep               }   from '../modules/wes_bam_prep.nf'
include { picard_rrg                 }   from '../modules/picard_rrg.nf'
include { wes_alignment_filtration   }   from "../modules/wes_alignment_filtration.nf"
include { wes_post_alignment         }   from "../modules/wes_post_alignment.nf"
include { wes_annot_index            }   from "../modules/wes_annot_index.nf"

workflow normal_prep_wes {
    take: 
    normal                // Normal sample BAM
    ref_genome            // Reference genome fasta
    ref_files             // Reference genome related files
    bed_files             // Targets BED related files
    targets_bed           // BED file that contains interesting regions coordinates
    ref_dict              // Cleaned reference genome dict file
    cleaned_ref_genome    // Cleaned reference genome fasta

    main: 
    if (params.seq == 'p') {
        wes_trimming(normal)
        alignment(wes_trimming.out, ref_files, ref_genome)
        wes_bam_prep(alignment.out)
        picard_rrg(wes_bam_prep.out)
        wes_alignment_filtration(picard_rrg.out, bed_files, targets_bed, ref_dict)
        wes_post_alignment(wes_alignment_filtration.out, ref_files, cleaned_ref_genome, bed_files, targets_bed, ref_dict)
        wes_annot_index(wes_post_alignment.out.mdup_bam_bai)
    } else if (params.seq == 's') {
        wes_trimming_s(normal)
        alignment_s(wes_trimming_s.out.trimmed_fastq_s, ref_files, ref_genome)
        wes_bam_prep(alignment_s.out.aligned_sam_s)
        picard_rrg(wes_bam_prep.out)
        wes_alignment_filtration(picard_rrg.out, bed_files, targets_bed, ref_dict)
        wes_post_alignment(wes_alignment_filtration.out.filtrated_bam_bai, ref_files, cleaned_ref_genome, bed_files, targets_bed, ref_dict)
        wes_annot_index(wes_post_alignment.out.mdup_bam_bai)
    }
    emit: 
    normal_prep_wes = wes_annot_index.out.mdup_bam_bai
}