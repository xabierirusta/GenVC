// Workflow for the preparation of wgs samples

include { wes_trimming               }   from '../modules/wes_trimming.nf'
include { wes_trimming_s             }   from '../modules/wes_trimming_s.nf'
include { alignment                  }   from '../modules/alignment.nf'
include { alignment_s                }   from '../modules/alignment_s.nf'
include { wes_bam_prep               }   from '../modules/wes_bam_prep.nf'
include { picard_rrg                 }   from '../modules/picard_rrg.nf'
include { wes_alignment_filtration   }   from "../modules/wes_alignment_filtration.nf"
include { wes_post_alignment         }   from "../modules/wes_post_alignment.nf"
include { wes_annot_index            }   from "../modules/wes_annot_index.nf"

workflow prep_wes{
    take:
    fastq                // Fastq file for normal sample
    ref_genome           // Reference genome fasta
    ref_files            // Reference genome related files
    bed_files            // Target BED file related files
    targets_bed          // BED file with interesting region coordinates
    ref_dict             // Reference genome dict 
    cleaned_ref_genome   // Cleaned reference genome

    main: 
    if (params.seq == 'p') {
        wes_trimming(fastq)
        alignment(wes_trimming.out, ref_files, ref_genome)
        wes_bam_prep(alignment.out)
        picard_rrg(wes_bam_prep.out)
        wes_alignment_filtration(picard_rrg.out, bed_files, targets_bed, ref_dict)
        wes_post_alignment(wes_alignment_filtration.out, ref_files, cleaned_ref_genome, bed_files, targets_bed, ref_dict)
        wes_annot_index(wes_post_alignment.out.mdup_bam_bai)
    } else if (params.seq == 's') {
        wes_trimming_s(fastq)
        alignment_s(wes_trimming_s.out.trimmed_fastq_s, ref_files, ref_genome)
        wes_bam_prep(alignment_s.out.aligned_sam_s)
        picard_rrg(wes_bam_prep.out)
        wes_alignment_filtration(picard_rrg.out, bed_files, targets_bed, ref_dict)
        wes_post_alignment(wes_alignment_filtration.out.filtrated_bam_bai, ref_files, cleaned_ref_genome, bed_files, targets_bed, ref_dict)
        wes_annot_index(wes_post_alignment.out.mdup_bam_bai)
    }
    emit: 
    prep_wes = wes_annot_index.out.mdup_bam_bai
}