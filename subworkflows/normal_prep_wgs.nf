// Workflow for the preparation of normal wgs samples
// Import processes
include { wgs_trimming               }   from '../modules/wgs_trimming.nf'
include { wgs_trimming_s             }   from '../modules/wgs_trimming_s.nf'
include { alignment                  }   from '../modules/alignment.nf'
include { alignment_s                }   from '../modules/alignment_s.nf'
include { wgs_bam_prep               }   from '../modules/wgs_bam_prep.nf'
include { picard_rrg                 }   from '../modules/picard_rrg.nf'
include { wgs_post_alignment         }   from "../modules/wgs_post_alignment.nf"
include { wgs_annot_index            }   from "../modules/wgs_annot_index.nf"
workflow normal_prep_wgs {
    take:
    normal                // Fastq file for normal sample
    ref_genome            // Reference genome fasta
    ref_files             // Reference genome related files
    ref_dict              // Reference genome dict file

    main:
    if (params.seq == 'p') {
        wgs_trimming(normal)
        alignment(wgs_trimming.out, ref_files, ref_genome)
        wgs_bam_prep(alignment.out, ref_dict)
        picard_rrg(wgs_bam_prep.out)
        wgs_post_alignment(picard_rrg.out, ref_files, ref_genome)
        wgs_annot_index(wgs_post_alignment.out.mdup_bam_bai)
    } else if (params.seq == 's') {
        wgs_trimming_s(normal)
        alignment_s(wgs_trimming_s.out.trimmed_fastq_s, ref_files, ref_genome)
        wgs_bam_prep(alignment_s.out, ref_dict)
        picard_rrg(wgs_bam_prep.out)
        wgs_post_alignment(picard_rrg.out, ref_files, ref_genome)
        wgs_annot_index(wgs_post_alignment.out.mdup_bam_bai)
    }

    emit:
    normal_prep_wgs = wgs_annot_index.out.mdup_bam_bai

}