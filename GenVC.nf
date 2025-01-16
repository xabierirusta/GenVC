#!/bin/env nextflow

// Import processes
include { wes_trimming               }   from './modules/wes_trimming.nf'
include { wes_trimming_s             }   from './modules/wes_trimming_s.nf'
include { wgs_trimming               }   from './modules/wgs_trimming.nf'
include { wgs_trimming_s             }   from './modules/wgs_trimming_s.nf'
include { alignment                  }   from './modules/alignment.nf'
include { alignment_s                }   from './modules/alignment_s.nf'
include { wgs_bam_prep               }   from './modules/wgs_bam_prep.nf'
include { wes_bam_prep               }   from './modules/wes_bam_prep.nf'
include { picard_rrg                 }   from './modules/picard_rrg.nf'
include { wes_alignment_filtration   }   from "./modules/wes_alignment_filtration.nf"
include { wes_post_alignment         }   from "./modules/wes_post_alignment.nf"
include { wgs_post_alignment         }   from "./modules/wgs_post_alignment.nf"
include { wgs_annot_index            }   from "./modules/wgs_annot_index.nf"
include { wes_annot_index            }   from "./modules/wes_annot_index.nf"
include { vc_haplotypecaller         }   from "./modules/vc_haplotypecaller.nf"
include { filter_haplotypecaller     }   from "./modules/filter_haplotypecaller.nf"
include { vc_mutect2                 }   from "./modules/vc_mutect2.nf"
include { filter_mutect2             }   from "./modules/filter_mutect2.nf"
include { vc_delly                   }   from "./modules/vc_delly.nf"
include { filter_delly               }   from "./modules/filter_delly.nf"
include { vc_cnvkit_wes              }   from "./modules/vc_cnvkit_wes.nf"
include { vc_cnvkit_wes_noref        }   from "./modules/vc_cnvkit_wes_noref.nf"
include { vc_cnvkit_wgs              }   from "./modules/vc_cnvkit_wgs.nf"
include { visualize_cnvkit           }   from "./modules/visualize_cnvkit.nf"
include { heatmap_cnvkit             }   from "./modules/heatmap_cnvkit.nf"
include { variant_annotation         }   from "./modules/variant_annotation.nf"
include { final_index                }   from "./modules/final_index.nf"

// Import subworkflows 
include { normal_prep_wes            }   from "./subworkflows/normal_prep_wes.nf"
include { prep_wes                   }   from "./subworkflows/prep_wes.nf"
include { vc_germline_wes_noref      }   from "./subworkflows/vc_germline_wes_noref.nf"
include { vc_germline_wes            }   from "./subworkflows/vc_germline_wes.nf"
include { vc_somatic_wes             }   from "./subworkflows/vc_somatic_wes.nf"
include { normal_prep_wgs            }   from "./subworkflows/normal_prep_wgs.nf"
include { prep_wgs                   }   from "./subworkflows/prep_wgs.nf"
include { vc_germline_wgs            }   from "./subworkflows/vc_germline_wgs.nf"
include { vc_somatic_wgs             }   from "./subworkflows/vc_somatic_wgs.nf"

// Print information about the workflow of the pipeline
def intro(){	
log.info """\
=============================================================================================================================================== 
GenVC - Germline and somatic variant calling pipeline
=============================================================================================================================================== 

Author: Xabier Irusta
Version: 1.0
Last Update: 16/01/2025

Usage: nextflow run GenVC (--input_type) (--input) (--seq) (--v) (--vep_cache) (--outDir)
Optional arguments:
    --help          Shows a help message and exits the pipeline.
    --ref           Directory where reference files are stored, default "./ref/".
    --vcf           Directory where vcf files are stored, default "./vcf/".
    --tools         Directory where files by tools are stored, default "./tools/".
    --t             Number of threads used in the pipeline.
    --normal        Path to the fastq file that is going to be used as normal control (Compulsory when using --v somatic or --input_type WGS).
    --bed           Path to the directory where the BED files are located (For WES, target Bed file HAS to be named targets.bed).

=============================================================================================================================================== 
Parameters Used in the workflow
=============================================================================================================================================== 

input                   = ${params.input}
ref                     = ${params.ref}  
bed                     = ${params.bed}  
tools                   = ${params.tools}  
seq                     = ${params.seq} 
t                       = ${params.t}   
input_type              = ${params.input_type}            
normal                  = ${params.normal} 
v                       = ${params.v} 
vcf                     = ${params.vcf}      
vep_cache               = ${params.vep_cache}   
outDir                  = ${params.outDir}

=============================================================================================================================================== 
"""
}

// Help function
def help_message() {
	log.info"""
	================================================================================================================================= 
	HELP OPTIONS
	=================================================================================================================================
	Pipeline usage: nextflow run GenVC.nf (--input) (--input_type) (--vep_cache) (--outDir)

	Arguments:
	    --input                 Path to the directory where the input FASTQ files are located.
	    --input_type            Type of samples (WGS, WES), default: WGS.
	    --v                     Type of variants to detect (germline, somatic), default: germline.
                                    - Usage of --v somatic means that --normal has to be providaded.
	    --vep_cache             Path to the directory where VEP cache files are located (Previous installation needed).
	    --outDir                Output directory where all the files will be stored.
	       
	Optional arguments:
	    --h                     Shows a help message and exits the pipeline.
	    --t                     Number of threads used in the pipeline.
	    --normal                Path to the location directory of fastq file to be used as normal control.
		                        -Not compulsory for germline WES but REQUIRED for WGS and somatic WES.
	    --seq                   Type of sequencing of the input FASTQ files (p for paired-ended, s for single-ended).
	    --ref                   Path to the directory where the reference files are located.
	    --bed                   Path to the directory where the BED files are located (For WES, target Bed file HAS to be named targets.bed).
	    --tools                 Path to the directory where files needed for Trimmomatic and refFlat.txt are located.
	    --vcf                   Path to the directory where the reference vcf files are located.
	""".stripIndent
}

workflow {
    if(params.help && !params.vep_cache && !params.input) {
        // Invoke help function and exit
        help_message()
        exit 1
    } else if(!params.help && params.vep_cache && params.input) {
        // Channel Generation
        // Channels for each reference file
        Channel
            .fromPath("${params.ref}/*.{fai,amb,ann,bwt,gz,pac,sa}")  // Capture all reference files
            .toList()                                                 // Collect all matching files into a list
            .set { ref_files }                                        
        Channel
            .fromPath("${params.ref}/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
            .first()                                                  // Convert to value channel
            .set { ref_genome }
        Channel
            .fromPath("${params.ref}/cleaned_reference.fa")
            .first() 
            .set { cleaned_ref_genome }
        Channel
            .fromPath("${params.ref}/cleaned_reference.dict")
            .first() 
            .set { ref_dict }
        Channel
            .fromPath("${params.vcf}/00-All_modified.vcf.gz")
            .first()
            .set { dbsnp_ann }
        Channel
            .fromPath("${params.vcf}/af-only-gnomad.hg38.vcf.gz")
            .first()
            .set { gnomad_ann }
        Channel
            .fromPath("${params.vcf}/1000g_pon.hg38.vcf.gz")
            .first()
            .set { pon_hg38 }
        Channel
            .fromPath("${params.vcf}/*.vcf.gz.tbi") 
            .toList() 
            .set { vcf_files } 
        Channel
            .fromPath("${params.bed}/*.{bed.gz,bed.gz.tbi}") 
            .toList() 
            .set { bed_files } 
        Channel
            .fromPath("${params.bed}/targets.bed") 
            .first() 
            .set { targets_bed }  
        Channel
            .fromPath("${params.bed}/access-10kb.hg38.bed")
            .first() 
            .set { access_bed }
        // Channel for input fastq files
        if (params.seq == "p") {
            // Match *_1.fastq.gz and *_2.fastq.gz files for paired-end data
            fastq_files = Channel.fromFilePairs("${params.input}/*_{1,2}.fastq.gz", flat: false, checkIfExists: true)
        } else {
            // For single-end data, match *.fastq.gz files
            fastq_files = Channel.fromPath("${params.input}/*.fastq.gz", checkIfExists: true)
                .map { file -> [file] } // Wrap single file in a list
        }
        // Channel for normal samples
        if (params.normal) {
            if (params.seq == "p") {
                normal = Channel.fromFilePairs("${params.normal}/*_{1,2}.fastq.gz", flat: false, checkIfExists: true)
            } else {
                normal = Channel.fromPath("${params.normal}/*.fastq.gz", checkIfExists: true)
                .map { file -> [file] } // Wrap single file in a list
            }
        }
        cache_dir = Channel.fromPath("${params.vep_cache}").first()
        Channel
            .fromPath("${params.tools}/refFlat.txt")
            .first()
            .set { refFlat }

        Channel
            .fromPath("${params.tools}/contigs.list")
            .first()
            .set { contigs }
				
        // INITIALIZE THE PROCESSES
        intro() // Print the introduction message
        if (params.input_type == 'WES') {
            fastq_files.view { "Input files: ${it}" }
            // WES workflow with normal control
            if (params.normal) {
                normal.view { "Input normal: ${it}" }
                normal_prep_wes(normal, ref_genome, ref_files, bed_files, targets_bed, ref_dict, cleaned_ref_genome)
                prep_wes(fastq_files, ref_genome, ref_files, bed_files, targets_bed, ref_dict, cleaned_ref_genome)
                if (params.v == 'germline') {
                    vc_germline_wes(prep_wes.out.prep_wes, ref_files, ref_dict, cleaned_ref_genome, dbsnp_ann, refFlat, bed_files, targets_bed, access_bed, cache_dir, normal_prep_wes.out.normal_prep_wes, vcf_files)
                } else if (params.v == 'somatic') {
                    vc_somatic_wes(prep_wes.out.prep_wes, ref_files, ref_dict, gnomad_ann, pon_hg38, cleaned_ref_genome, dbsnp_ann, refFlat, bed_files, targets_bed, access_bed, cache_dir, normal_prep_wes.out.normal_prep_wes, vcf_files, contigs)       
                }
            // WES worflow but without reference
            } else {
                prep_wes(fastq_files, ref_genome, ref_files, bed_files, targets_bed, ref_dict, cleaned_ref_genome)
                if (params.v == 'germline') {
                    vc_germline_wes_noref(prep_wes.out.prep_wes, ref_files, ref_dict, cleaned_ref_genome, dbsnp_ann, refFlat, bed_files, targets_bed, access_bed, cache_dir, vcf_files)
                }
            }
        } else if (params.input_type == 'WGS') {
            fastq_files.view { "Input files: ${it}" }
            // WGS workflow with normal control
            if (params.normal) {
                normal.view { "Input normal: ${it}" }
                normal_prep_wgs(normal, ref_genome, ref_files, ref_dict)
                prep_wgs(fastq_files, ref_genome, ref_files, ref_dict)
                if (params.v == "germline"){
                    vc_germline_wgs(prep_wgs.out.prep_wgs, normal_prep_wgs.out.normal_prep_wgs, ref_files, dbsnp_ann, bed_files, refFlat, cleaned_ref_genome, cache_dir, targets_bed, ref_dict, vcf_files)
                } else if (params.v == 'somatic') {
                    vc_somatic_wgs(prep_wgs.out.prep_wgs, normal_prep_wgs.out.normal_prep_wgs, ref_files, dbsnp_ann, gnomad_ann, pon_hg38, ref_dict, refFlat, bed_files, targets_bed, cleaned_ref_genome, cache_dir, vcf_files, contigs)
                }
            }   
        } 
    } else {
        // Invoke help function and exit
        help_message
        exit 1
    }
}

workflow.onComplete {

summary = """
=============================================================================================================================================== 
Workflow summary
=============================================================================================================================================== 

Duration		: ${workflow.duration}
Success			: ${workflow.success}
workDir			: ${workflow.workDir}
Exit status	    : ${workflow.exitStatus}
outDir			: ${params.outDir}

=============================================================================================================================================== 

"""
println summary
}
