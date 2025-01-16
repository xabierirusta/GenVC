// Workflow for germline variant calling with a reference sample
// Import processes
include { vc_haplotypecaller         }   from "../modules/vc_haplotypecaller.nf"
include { filter_haplotypecaller     }   from "../modules/filter_haplotypecaller.nf"
include { vc_delly                   }   from "../modules/vc_delly.nf"
include { filter_delly               }   from "../modules/filter_delly.nf"
include { vc_cnvkit_wgs              }   from "../modules/vc_cnvkit_wgs.nf"
include { visualize_cnvkit           }   from "../modules/visualize_cnvkit.nf"
include { heatmap_cnvkit             }   from "../modules/heatmap_cnvkit.nf"
include { variant_annotation         }   from "../modules/variant_annotation.nf"
include { final_index                }   from "../modules/final_index.nf"

workflow vc_germline_wgs {
    take:
    mdup                // Markduplicates BAM and BAI files for the samples
    norm_mdup           // Markduplicates BAM and BAI files for the normal control
    ref_files           // Reference genome related files
    dbsnp_ann           // dbSNP database VCF file
    bed_files           // Targets BED related files
    refFlat             // refFlat.txt file
    cleaned_ref_genome  // Cleaned reference genome fasta
    cache_dir           // Anntation cache directory for Ensembl VEP
    targets_bed         // BED file with regions of interest
    ref_dict            // Reference genome dict
    vcf_files           // VCF files needed in the workflow

    main:
    normal_channel = norm_mdup.first()
    vc_haplotypecaller(mdup, ref_files, dbsnp_ann, cleaned_ref_genome, ref_dict, vcf_files)
    filter_haplotypecaller(vc_haplotypecaller.out.htcaller_vcf, ref_files, cleaned_ref_genome, ref_dict, vcf_files)
    vc_delly(mdup, ref_files, cleaned_ref_genome)
    filter_delly(vc_delly.out.delly_vcf, ref_files, cleaned_ref_genome, ref_dict)
    vc_cnvkit_wgs(mdup, ref_files, refFlat, bed_files, cleaned_ref_genome, targets_bed, normal_channel)
    visualize_cnvkit(vc_cnvkit_wgs.out.cnvkit_output)
    // Intermediate Channels needed in the heatmap generation
    wgs_cns_files = vc_cnvkit_wgs.out.cnvkit_output
        .map { it[1] }  // it[1] extracts the file part of the tuple
        .collect()      // Collect all the files into a list
        .view { "CNVKit Heatmap files: ${it}" }          // View the collected files
    heatmap_cnvkit(wgs_cns_files)
    vcf_mix = filter_haplotypecaller.out.filt_htc_vcf
        | mix(vc_cnvkit_wgs.out.cnvkit_vcf, filter_delly.out.filt_delly_vcf)
        | groupTuple(size: 3)
        | map { tuple ->
            // Access the ID (first element of the tuple) and sort the files (second element)
            def id = tuple[0]
            def sorted_files = tuple[1].sort { a, b -> 
                a.getFileName().toString() <=> b.getFileName().toString()  // Sort by file name
            }
            return [id, sorted_files]  // Return the ID and the sorted list of files   
        }
        | view { "Filtered VCF files: ${it}" }  // For debugging, view the sorted result

    variant_annotation(vcf_mix, cache_dir)
    ann_mix = variant_annotation.out.gatk_ann
        | mix(variant_annotation.out.delly_ann, variant_annotation.out.cnvkit_ann)
        | groupTuple(size: 3)
        | map { tuple ->
            // Access the ID (first element of the tuple) and sort the files (second element)
            def id = tuple[0]
            def sorted_files = tuple[1].sort { a, b -> 
                a.getFileName().toString() <=> b.getFileName().toString()  // Sort by file name
            }
            return [id, sorted_files]  // Return the ID and the sorted list of files   
        }
        | view { "Annotated VCF files: ${it}" }   // For debugging, view the sorted result
    final_index(ann_mix) 

    emit:
    gatk_vcf = final_index.out.gatk_ann
    delly_vcf = final_index.out.delly_ann
    cnvkit_vcf = final_index.out.cnvkit_ann

}