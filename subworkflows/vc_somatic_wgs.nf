// Workflow for somatic variant calling with a reference sample
// Import processes
include { vc_mutect2                 }   from "../modules/vc_mutect2.nf"
include { filter_mutect2             }   from "../modules/filter_mutect2.nf"
include { vc_delly                   }   from "../modules/vc_delly.nf"
include { filter_delly               }   from "../modules/filter_delly.nf"
include { vc_cnvkit_wgs              }   from "../modules/vc_cnvkit_wgs.nf"
include { visualize_cnvkit           }   from "../modules/visualize_cnvkit.nf"
include { heatmap_cnvkit             }   from "../modules/heatmap_cnvkit.nf"
include { variant_annotation         }   from "../modules/variant_annotation.nf"
include { final_index                }   from "../modules/final_index.nf"

workflow vc_somatic_wgs {
    take:
    mdup                // Markduplicates BAM and BAI files for the samples
    normal              // Markduplicates BAM and BAI files for the normal control
    ref_files           // Reference genome related files
    dbsnp_ann           // dbSNP database annotation file
    gnomad_ann          // genomAD database annotation file
    pon_hg38            // Panel of Normals for hg38 for mutect2
    ref_dict            // Reference genome dict
    refFlat             // refFlat.txt file
    bed_files           // Targets BED related files
    targets_bed         // BED file containing regions of interest
    cleaned_ref_genome  // Cleaned reference genome
    cache_dir           // Annotations cache directory for Ensembl VEP
    vcf_files           // VCF files needed in the workflow
    contigs             // Contigs list contianing chromosomes of interest

    main:
    normal_channel = normal.first()
    vc_mutect2(mdup, normal_channel, ref_files, cleaned_ref_genome, dbsnp_ann, gnomad_ann, pon_hg38, ref_dict, vcf_files) 
    filtering_input = mdup
        .mix(vc_mutect2.out.mutect2_vcf)  // Mix the channels first
        .groupTuple()  // Group the items by the ID
        .map { tuple -> 
            // Extract the ID and the file lists from the tuple
            def id = tuple[0]  // The ID (e.g., SRR7890851)
            def bam_vcf = tuple[1]  // The list with BAM and VCF files
            def bai = tuple[2]  // The list with the BAI file

            // Extract BAM and VCF files
            def bam = bam_vcf[0]  // The BAM file
            def vcf = bam_vcf[1]  // The VCF file

            // Flatten the BAM and BAI files into one list
            def sorted_files = [bam, bai[0]].sort { a, b -> 
                a.getFileName().toString() <=> b.getFileName().toString()  // Sort by file name
            }

            // Return the final structure: [id, [bam, bai, vcf]]
            return [id, sorted_files + [vcf]]  // Concatenate sorted BAM/BAI with VCF
        }
        .view { "Filtered mutect2: ${it}" }  // For debugging, view the sorted result
    filter_mutect2(filtering_input, gnomad_ann, ref_files, cleaned_ref_genome, ref_dict, vcf_files, contigs, vc_mutect2.out.mutect2_stats)
    vc_delly(mdup, ref_files, cleaned_ref_genome)
    filter_delly(vc_delly.out.delly_vcf, ref_files, cleaned_ref_genome, ref_dict)
    vc_cnvkit_wgs(mdup, ref_files, refFlat, bed_files, cleaned_ref_genome, targets_bed, normal_channel)
    visualize_cnvkit(vc_cnvkit_wgs.out.cnvkit_output)
    // Intermediate Channels needed in the heatmap generation
    cns_files = vc_cnvkit_wgs.out.cnvkit_output
        .map { it[1] }  // it[1] extracts the file part of the tuple
        .collect()      // Collect all the files into a list
        .view { "CNVKit Heatmap files: ${it}" }        // View the collected files
    heatmap_cnvkit(cns_files)
    vcf_mix = filter_mutect2.out.filt_mutect2_vcf
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