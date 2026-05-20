version 1.0

import "../../tasks/wgmlst/task_chewbbaca.wdl" as task_chewbbaca
import "../../tasks/task_fail.wdl" as task_fail

workflow wf_chewbbaca {
    input {
        Array[File] assemblies
        File   schema_adapted
        File?  schema_adapted_zip
        File? assemblies_zipped


        ## Optionals ##
        
        ## ── PrepSchema ──────────────────────────────────────
        # Boolean run_prep_schema   = false
        # File?   external_schema_zip
        # File?   prodigal_training_file   

        ## ── AlleleCall ──────────────────────────────────────
        Int  allele_call_mode                   = 4
        Float   allele_call_bsr                 = 0.6
        Int     allele_call_minimum_length      = 0
        Int     allele_call_translation_table   = 11
        Float   allele_call_size_threshold      = 0.2
        Boolean allele_call_cds_input           = false
        Boolean allele_call_no_inferred         = false

        ## ── ExtractCgMLST ───────────────────────────────────
        Array[Float] cgmlst_thresholds = [0.95, 0.99, 1]

        # ## ── Modules optionnels ──────────────────────────────
        # Boolean run_schema_evaluator       = false
        # Boolean run_allele_call_evaluator  = false

        # Boolean run_remove_genes = false
        # File?   genes_to_remove          # liste de loci à exclure
        # Boolean remove_genes_inverse = false  # si true: garder les loci listés

        # Boolean run_get_alleles    = false
        # Boolean get_alleles_distinct   = false
        # Boolean get_alleles_translate  = false

        # Boolean run_join_profiles  = false
        # Array[File] extra_profiles = []  # TSV de runs précédents à joindre

        # Boolean run_uniprot_finder = false
        # File?   protein_table               # cds_coordinates.tsv
        # String? uniprot_taxa                # ex: "Escherichia coli"

        # Boolean run_compute_msa    = false
    }

    Boolean assemblies_valid = 
        (defined(assemblies) && !defined(assemblies_zipped)) ||
        (!defined(assemblies) && defined(assemblies_zipped))

    if (!assemblies_valid) {
        call task_fail.fail {
            input:
                message = "You must provide EITHER 'assemblies' OR 'assemblies_zipped', but not both."
        }
    } 

    # TODO : check if needed, for now, dont use it
    # if (run_prep_schema) {
    #     call chewbacca.prep_external_schema {
    #         input:
    #             external_schema_zip = select_first([external_schema_zip]),
    #             training_file       = prodigal_training_file,
    #             bsr                 = bsr,
    #             minimum_length      = minimum_length,
    #             translation_table   = translation_table,
    #             size_threshold      = size_threshold,
    #             cpu                 = cpu,
    #             disk_gb             = 100,
    #             memory_gb           = 16
    #     }
    
    # }


    call task_chewbbaca.allele_calling {
        input:
            schema_adapted = schema_adapted
            schema_adapted_zip = schema_adapted_zip,
            assemblies = assemblies,
            assemblies_zipped = assemblies_zipped,
            allele_call_mode          = allele_call_mode,
            allele_call_cds_input     = allele_call_cds_input,
            allele_call_no_inferred    = allele_call_no_inferred,
            allele_call_bsr               = allele_call_bsr,
            allele_call_minimum_length    = allele_call_minimum_length,
            allele_call_translation_table = allele_call_translation_table,
            allele_call_size_threshold    = allele_call_size_threshold
    }

    call task_chewbbaca.extract_cgmlst {
        input:
            cleaned_results = allele_calling.calling_alleles_cleaned,
            cgmlst_thresholds = cgmlst_thresholds
    }

    output {
        String chew_version = "3.3.10"
        File? fail_logs = fail.fail_logs
        File chew_alleles = allele_calling.calling_alleles_cleaned
        File chew_visualization_zip = extract_cgmlst.visualization_zip
    }
    
}
