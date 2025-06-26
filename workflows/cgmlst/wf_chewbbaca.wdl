version 1.0

import "../../tasks/cgmlst/task_chewbbaca.wdl" as task_chewbbaca

workflow wf_chewbbaca {
    input {
        File schema_zip
        File prodigal_zip
        Array[File] input_assemblies
    }

    call task_chewbbaca.allele_calling {
        input:
            schema_zip = schema_zip,
            prodigal_zip = prodigal_zip,
            assemblies = input_assemblies
    }

    call task_chewbbaca.extract_cgmlst {
        input:
            cleaned_results = allele_calling.alleles_cleaned,
    }

    output {
        File chew_alleles = allele_calling.alleles_cleaned
        File chew_visualization = extract_cgmlst.visualization
    }
}