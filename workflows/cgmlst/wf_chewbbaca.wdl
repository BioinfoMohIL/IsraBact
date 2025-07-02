version 1.0

import "../../tasks/cgmlst/task_chewbbaca.wdl" as task_chewbbaca
import "../../tasks/task_fail.wdl" as task_fail

workflow wf_chewbbaca {
    input {
        File schema_zip
        File prodigal_zip
        Array[File]? input_assemblies
        File? assemblies_zipped
    }

    Boolean assemblies_valid = 
        (defined(input_assemblies) && !defined(assemblies_zipped)) ||
        (!defined(input_assemblies) && defined(assemblies_zipped))

    if (!assemblies_valid) {
        String error = "You must provide EITHER 'input_assemblies' OR 'assemblies_zipped', but not both."
        call task_fail.fail {
            message = error
        }
    } 
    else {
        call task_chewbbaca.allele_calling {
            input:
                schema_zip = schema_zip,
                prodigal_zip = prodigal_zip,
                assemblies = input_assemblies,
                assemblies_zipped = assemblies_zipped
        }

        call task_chewbbaca.extract_cgmlst {
            input:
                cleaned_results = allele_calling.alleles_cleaned
        }

        output {
            File chew_alleles = allele_calling.alleles_cleaned
            File chew_visualization = extract_cgmlst.visualization
        }
    }
}
