version 1.0

import "../../tasks/metagenomics/task_detect_species.wdl" as metagenomics

workflow species_detection {
    meta {
        description: "Detect your species (Kraken2)."
        author: "David Maimoun (The Codon Bleu)"
        email: "thecodonbleu@outlook.com"
    }

    input {
        File read1
        File? read2
        String samplename
    }


    call metagenomics.detect_species  {
        input:
            read1 = read1,
            read2 = read2,
            sample_id = samplename
    }
        

    output {
        String meta_version = detect_species.version
        File meta_kraken_report = detect_species.report
        String meta_species_detected = detect_species.species_detected
    
    }

}
  