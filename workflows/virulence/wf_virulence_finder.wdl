version 1.0

workflow wf_virulence_finder {
    input {
        File fasta_file
        String species        # ecoli  or salmonella 
        String samplename     
        Boolean plasmid_check = false
    }

    meta {
		description: "The code analyzes virulence gene presence from genomic data using the FDA database (https://virulence.fda.gov/) and generates an HTML report with styled tables."
		author: "David Maimoun"
        email: "david.maimoun@moh.gov.il"
        version: "1.0"
	}

    parameter_meta {
        fasta_file: "Fasta assembly file"
        species: "Two options: ecoli or salmonella"
        plasmid_check: "Check plasmid virulence factors. Boolean by default 'false'"
    }

    call run_virulence_finder {
        input:
        fasta_file = fasta_file,
        species = species,
        samplename = samplename,
        plasmid_check = plasmid_check
    }

    output {
            File virfind_hits_csv     = run_virulence_finder.hits_csv
            File virfind_summary_csv  = run_virulence_finder.summary_csv
            File virfind_summary_html = run_virulence_finder.summary_html

            File? virfind_plasmid_hits_csv     = run_virulence_finder.plasmid_hits_csv
            File? virfind_plasmid_summary_csv  = run_virulence_finder.plasmid_summary_csv
            File? virfind_plasmid_summary_html = run_virulence_finder.plasmid_summary_html
    }
}

task run_virulence_finder {
    input {
        File fasta_file
        String species
        String samplename
        Boolean plasmid_check
    }

    command <<<
        virulence_finder \
            --species ~{species} \
            --samplename ~{samplename} \
            --fasta ~{fasta_file} \
            ~{if plasmid_check then "--plasmid_check" else ""}
    >>>

    output {
        File hits_csv     = "~{samplename}_virulence_hits.csv"
        File summary_csv  = "~{samplename}_virulence_summary.csv"
        File summary_html = "~{samplename}_virulence_summary.html"

        File? plasmid_hits_csv     = "~{samplename}_plasm_virulence_hits.csv"
        File? plasmid_summary_csv  = "~{samplename}_plasm_virulence_summary.csv"
        File? plasmid_summary_html = "~{samplename}_plasm_virulence_summary.html"
    }

  runtime {
    docker: "bioinfomoh/virulence_finder:1"
    memory: "4G"
    cpu: "12"
  }
}
