version 1.0  

task neisseria_typing {
    input {
        File assembly
    }
    
    command <<<
        mkdir results

        echo "NeisseriaTyping v1.0" > VERSION

        val=("st" "cc" "FetA_VR" "bast_type" "fHbp_peptide" "NHBA_peptide" "NadA_peptide" "PorA_VR1" "PorA_VR2" "abcZ" "adk" "aroE" "fumC" "gdh" "pdhC" "pgm" "bexsero_cross_reactivity" "trumenba_cross_reactivity")

        # Join them into a comma-separated string for argparse
        files_arg=$(IFS=, ; echo "${val[*]}")

        neisseria_typing --input ~{assembly} --output results --output_csv_file typing_report.csv --files "$files_arg" --split
    >>>

    output {
        String version = read_string("VERSION")
        File typing_results = 'results/typing_report.csv'
        String fHbp_peptide = read_string("results/fHbp_peptide.txt")
        String NHBA_peptide = read_string("results/NHBA_peptide.txt") 
        String NadA_peptide = read_string("results/NadA_peptide.txt") 
        String PorA_VR1 = read_string("results/PorA_VR1.txt") 
        String PorA_VR2 = read_string("results/PorA_VR2.txt") 
        String abcZ = read_string("results/abcZ.txt") 
        String adk = read_string("results/adk.txt") 
        String aroE = read_string("results/aroE.txt") 
        String gdh = read_string("results/gdh.txt") 
        String fumC = read_string("results/fumC.txt") 
        String pdhC = read_string("results/pdhC.txt") 
        String pgm = read_string("results/pgm.txt") 
        String st = read_string("results/st.txt") 
        String clonal_complex = read_string("results/clonal_complex.txt") 
        String FetA_VR = read_string("results/FetA_VR.txt") 
        String bast_type = read_string("results/bast_type.txt") 
        String bexsero_cross_reactivity = read_string("results/bexsero_cross_reactivity.txt") 
        String trumenba_cross_reactivity = read_string("results/trumenba_cross_reactivity.txt") 
        
    }

    runtime {
        docker: "bioinfomoh/neisseria_typing:1"
        memory: "8 GB"
        maxRetries: 2

    }
}

task serogrouping {
    input {
        String samplename
        File assembly
    }

    command <<<
        python3 /app/nm_completion.py \
        --input ~{assembly} \
        --serogroup_filename ~{samplename}_serogroup.txt \
        --genogroup_filename ~{samplename}_genogroup.txt \
    >>>

    output {
        String serogroup = read_string("~{samplename}_serogroup.txt")
        String genogroup = read_string("~{samplename}_genogroup.txt")
    }

    runtime {
        docker: "bioinfomoh/nm_completion_analysis:1"
        cpu: 8
        maxRetries: 3

    }
}
