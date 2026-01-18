version 1.0

workflow wf_bordetella_typing {
  input {
    File assembly_file
    String samplename
  }

  call typing {
    input:
        assembly_file = assembly_file,
        samplename    = samplename
  }

  output {      
    File alleles_ids_files = typing.alleles_ids_files
    File alleles_ids_sum   = typing.alleles_ids_sum

    # Macrolide_resistance
    String rRNA_23S_mr          = typing.rRNA_23S
    String fhaB_mr              = typing.fhaB
    String fhaB_2400_5550_mr    = typing.fhaB_2400_5550
    String fhaBx_1_3193_mr      = typing.fhaBx_1_3193
    String fhaBy_3190_7183_mr   = typing.fhaBy_3190_7183
    String fhaBz_7180_10773_mr  = typing.fhaBz_7180_10773
    # Phase
    String bipA_phase              = typing.bipA
    String bvgS_phase              = typing.bvgS
    String bvgA_phase              = typing.bvgA
    # Locus
    String cyaA_toxins               = typing.cyaA
    String dnt_toxins               = typing.dnt
    # T3SS
    String T3SStype          = typing.T3SStype
    String bopB              = typing.bopB
    String bopN              = typing.bopN
    String bteA              = typing.bteA
    String bsp22             = typing.bsp22
    String bopD              = typing.bopD
    # Autotransporters:
    String brkA_autotransporters              = typing.brkA
    String prn_autotransporters              = typing.prn
    String bapC_autotransporters              = typing.bapC
    String vag8_autotransporters              = typing.vag8
    String tcfA_autotransporters              = typing.tcfA
    # Bp_vaccine_antigens
    String ptxA              = typing.ptxA
    String ptxB              = typing.ptxB
    String ptxC              = typing.ptxC
    String ptxD              = typing.ptxD
    String ptxE              = typing.ptxE
    String ptxP              = typing.ptxP
    String fim3              = typing.fim3
    String fim2              = typing.fim2
    String BPagST            = typing.BPagST
    # PRN-test-Bp
    # String PRN_Bp            = typing.PRN_Bp

  }


}

task typing {
    input {
        File assembly_file
        String samplename
    }

    command <<< 
        
        bordetella_typing \
            --assembly ~{assembly_file} \
            --output_dir ~{samplename} \
            --output_ids_file 'alleles_ids.txt'

        tar -czf ~{samplename}.tar.gz ~{samplename}
    
        # Alleles processed:
        # Macrolide_resistance:
        #   23S_rRNA, fhaBz-7180_10773,fhaBy-3190_7183, fhaBx-1_3193, fhaB-2400_5550, fhaB
        # Phase:
        #   bipA, bvgS, bvgA
        # Locus:
        #   cyaA, dnt
        # T3SS:
        #   T3SStype, bopN, bopB, bteA, bsp22, bopD
        # Autotransporters:
        #   brkA, prn, vag8, bapC, tcfA
        # Bp_vaccine_antigens:
        #   ptxC, ptxA, ptxP, fim3, ptxE ptxB, ptxD, BPagST
        #   (fhaB-2400_5550 - already in Macrolide_resistance)
        # PRN-test-Bp:
        #    PRN-Bp


        # All expected alleles
        alleles=(
            23S_rRNA fhaBz-7180_10773 fhaBy-3190_7183 fhaBx-1_3193 fhaB-2400_5550 fhaB
            bipA bvgS bvgA
            cyaA dnt
            T3SStype bopN bopB bteA bsp22 bopD
            brkA prn vag8 bapC tcfA
            ptxC ptxA ptxP fim3 ptxE ptxB ptxD BPagST
        )

        # Track alleles we processed
        processed=()

        while IFS=':' read -r allele_name allele_id; do
            allele_name=$(echo "$allele_name" | xargs)
            allele_id=$(echo "$allele_id" | xargs)

            if [[ -z "$allele_name" ]]; then
                continue
            fi

            safe_name=$(echo "$allele_name" | sed 's#[^a-zA-Z0-9._-]#_#g')
            echo "$allele_id" > "${safe_name}.txt"

            processed+=("$allele_name")
        done < "alleles_ids.txt"

        # Create empty files for missing alleles
        for allele in "${alleles[@]}"; do
            if [[ ! " ${processed[@]} " =~ " ${allele} " ]]; then
                safe_name=$(echo "$allele" | sed 's#[^a-zA-Z0-9._-]#_#g')
                echo "" > "${safe_name}.txt"
            fi
        done
    >>>

    output {
        File alleles_ids_files = "~{samplename}.tar.gz"
        File alleles_ids_sum   = "alleles_ids.txt"
        
        # Macrolide_resistance
        String rRNA_23S          = read_string("23S_rRNA.txt")
        String fhaB              = read_string("fhaB.txt")
        String fhaB_2400_5550    = read_string("fhaB-2400_5550.txt")
        String fhaBx_1_3193      = read_string("fhaBx-1_3193.txt")
        String fhaBy_3190_7183   = read_string("fhaBy-3190_7183.txt")
        String fhaBz_7180_10773  = read_string("fhaBz-7180_10773.txt")

        # Phase
        String bipA           = read_string("bipA.txt")
        String bvgS           = read_string("bvgS.txt")
        String bvgA           = read_string("bvgA.txt")

        # Locus
        String cyaA         = read_string("cyaA.txt")
        String dnt          = read_string("dnt.txt")

        # T3SS
        String T3SStype             = read_string("T3SStype.txt")
        String bopB                 = read_string("bopB.txt")
        String bopN                 = read_string("bopN.txt")
        String bteA                 = read_string("bteA.txt")
        String bsp22                = read_string("bsp22.txt")
        String bopD                 = read_string("bopD.txt")

        # Autotransporters
        String brkA    = read_string("brkA.txt")
        String prn     = read_string("prn.txt")
        String bapC    = read_string("bapC.txt")
        String vag8    = read_string("vag8.txt")
        String tcfA    = read_string("tcfA.txt")
        
        # Bp_vaccine_antigens
        String ptxA              = read_string("ptxA.txt")
        String ptxB              = read_string("ptxB.txt")
        String ptxC              = read_string("ptxC.txt")
        String ptxD              = read_string("ptxD.txt")
        String ptxE              = read_string("ptxE.txt")
        String ptxP              = read_string("ptxP.txt")
        String fim3              = read_string("fim3.txt")
        String fim2              = read_string("fim2.txt")
        String BPagST            = read_string("BPagST.txt")
        
        # PRN-test-Bp
        # String PRN_Bp            = read_string("PRN-Bp.txt")

    }

    runtime {
        docker: "bioinfomoh/bordetella_typing:1"
        maxRetries: 2

    }

}


