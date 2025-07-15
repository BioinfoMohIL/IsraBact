version 1.0


workflow listeria_pred_reads {

    meta {
        author: "David Maimoun"
        organization : "MOH"
        description: "Runs ListPred (https://github.com/genomicepidemiology/ListPred.git)on raw data to predict Listeria serotype."
        version: "1.0"
    }
    
    input {
        File read1
        File read2          
        Int cpu = 12          
    }

    call listeria_pred {
        input:
            read1 = read1,
            read2 = read2,
            cpu = cpu,
    }

    output {
        File test = listeria_pred.test
        # String  listpred_virulence              = listeria_pred.virulence_class
        # String  listpred_disinftolerance        = listeria_pred.disinfectant_phenotype
        # File    listpred_virulence_pred         = listeria_pred.virulence_pred
        # File    listpred_combined_pred_out_cat  = listeria_pred.combined_pred_out_cat
        # File    listpred_combined_pred_out_num  = listeria_pred.combined_pred_out_num
        # File    listpred_disinftolerance_pred   = listeria_pred.disinftolerance_pred
        # File?   listpred_disinf_align           = listeria_pred.disinf_align
        # File?   listpred_vir_align              = listeria_pred.vir_align
    }
}

task listeria_pred {
    input {
        File read1
        File read2
        Int cpu

    }

    command <<<
        
 
        cd /ListPred
        cat ~{read1} > read1.fastq.gz
        cat ~{read2} > read2.fastq.gz
        
        snakemake -s /ListPred/workflow/Snakefile --cores ~{cpu} --use-conda --config ipe="read1.fastq.gz read2.fastq.gz" outd="pred_results"
        echo 'her' >> t.txt

    >>>


    output {
        File test = 't.txt'
        # String virulence_class        = read_string('virulence_class.txt')
        # String disinfectant_phenotype = read_string('disinfectant_phenotype.txt')
        # File virulence_pred           = "virulence_prediction_out.csv"
        # File combined_pred_out_cat    = "combined_predictions_out_categorical.csv"  
        # File combined_pred_out_num    = "combined_predictions_out_numerical.csv"  
        # File disinftolerance_pred     = "disinftolerance_prediction_out.csv"
        # File? disinf_align            = "disinf_align_out.tar.gz"
        # File? vir_align               = "vir_align_out.tar.gz"


    }

    runtime {
        docker: "genomicepidemiology/listpred:0.2.0"
        cpu: cpu 
        memory: "8 GB" 
    }

}
