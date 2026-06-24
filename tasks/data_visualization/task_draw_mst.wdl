version 1.0

## ════════════════════════════════════════════════════════════
## task_draw_mst.wdl
## Génère un MST style GrapeTree (PNG) depuis un .nwk
## + métadonnées optionnelles (TSV).
## ════════════════════════════════════════════════════════════

task draw_mst_network {
    input {
        File    nwk_file                            # .nwk produit par GrapeTree

        File?   metadata_tsv                        # TSV optionnel (sample_id, year, ...)
        String  color_col       = "year"            # colonne de couleur
        String  id_col          = "sample_id"       # colonne ID dans le TSV
        String? count_col                           # colonne pour la taille des nœuds

        String  title           = "Minimum Spanning Tree — cgMLST"
        String  output_name     = "mst_network.png"
        Int     dpi             = 200

        Int     cpu             = 2
        Int     memory_gb       = 4
        Int     disk_gb         = 20
    }

    command <<<
        set -euo pipefail

        META_ARG=""
        if [ -s "~{metadata_tsv}" ]; then
            META_ARG="--meta ~{metadata_tsv} --id-col ~{id_col} --color-col ~{color_col}"
        fi

        COUNT_ARG=""
        if [ -n "~{count_col}" ]; then
            COUNT_ARG="--count-col ~{count_col}"
        fi

        echo "[run] draw_mst_network"
        python /app/draw_mst_network.py \
            --input  "~{nwk_file}" \
            --output "~{output_name}" \
            --title  "~{title}" \
            --dpi    ~{dpi} \
            $META_ARG \
            $COUNT_ARG

        echo "[done] $(ls -lh ~{output_name})"
    >>>

    output {
        File mst_png = output_name
    }

    runtime {
        docker:     "bioinfomoh/draw_mst_network:1.0"   
        cpu:        cpu
        memory:     "~{memory_gb} GB"
        disks:      "local-disk ~{disk_gb} HDD"
        maxRetries: 1
    }
}
