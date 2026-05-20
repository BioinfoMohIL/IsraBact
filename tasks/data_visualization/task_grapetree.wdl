version 1.0

task grapetree_mst {
    input {
        File cgmlst_matrix  # TSV file from chewBBACA extract_cgmlst
        String tree_method = "MSTreeV2"  # MSTreeV2 (default), MSTree, or NJ
        String output_prefix = "cgmlst_tree"
        
        # Docker image with GrapeTree installed
        String docker_image = "bioinfomoh/grapetree:1"
        Int cpu = 2
        Int memory_gb = 4
        Int disk_size_gb = 10
    }

    command <<<
        set -euxo pipefail

        # Run GrapeTree to generate MST
        # Input: cgMLST matrix (TSV format)
        # Output: Newick tree file
        grapetree \
            --method ~{tree_method} \
            --profile ~{cgmlst_matrix} > ~{output_prefix}.nwk

        if [[ ! -f "~{output_prefix}.nwk" ]]; then
            echo "Error: Newick file not generated"
            exit 1
        fi

     
        echo "✓ GrapeTree nwk generation complete"
    >>>

    output {
        File tree_newick = "~{output_prefix}.nwk"
    }

    runtime {
        docker: docker_image
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_size_gb} HDD"
    }

    meta {
        author: "David Maimoun"
        description: "Generate Minimum Spanning Tree (MST) visualization from cgMLST matrix using GrapeTree"
    }
}
