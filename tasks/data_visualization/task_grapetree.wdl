version 1.0

task grapetree {
    input {
        File matrix

        String tree_method = "MSTreeV2"

        String docker_image = "bioinfomoh/grapetree:1"

        Int cpu = 2
        Int memory_gb = 4
        Int disk_size_gb = 10
    }

    command <<<
        set -euxo pipefail

        BASENAME=$(basename "~{matrix}" .tsv)

        grapetree \
            --method ~{tree_method} \
            --profile ~{matrix} \
            > "${BASENAME}.nwk"

        test -s "${BASENAME}.nwk"
    >>>

    output {
        File tree_newick = glob("*.nwk")[0]
    }

    runtime {
        docker: docker_image
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_size_gb} HDD"
    }
}