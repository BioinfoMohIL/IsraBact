version 1.0

task run_diphtOscan {
    input {
        Array[File] assemblies
        Boolean mlst
        Boolean tox
        Boolean res_vir
        Boolean extend_genotyping
        Boolean integron
        Boolean tree
        Int threads
        String docker_image = "bioinfomoh/diphtoscan:latest"
    }

    command {
        diphtoscan \
        -a ~{sep=" " assemblies} \
        ~{if mlst then "-st" else ""} \
        ~{if tox then "-t" else ""} \
        ~{if res_vir then "-res_vir" else ""} \
        ~{if extend_genotyping then "-plus" else ""} \
        ~{if integron then "-integron" else ""} \
        ~{if tree then "-tree" else ""} \
        -o results \
        --threads ~{threads} \
    }

    output {
        Array[File] results = glob("results/*")
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: "8G"
    }
}
