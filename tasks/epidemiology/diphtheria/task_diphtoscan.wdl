version 1.0

task diphtoscan {
    input {
        Array[File] my_input

        Boolean mlst
        Boolean tox
        Boolean res_vir
        Boolean extend_genotyping
        Boolean integron
        Boolean tree

        Int threads 
        Boolean update_db
        String docker_image = docker_image
    }

    command <<<
        # Update the db before
        ~{if update_db then "diphtoscan -u" else ""}

        diphtoscan \
        -a ~{sep=" " my_input} \
        ~{if mlst then "-st" else ""} \
        ~{if tox then "-t" else ""} \
        ~{if res_vir then "-res_vir" else ""} \
        ~{if extend_genotyping then "-plus" else ""} \
        ~{if integron then "-integron" else ""} \
        ~{if tree then "-tree" else ""} \
        -o results \
        --threads ~{threads} 

        if [[ -d results && "$(ls -A results)" ]]; then
            zip -r results.zip results
        fi


    >>>


    output {
        File results = "results.zip"
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: "8G"
    }
}
