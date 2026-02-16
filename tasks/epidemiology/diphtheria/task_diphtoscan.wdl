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
    }

    command <<<
        # Update the db before
        # ~{if update_db then "diphtoscan -u" else ""}
        output_name="results"

        diphtoscan \
        -a ~{sep=" " my_input} \
        ~{if mlst then "-st" else ""} \
        ~{if tox then "-t" else ""} \
        ~{if res_vir then "-res_vir" else ""} \
        ~{if extend_genotyping then "-plus" else ""} \
        ~{if integron then "-integron" else ""} \
        ~{if tree then "-tree" else ""} \
        -o ${output_name} \
        --threads ~{threads} 

        if [[ -d results && "$(ls -A results)" ]]; then
            zip -r results.zip results
        fi

        # Run this block only if tree = true
        ~{if tree then "
            mkdir -p jolytree_results
            echo Gathering tree results ...
            mv *jolytree.* jolytree_results/ 2>/dev/null

            if [[ -d jolytree_results && \"\$(ls -A jolytree_results)\" ]]; then
                zip -r jolytree_results.zip jolytree_results
            fi
        " else ""}
    >>>


    output {
        File results = "results.zip"
        File? tree_results = "jolytree_results.zip"

    }

    runtime {
        docker: "bioinfomoh/diphtoscan:latest"
        cpu: threads
        memory: "8G"
    }
}
