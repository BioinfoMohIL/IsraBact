version 1.0

task unicycler {
    input {
        File read1
        File? read2
        String samplename

        File? illumina_unpaired
        File? long_reads

        Int disk_size = 100
        Int cpu = 16
        Int memory = 16
        String docker = "quay.io/staphb/unicycler:0.5.1"
    }

    command <<<
        unicycler --version 2>&1 | sed -E 's/.*Unicycler v([0-9.]+).*/\1/' | tee VERSION

        unicycler --threads ~{cpu} \
            -1 ~{read1} -2 ~{read2} \
            ~{'-s ' + illumina_unpaired} \
            ~{'-l ' + long_reads} \
            -o output

        cp "output/assembly.fasta" "~{samplename}.fasta"
        cp "output/assembly.gfa" "~{samplename}.gfa"
    >>>

    output {
        File assembly_fasta = "~{samplename}.fasta"
        File assembly_gfa = "~{samplename}.gfa"
        File unicycler_log = "output/unicycler.log"
        String unicycler_version = read_string("VERSION")
        Array[File] unicycler_intermediate_graphs = glob("output/0*.gfa")
    }

    runtime {
        cpu: "~{cpu}"
        docker: "~{docker}"
        memory: "~{memory} GiB"
        disks: "local-disk " + "~{disk_size}" + " SSD"
        preemptible: 1
        maxRetries:  2
    }
}
