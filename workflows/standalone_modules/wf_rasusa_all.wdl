version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/utilities/task_rasusa.wdl" as rasusa
import "../../tasks/quality_control/comparisons/task_screen.wdl" as task_screen

workflow wf_rasusa_workflow {
  input {
    File read1
    String samplename
    File read2
    Float? coverage
    Int? genome_length

    Int min_reads = 7472
    Int min_basepairs = 2241820
    Int min_genome_length = 100000
    Int max_genome_length = 18040666
    Int min_coverage = 10
    Int min_proportion = 40

    Boolean force = false


  }

  call task_screen.check_reads as raw_check_reads{
      input:
        read1 = read1,
        read2 = read2,
        min_reads = min_reads,
        min_basepairs = min_basepairs,
        min_genome_length = min_genome_length,
        max_genome_length = max_genome_length,
        min_coverage = min_coverage,
        min_proportion = min_proportion,
        expected_genome_length = genome_length,
        workflow_series = "theiaprok"
    
  }

  if (raw_check_reads.read_screen != "PASS" || force) {

    call fetch_genome_length {
        input:
            samplename = samplename

    }

    call rasusa.rasusa as rasusa_task {
        input:
            read1 = read1,
            read2 = read2,
            samplename = samplename,
            genome_length = fetch_genome_length.genome_length,
            coverage = coverage
    }

  }

  call versioning.version_capture {
    input:
  }

  output {
    String rasusa_wf_version_IGNORE = version_capture.version
    String rasusa_wf_analysis_date_IGNORE = version_capture.date
    String rasusa_screen_pass_IGNORE = raw_check_reads.read_screen

    String? rasusa_version_IGNORE = rasusa_task.rasusa_version
    File? read1_subsampled = rasusa_task.read1_subsampled
    File? read2_subsampled = rasusa_task.read2_subsampled
  }
}

task fetch_genome_length {
    input {
        String samplename
    }

    command <<<
        declare -A genome_lengths
        genome_lengths["NM"]=2300000
        genome_lengths["M"]=2300000
        genome_lengths["NG"]=2200000
        genome_lengths["HI"]=1800000
        genome_lengths["SH"]=4800000
        genome_lengths["SO"]=4800000
        genome_lengths["LC"]=2900000
        genome_lengths["LF"]=2900000
        genome_lengths["SG"]=4500000
        genome_lengths["CA"]=1700000
        genome_lengths["VIB"]=5000000
        genome_lengths["V"]=5000000
        genome_lengths["EC"]=5000000
        genome_lengths["SA"]=2800000
        genome_lengths["BP"]=4100000
        genome_lengths["SP"]=2100000
        genome_lengths["ST"]=2000000
        genome_lengths["LG"]=3400000
        genome_lengths["LW"]=3400000
        genome_lengths["CB"]=2500000

        
    #    prefix=$(echo "~{samplename}" | sed -E 's/^([A-Z]+)[0-9].*$/\1/')

        for key in "${!genome_lengths[@]}"; do
            if [[ "~{samplename}" == *"$key"* ]]; then
                prefix="$key"
                break
            fi
        done

        if [ -z "$prefix" ]; then
            echo "Unknown prefix in ~{samplename}" >&2
            exit 1
        fi


        genome_length=${genome_lengths[$prefix]}

        echo "Filename: ~{samplename}"
        echo "Prefix: $prefix"
        echo "Genome length: ${genome_length}"

        echo ${genome_length} > genome_length.txt
    >>>

    output {
        String genome_length = read_string("genome_length.txt")
    }

    runtime { 
        docker: "ubuntu:22.04"
        cpu: 4
    }
}