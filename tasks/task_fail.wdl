version 1.0

task fail {
    input {
        String message
    }
    command <<< 
        echo "~{message}" | tee "ERROR"
        exit 1
    >>>

    output {
        File fail_logs = 'ERROR'
    }

    runtime {
        docker: "alpine" 
    }
}