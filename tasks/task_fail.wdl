version 1.0

task fail {
    input {
        String message
    }
    command <<< 
        echo "~{message}"
        exit 1
    >>>

    output {
        String dummy = "Fail triggered"
    }
    
    runtime {
        docker: "alpine" 
    }
}