version 1.0

task version_capture {
  input {
    String? timezone
    String docker = "debian:12-slim"
  }

  meta {
    volatile: true
  }


  command <<<
    VERSION="IsraBact v1.0"
    ~{default='' 'export TZ=' + timezone}
    date +"%Y-%m-%d" > TODAY
    echo "$VERSION" > VERSION
  >>>

  output {
    String date = read_string("TODAY")
    String version = read_string("VERSION")
  }

  runtime {
    memory: "1 GB"
    cpu: 1
    docker: docker
    disks: "local-disk 10 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2" 
    preemptible: 1
  }
  
}