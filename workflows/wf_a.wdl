version 1.0

import '../tasks/tasks_versioning.wdl' as v


workflow sbt_analysis {
    

    call v.version_capture {
        input:
    }

  

    output {
        String test = version_capture.version
        String sbt_elgato_version = elgato_reads.elgato_version
      
    }
}



