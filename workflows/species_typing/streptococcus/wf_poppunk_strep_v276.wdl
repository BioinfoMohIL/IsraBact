version 1.0

import "../../../tasks/species_typing/streptococcus/task_poppunk_streppneumo_v2-7-6.wdl" as st_poppunk

workflow wf_strep_typing_poppunk {
    meta {
        description: "Streptococcus Serotying using popPUNK v2.7.6 ."
        author: "David Maimoun"
        organization: "MOH Jerusalem"
    }

    input {
        File assembly
        String samplename
        File? GPS_db = "gs://fc-5d4556f8-3de6-4709-85da-11445772644d/datasets/sp/GPS/GPS_v10.tar.gz"
        File? GPS_external_cluster = "gs://fc-5d4556f8-3de6-4709-85da-11445772644d/datasets/sp/GPS/GPS_v10_external_clusters.csv"
    }

    call st_poppunk.poppunk  {
        input:
            assembly = assembly,
            samplename = samplename,
            GPS_db = GPS_db,
            GPS_external_cluster = GPS_external_cluster
    }

    output {
        String  poppunk_v276_gps_cluster = poppunk.poppunk_gps_cluster
        String  poppunk_v276_version = poppunk.poppunk_version
        String  poppunk_v276_GPS_db_version = poppunk.poppunk_GPS_db_version
        File?   poppunk_v276_dists_npy = poppunk.poppunk_dists_npy 
        File?   poppunk_v276_dists_pkl = poppunk.poppunk_dists_pkl 
        File?   poppunk_v276_h5 = poppunk.poppunk_h5
        File?   poppunk_v276_gps_external_cluster_csv = poppunk.poppunk_gps_external_cluster_csv
    }

}