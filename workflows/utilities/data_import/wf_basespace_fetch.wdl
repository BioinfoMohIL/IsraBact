version 1.0

import "../../../tasks/task_versioning.wdl" as versioning
import "../../../tasks/utilities/data_import/task_basespace_cli.wdl" as basespace

workflow basespace_fetch {
  meta {
    description: "Fetch reads from Basespace via Basespace CLI."
    author: "David Maimoun (The Codon Bleu)"
    email: "thecodonbleu@outlook.com"
  }

  input {
    String sample_name
    String basespace_sample_name
    String? basespace_sample_id
    String basespace_collection_id
    String api_server
    String access_token
  }

  call basespace.fetch_bs as fetch_bs {
    input:
      sample_name = sample_name,
      basespace_sample_id = basespace_sample_id,
      basespace_sample_name = basespace_sample_name,
      basespace_collection_id = basespace_collection_id,
      api_server = api_server,
      access_token = access_token
  }

  call versioning.version_capture {
    input:
  }

  output {
    String basespace_fetch_version = version_capture.version
    String basespace_fetch_analysis_date = version_capture.date
    File read1 = fetch_bs.read1
    File? read2 = fetch_bs.read2
    Float read1_file_size_mb = fetch_bs.fwd_file_size
    Float read2_file_size_mb = fetch_bs.rev_file_size
    Float bs_read1_file_size_mb = fetch_bs.bs_fwd_file_size
    Float bs_read2_file_size_mb = fetch_bs.bs_rev_file_size
  }
}


