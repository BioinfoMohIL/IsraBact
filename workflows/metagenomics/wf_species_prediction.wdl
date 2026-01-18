version 1.0

workflow strain_prediction {
  input {
    File read1
    File read2
    File model_pickle
    Int max_reads = 5000
  }

    call predict_strain {
      input:
        read1 = read1,
        read2 = read2,
        model_pickle = model_pickle,
        max_reads = max_reads
    }
  

  output {
    File all_predictions = predict_strain.prediction_csv
  }

}

task predict_strain {
  input {
    File read1
    File read2
    File model_pickle
    Int max_reads = 5000
  }

  command {
    python3 <<CODE
        import os
        import csv
        import gzip
        from collections import Counter
        import pickle
        import sys

        def read_fastq(file, max_reads=2000):
            reads = []
            with gzip.open(file, "rt") as f:
                for i, line in enumerate(f):
                    if i % 4 == 1:
                        reads.append(line.strip())
                        if len(reads) >= max_reads:
                            break
            return reads

        def read_paired_fastq(file1, file2, max_reads=2000):
            reads = []
            reads.extend(read_fastq(file1, max_reads=max_reads//2))
            reads.extend(read_fastq(file2, max_reads=max_reads//2))
            return reads

        def kmers(seq, k=6):
            return [seq[i:i+k] for i in range(len(seq)-k+1)]

        def extract_features(reads, k=6):
            counter = Counter()
            for seq in reads:
                counter.update(kmers(seq, k))
            return dict(counter)

        # Load model
        with open("${model_pickle}", "rb") as f:
            model, vec = pickle.load(f)

        # Process files
        reads = read_paired_fastq("${read1}", "${fq_R2}", max_reads=${max_reads})
        features = extract_features(reads, k=6)
        X_new = vec.transform([features])

        predicted_strain = model.predict(X_new)
        prob = model.predict_proba(X_new)

        # Save CSV
        output_csv = "prediction_${read1.basename}.csv"
        with open(output_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["sample", "predicted"] + list(model.classes_))
            writer.writeheader()
            writer.writerow({
                "sample": os.path.basename("${read1}").split('_')[0],
                "predicted": predicted_strain[0],
                **dict(zip(model.classes_, prob[0]))
            })
        CODE
  }

  output {
    File prediction_csv = "prediction_${basename(read1)}.csv"
  }

  runtime {
    docker: "python:3.11-slim"
    memory: "4G"
    cpu: 2
  }
}

