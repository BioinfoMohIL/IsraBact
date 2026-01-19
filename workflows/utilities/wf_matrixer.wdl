version 1.0

workflow wf_matrixer {
    input {
        File input_tsv
    }

    call build_matrix {
        input:
            amr_tsv = input_tsv
    }

    output {
        File matrixer_matrix = build_matrix.output_matrix
    }
}

task build_matrix {
    input {
        File amr_tsv
    }

    command <<<
        python3 << 'EOF'
        import pandas as pd

        df = pd.read_csv("~{amr_tsv}", sep="\t", dtype=str)

        # Clean cols
        df.columns = [c.strip().lower() for c in df.columns]

        # Del rows with header (new header with each sample appears when concate)
        df = df[df['name'].str.lower() != 'name']

        required_cols = ['name', 'element symbol', '% identity to reference']
        for col in required_cols:
            if col not in df.columns:
                raise ValueError(f"Colonne manquante : {col}")

        df = df.rename(columns={
            'name': 'Name',
            'element symbol': 'Element',
            '% identity to reference': 'Identity'
        })

        matrix = df.pivot_table(
            index='Name',
            columns='Element',
            values='Identity',
            aggfunc='first'  # if duplicate, take first
        )

        matrix = matrix.reindex(sorted(matrix.columns), axis=1)

        matrix = matrix.fillna("0")

        matrix.to_csv("matrix.tsv", sep="\t")
        EOF
    >>>

    output {
        File output_matrix = "matrix.tsv"
    }

    runtime {
        docker: "bioinfomoh/data_analysis_pytools:1"
        memory: "2G"
        cpu: 1
    }
}


