version 1.0

task build_matrix {
    input {
        File input_tsv
        String scope
    }

    command <<<
        python3 << 'EOF'
        import pandas as pd

        df = pd.read_csv("~{input_tsv}", sep="\t", dtype=str)

        # Clean cols
        df.columns = [c.strip().lower() for c in df.columns]

        # Remove row when "name" (appears because the Concatenate_Columns_Content_PHB concatenation )
        df = df[df['name'].str.lower() != 'name']

        required_cols = ['name', 'element symbol', '% identity to reference', 'scope']
        for col in required_cols:
            if col not in df.columns:
                raise ValueError(f"Column missing : {col}")

        df = df.rename(columns={
            'name': 'Name',
            'element symbol': 'Element',
            '% identity to reference': 'Identity',
            'scope': 'Scope'
        })

        if "~{scope}".lower() == "core":
            df = df[df['Scope'].str.lower() == 'core']

        print(df['Name'],df['Scope'])

        matrix = df.pivot_table(
            index='Name',
            columns='Element',
            values='Identity',
            aggfunc='first'  # if duplicate, pick first
        )

        matrix = matrix.reindex(sorted(matrix.columns), axis=1)
        matrix = matrix.fillna("0")
        matrix.to_csv(f"matrix_~{scope}_genes.tsv", sep="\t")

        EOF
    >>>

    output {
        File output_matrix = "matrix_~{scope}_genes.tsv"
    }

    runtime {
        docker: "bioinfomoh/data_analysis_pytools:1"
        memory: "2G"
        cpu: 1
    }
}
