version 1.0

task build_matrix {
    input {
        File amr_tsv
        String genes
    }

    command <<<
        python3 << 'EOF'
        import pandas as pd

        df = pd.read_csv("~{amr_tsv}", sep="\t", dtype=str)

        # Clean column names
        df.columns = [c.strip().lower() for c in df.columns]

        # Remove duplicated headers after concat
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

        # Build base matrix
        matrix = df.pivot_table(
            index='Name',
            columns='Element',
            values='Identity',
            aggfunc='first'
        )

        # If genes list is provided, restrict columns to it
        genes_str = "~{genes}".strip()
        if genes_str and genes_str.lower() != "na":
            genes_list = [g.strip() for g in genes_str.split(",") if g.strip()]

            # Add missing genes as empty columns
            for g in genes_list:
                if g not in matrix.columns:
                    matrix[g] = "0"

            # Reorder columns exactly as provided
            matrix = matrix[genes_list]
        else:
            # Default behavior: sort all genes alphabetically
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
