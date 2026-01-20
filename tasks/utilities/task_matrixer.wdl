version 1.0

task build_matrix {
    input {
        File input_tsv
        String scope
    }

    command <<<
        python3 << 'EOF'
        import pandas as pd

        # Lire le TSV
        df = pd.read_csv("~{input_tsv}", sep="\t", dtype=str)

        # Nettoyer les colonnes
        df.columns = [c.strip().lower() for c in df.columns]

        # Supprimer les lignes où "name" apparaît comme ligne d'en-tête
        df = df[df['name'].str.lower() != 'name']

        # Colonnes requises
        required_cols = ['name', 'element symbol', '% identity to reference', 'scope']
        for col in required_cols:
            if col not in df.columns:
                raise ValueError(f"Colonne manquante : {col}")

        df = df.rename(columns={
            'name': 'Name',
            'element symbol': 'Element',
            '% identity to reference': 'Identity',
            'scope': 'Scope'
        })



        if "~{scope}".lower() == "core":
            df = df[df['Scope'].str.lower() == 'core']

        print(df['Name'],df['Scope'])
        # Construire la matrice pivot
        matrix = df.pivot_table(
            index='Name',
            columns='Element',
            values='Identity',
            aggfunc='first'  # if duplicate, pick first
        )

        matrix = matrix.reindex(sorted(matrix.columns), axis=1)
        matrix = matrix.fillna("0")

        # Sauvegarder
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
