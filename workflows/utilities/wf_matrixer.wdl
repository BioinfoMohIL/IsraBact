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
        File matrix_tsv = build_matrix.output_matrix
    }
}

task build_matrix {
    input {
        File amr_tsv
    }

    command <<<
        python3 << 'EOF'
        import pandas as pd

        # Lecture du TSV
        df = pd.read_csv("~{amr_tsv}", sep="\t", dtype=str)

        # Nettoyage des colonnes : minuscules et strip
        df.columns = [c.strip().lower() for c in df.columns]

        # Supprimer les lignes qui sont des headers répétés
        df = df[df['name'].str.lower() != 'name']

        # Colonnes requises
        required_cols = ['name', 'element symbol', '% identity to reference']
        for col in required_cols:
            if col not in df.columns:
                raise ValueError(f"Colonne manquante : {col}")

        # Renommage pour uniformité
        df = df.rename(columns={
            'name': 'Name',
            'element symbol': 'Element',
            '% identity to reference': 'Identity'
        })

        # Construction de la matrice pivot
        matrix = df.pivot_table(
            index='Name',
            columns='Element',
            values='Identity',
            aggfunc='first'  # si doublons, prend le premier
        )

        # Tri des colonnes par ordre alphabétique
        matrix = matrix.reindex(sorted(matrix.columns), axis=1)

        # Remplacer NaN par vide
        matrix = matrix.fillna("0")

        # Sauvegarde finale
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


