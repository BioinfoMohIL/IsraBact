import pandas as pd

# Charger les fichiers CSV
df_dna = pd.read_csv("/home/user1/Desktop/programs/salmonella/plasmid_vir_dna_seq.csv")
df_virulence = pd.read_csv("plasmid_vir_factors.csv")

# Renommer colonne 'genome_id' → 'ID' pour merger

# S'assurer que ID est de type string
df_dna["id"] = df_dna["id"].astype(str)
df_virulence["id"] = df_virulence["id"].astype(str)

# Fusionner
merged_df = pd.merge(df_dna, df_virulence, on="id", how="left")

# Sauvegarder
merged_df.to_csv("plasmid_all.csv", index=False)

print("✅ CSV fusionné dans merged_salmonella.csv")
