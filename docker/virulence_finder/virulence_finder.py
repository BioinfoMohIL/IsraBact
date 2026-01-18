import os
import json
import argparse
import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from tempfile import NamedTemporaryFile

def color_presence(val):
    if val == "1":
        return 'background-color: lightgreen; text-align: center;'
    elif val == "0":
        return 'background-color: lightcoral; text-align: center;'
    return ''

# Create styled table
def create_html_table(df, samplename, is_plamid=False):
    styled = df.style

    styled = styled\
        .applymap(color_presence, subset=[samplename])\
        .set_table_styles([
            {'selector': 'table', 'props': [('margin-left', 'auto'), ('margin-right', 'auto'), ('border-collapse', 'collapse')]},  # center table
            {'selector': 'tr:nth-child(even)', 'props': [('background-color', '#f9f9f9')]},  # striped even rows
            {'selector': 'tr:hover', 'props': [('background-color', '#e0e0e0')]},           # hover effect
            {'selector': 'th', 'props': [('background-color', '#333'), ('color', 'white'), ('padding', '8px 12px')]},
            {'selector': 'td', 'props': [('padding', '12px 12px')]},

            # Column widths
            {'selector': 'th:nth-child(1), td:nth-child(1)', 'props': [('width', '40%')]},
            {'selector': 'th:nth-child(2), td:nth-child(2)', 'props': [('width', '30%')]},
            {'selector': 'th:nth-child(3), td:nth-child(3)', 'props': [('width', '30%')]},
        ])\
        .hide(axis='index')

    # Render HTML table as string
    html_table = styled.to_html()

    title = f"{samplename} {'Plasmid' if is_plamid else ''}"


    # Compose full HTML with centered title and table
    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <title>{title} Virulence Comparison Results</title>
        <style>
            body {{
                font-family: Arial, sans-serif;
                margin: 40px;
                background-color: #f5f5f5;
                text-align: center;
            }}
            h1 {{
                color: #333;
                margin-bottom: 30px;
            }}
            table {{
                margin-left: auto;
                margin-right: auto;
                border-collapse: collapse;
                box-shadow: 0 2px 5px rgba(0,0,0,0.15);
                background: white;
                width: 70%;  
            }}
            th, td {{
                border: 1px solid #ccc;
                background-color: #f0f0f0;
                padding: 8px 12px;
            }}
            tr:nth-child(even) {{
                background-color: #f9f9f9;
            }}
            tr:hover {{
                background-color: #e0e0e0;
            }}
        </style>
    </head>
    <body>
        <h1>{title} Virulence Comparison Results</h1>
        {html_table}
    </body>
    </html>
    """

    return html_content

def parse_args():
    parser = argparse.ArgumentParser(description="Profile Assessment Tool")
    
    parser.add_argument(
        '--species',
        required=True,
        help='Species name (e.g., Salmonella enterica)'
    )

    parser.add_argument(
        '--fasta',
        required=True,
        help='Path to input FASTA file'
    )

    parser.add_argument(
        '--samplename',
        required=False,
        help='Here the samplename'
    )
    
    parser.add_argument(
        '--plasmid_check',
        action='store_true',
        help='Check the plasmid virulence'
    )

    return parser.parse_args()

def check_virulence(data, hits_csv, summary_csv, summary_html, is_plamid=False):
     # === Load virulence genes JSON and write to FASTA ===
    with open(data) as f:
        virulence_data = json.load(f)

    with NamedTemporaryFile("w+", delete=False, suffix=".fasta") as vfasta:
        for gene in virulence_data:
            seq = gene.get("dna_sequence")
            if seq:
                header = f">{gene['Virulence_Gene']}|{gene.get('Locus_Tag', '')}|{gene.get('Location', '')}"
                vfasta.write(f"{header}\n{seq}\n")
        vfasta_path = vfasta.name

    # === Create BLAST DB ===
    db_name = os.path.join(f"{samplename}_blastdb")
    os.system(f"makeblastdb -in {fasta_file} -dbtype nucl -out {db_name}")

    # === Run BLASTn with megablast ===
    blast_output = os.path.join(f"{samplename}_blast_results.xml")
    blastn_cline = NcbiblastnCommandline(
        query=vfasta_path,
        db=db_name,
        evalue=1e-3,
        outfmt=5,  # XML
        out=blast_output,
        task="megablast"
    )
    blastn_cline()

    # === Parse BLAST XML and collect best hits ===
    match_rows = []
    presence_map = {}

    with open(blast_output) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        for record in blast_records:
            qinfo = record.query.split("|")
            gene = qinfo[0]
            locus = qinfo[1] if len(qinfo) > 1 else ""
            location = qinfo[2] if len(qinfo) > 2 else ""

            best_hsp = None
            best_score = -1

            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect <= 1e-3 and hsp.bits > best_score:
                        best_score = hsp.bits
                        best_hsp = {
                            "query_acc_ver": gene,
                            "subject_acc_ver": alignment.hit_def,
                            "location": location,
                            "% identity": round((hsp.identities / hsp.align_length) * 100, 3),
                            "alignment_length": hsp.align_length,
                            "mismatches": hsp.align_length - hsp.identities,
                            "gap_opens": hsp.gaps,
                            "q_start": hsp.query_start,
                            "q_end": hsp.query_end,
                            "s_start": hsp.sbjct_start,
                            "s_end": hsp.sbjct_end,
                            "evalue": hsp.expect,
                            "bit_score": hsp.bits,
                            "locus": locus,
                            "product": next((g["Product"] for g in virulence_data if g["Virulence_Gene"] == gene), "")
                        }

            if best_hsp:
                match_rows.append(best_hsp)
                presence_map[gene] = {
                    "Factor": gene,
                    "Location": location,
                    samplename: "1"
                }
            else:
                presence_map[gene] = {
                    "Factor": gene,
                    "Location": location,
                    samplename: "0"
                }

    # === Convert to DataFrames ===
    df_presence = pd.DataFrame(presence_map.values())
    df_presence = df_presence.sort_values("Factor", key=lambda col: col.str.lower())
    df_matches = pd.DataFrame(match_rows)
    # Sort so that rows with Location == 'genome' come first, then others alphabetically
    df_presence['is_genome'] = df_presence['Location'].apply(lambda x: 0 if x == 'genome' else 1)
    df_presence = df_presence.sort_values(by=['is_genome', 'Location', 'Factor'], key=lambda col: col.str.lower() if col.dtype == 'object' else col)
    df_presence = df_presence.drop(columns=['is_genome'])

    # Save CSV outputs as before
    df_presence.to_csv(summary_csv, index=False)
    df_matches.to_csv(hits_csv, index=False)


    html_content = create_html_table(df_presence, samplename, is_plamid)

    with open(summary_html, 'w') as f:
        f.write(html_content)


if __name__ == '__main__':
    args = parse_args()
    species = args.species.lower()
    fasta_file = args.fasta

    print("Species:", species)
    print("FASTA file:", fasta_file)

    # === Paths ===
    data = os.path.join("/app", "data", species, f"{species}_all_data.json")

    samplename = args.samplename if args.samplename else os.path.splitext(os.path.basename(fasta_file))[0]

    check_virulence(
        data=data, 
        hits_csv=f"{samplename}_virulence_hits.csv", 
        summary_csv=f"{samplename}_virulence_summary.csv", 
        summary_html=f"{samplename}_virulence_summary.html")
    
    if args.plasmid_check:
        check_virulence(
            data = os.path.join("/app", "data", 'plasmid', f"plasmid_all_data.json"), 
            hits_csv=f"{samplename}_plasm_virulence_hits.csv", 
            summary_csv=f"{samplename}_plasm_virulence_summary.csv", 
            summary_html=f"{samplename}_plasm_virulence_summary.html",
            is_plamid=True)


    print("✅ Virulence profile assessment completed, sorted, and cleaned up.")